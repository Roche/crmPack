#' @include Data-methods.R
#' @include Design-class.R
#' @include McmcOptions-class.R
#' @include Rules-methods.R
#' @include Simulations-class.R
#' @include helpers.R
#' @include mcmc.R
NULL

# nolint start

## ============================================================

##' Simulate outcomes from a CRM design
##'
##' @param object the \code{\linkS4class{Design}} object we want to simulate
##' data from
##' @param nsim the number of simulations (default: 1)
##' @param seed see \code{\link{set_seed}}
##' @param truth a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity. Additional arguments can be supplied
##' in \code{args}.
##' @param args data frame with arguments for the \code{truth} function. The
##' column names correspond to the argument names, the rows to the values of the
##' arguments. The rows are appropriately recycled in the \code{nsim}
##' simulations. In order to produce outcomes from the posterior predictive
##' distribution, e.g, pass an \code{object} that contains the data observed so
##' far, \code{truth} contains the \code{prob} function from the model in
##' \code{object}, and \code{args} contains posterior samples from the model.
##' @param firstSeparate enroll the first patient separately from the rest of
##' the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
##' in this patient.
##' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
##' giving the MCMC options for each evaluation in the trial. By default,
##' the standard options are used
##' @param parallel should the simulation runs be parallelized across the
##' clusters of the computer? (not default)
##' @param nCores how many cores should be used for parallel computing?
##' Defaults to the number of cores on the machine, maximum 5.
##' @param \dots not used
##' @param derive a named list of functions which derives statistics, based on the
##' vector of posterior MTD samples. Each list element must therefore accept
##' one and only one argument, which is a numeric vector, and return a number.
##'
##' @return an object of class \code{\linkS4class{Simulations}}
##'
##' @example examples/design-method-simulate-Design.R
##' @export
##' @importFrom parallel detectCores
##' @keywords methods
setMethod("simulate",
  signature =
    signature(
      object = "Design",
      nsim = "ANY",
      seed = "ANY"
    ),
  def =
    function(object, nsim = 1L, seed = NULL,
             truth, args = NULL, firstSeparate = FALSE,
             mcmcOptions = McmcOptions(),
             parallel = FALSE, nCores =
               min(parallel::detectCores(), 5), derive = list(),
             ...) {
      ## checks and extracts
      assert_function(truth)
      assert_flag(firstSeparate)
      assert_count(nsim, positive = TRUE)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      args <- as.data.frame(args)
      nArgs <- max(nrow(args), 1L)

      ## seed handling
      RNGstate <- set_seed(seed)

      ## from this,
      ## generate the individual seeds for the simulation runs
      simSeeds <- sample.int(n = 2147483647, size = as.integer(nsim))

      ## the function to produce the run a single simulation
      ## with index "iterSim"
      runSim <- function(iterSim) {
        ## set the seed for this run
        set.seed(simSeeds[iterSim])

        ## what is now the argument for the truth?
        ## (appropriately recycled)
        thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

        ## start the simulated data with the provided one
        thisData <- object@data

        # In case there are placebo
        if (thisData@placebo) {
          ## what is the probability for tox. at placebo?
          thisProb.PL <- h_this_truth(
            object@data@doseGrid[1],
            thisArgs,
            truth
          )
        }

        ## shall we stop the trial?
        ## First, we want to continue with the starting dose.
        ## This variable is updated after each cohort in the loop.
        stopit <- FALSE

        ## what is the next dose to be used?
        ## initialize with starting dose
        thisDose <- object@startingDose

        ## inside this loop we simulate the whole trial, until stopping
        while (!stopit) {
          ## what is the probability for tox. at this dose?
          thisProb <- h_this_truth(
            thisDose,
            thisArgs,
            truth
          )

          ## what is the cohort size at this dose?
          thisSize <- size(object@cohort_size,
            dose = thisDose,
            data = thisData
          )

          ## In case there are placebo
          if (thisData@placebo) {
            thisSize.PL <- size(object@pl_cohort_size,
              dose = thisDose,
              data = thisData
            )
          }

          thisData <- h_determine_dlts(
            data = thisData,
            dose = thisDose,
            prob = thisProb,
            prob_placebo = thisProb.PL,
            cohort_size = thisSize,
            cohort_size_placebo = thisSize.PL,
            dose_grid = object@data@doseGrid[1],
            first_separate = firstSeparate
          )

          ## what is the dose limit?
          doselimit <- maxDose(object@increments,
            data = thisData
          )

          ## generate samples from the model
          thisSamples <- mcmc(
            data = thisData,
            model = object@model,
            options = mcmcOptions
          )

          ## => what is the next best dose?
          thisDose <- nextBest(object@nextBest,
            doselimit = doselimit,
            samples = thisSamples,
            model = object@model,
            data = thisData
          )$value


          ## evaluate stopping rules
          stopit <- stopTrial(object@stopping,
            dose = thisDose,
            samples = thisSamples,
            model = object@model,
            data = thisData
          )

          stopit_results <- h_unpack_stopit(stopit)
        }

        ## get the fit
        thisFit <- fit(
          object = thisSamples,
          model = object@model,
          data = thisData
        )

        # Get the MTD estimate from the samples.

        target_dose_samples <- dose(
          mean(object@nextBest@target),
          model = object@model,
          samples = thisSamples
        )

        # Create a function for additional statistical summary.

        additional_stats <- lapply(derive, function(f) f(target_dose_samples))

        ## return the results
        thisResult <-
          list(
            data = thisData,
            dose = thisDose,
            fit =
              subset(thisFit,
                select = c(middle, lower, upper)
              ),
            stop =
              attr(
                stopit,
                "message"
              ),
            report_results = stopit_results,
            additional_stats = additional_stats
          )
        return(thisResult)
      }

      resultList <- get_result_list(
        fun = runSim,
        nsim = nsim,
        vars =
          c(
            "simSeeds",
            "args",
            "nArgs",
            "firstSeparate",
            "truth",
            "object",
            "mcmcOptions"
          ),
        parallel = parallel,
        n_cores = nCores
      )

      # format simulation output
      simulations_output <- h_simulations_output_format(resultList)

      ## return the results in the Simulations class object
      ret <- Simulations(
        data = simulations_output$dataList,
        doses = simulations_output$recommendedDoses,
        fit = simulations_output$fitList,
        stop_report = simulations_output$stop_matrix,
        stop_reasons = simulations_output$stopReasons,
        additional_stats = simulations_output$additional_stats,
        seed = RNGstate
      )

      return(ret)
    }
)




##' Simulate outcomes from a rule-based design
##'
##' @param object the \code{\linkS4class{RuleDesign}} object we want to simulate
##' data from
##' @param nsim the number of simulations (default: 1)
##' @param seed see \code{\link{set_seed}}
##' @param truth a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity. Additional arguments can be supplied
##' in \code{args}.
##' @param args data frame with arguments for the \code{truth} function. The
##' column names correspond to the argument names, the rows to the values of the
##' arguments. The rows are appropriately recycled in the \code{nsim}
##' simulations.
##' @param parallel should the simulation runs be parallelized across the
##' clusters of the computer? (not default)
##' @param nCores how many cores should be used for parallel computing?
##' Defaults to the number of cores on the machine, maximum 5.
##' @param \dots not used
##'
##' @return an object of class \code{\linkS4class{GeneralSimulations}}
##'
##' @example examples/design-method-simulate-RuleDesign.R
##' @export
##' @keywords methods
setMethod("simulate",
  signature =
    signature(
      object = "RuleDesign",
      nsim = "ANY",
      seed = "ANY"
    ),
  def =
    function(object, nsim = 1L, seed = NULL,
             truth, args = NULL,
             parallel = FALSE,
             nCores =
               min(parallel::detectCores(), 5L),
             ...) {
      ## checks and extracts
      assert_function(truth)
      assert_count(nsim, positive = TRUE)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      args <- as.data.frame(args)
      nArgs <- max(nrow(args), 1L)

      ## seed handling
      RNGstate <- set_seed(seed)

      ## from this,
      ## generate the individual seeds for the simulation runs
      simSeeds <- sample(x = seq_len(1e5), size = as.integer(nsim))

      ## the function to produce the run a single simulation
      ## with index "iterSim"
      runSim <- function(iterSim) {
        ## set the seed for this run
        set.seed(simSeeds[iterSim])

        ## what is now the argument for the truth?
        ## (appropriately recycled)
        thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

        ## so this truth is...
        thisTruth <- function(dose) {
          do.call(
            truth,
            ## First argument: the dose
            c(
              dose,
              ## Following arguments
              thisArgs
            )
          )
        }

        ## start the simulated data with the provided one
        thisData <- object@data

        ## shall we stop the trial?
        ## First, we want to continue with the starting dose.
        ## This variable is updated after each cohort in the loop.
        stopit <- FALSE

        ## what is the next dose to be used?
        ## initialize with starting dose
        thisDose <- object@startingDose

        ## inside this loop we simulate the whole trial, until stopping
        while (!stopit) {
          ## what is the probability for tox. at this dose?
          thisProb <- thisTruth(thisDose)

          ## what is the cohort size at this dose?
          thisSize <- size(object@cohort_size,
            dose = thisDose,
            data = thisData
          )

          ## simulate DLTs
          thisDLTs <- rbinom(
            n = thisSize,
            size = 1L,
            prob = thisProb
          )

          ## update the data with this cohort
          thisData <- update(
            object = thisData,
            x = thisDose,
            y = thisDLTs
          )

          ## evaluate the rule
          thisOutcome <- nextBest(object@nextBest,
            data = thisData
          )

          thisDose <- thisOutcome$value
          stopit <- thisOutcome$stopHere
        }

        ## return the results
        thisResult <-
          list(
            data = thisData,
            dose = thisDose
          )

        return(thisResult)
      }

      resultList <- get_result_list(
        fun = runSim,
        nsim = nsim,
        vars =
          c(
            "simSeeds",
            "args",
            "nArgs",
            "truth",
            "object"
          ),
        parallel = parallel,
        n_cores = nCores
      )

      ## put everything in the GeneralSimulations format:

      ## setup the list for the simulated data objects
      dataList <- lapply(resultList, "[[", "data")

      ## the vector of the final dose recommendations
      recommendedDoses <- as.numeric(sapply(resultList, "[[", "dose"))

      ## return the results in the GeneralSimulations class object
      ret <- GeneralSimulations(
        data = dataList,
        doses = recommendedDoses,
        seed = RNGstate
      )

      return(ret)
    }
)


##' Simulate outcomes from a dual-endpoint design
##'
##' @param object the \code{\linkS4class{DualDesign}} object we want to simulate
##' data from
##' @param nsim the number of simulations (default: 1)
##' @param seed see \code{\link{set_seed}}
##' @param trueTox a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity. Additional arguments can be supplied
##' in \code{args}.
##' @param trueBiomarker a function which takes as input a dose (vector) and
##' returns the true biomarker level (vector). Additional arguments can be
##' supplied in \code{args}.
##' @param args data frame with arguments for the \code{trueTox} and
##' \code{trueBiomarker} function. The column names correspond to the argument
##' names, the rows to the values of the arguments. The rows are appropriately
##' recycled in the \code{nsim} simulations.
##' @param sigma2W variance for the biomarker measurements
##' @param rho correlation between toxicity and biomarker measurements (default:
##' 0)
##' @param firstSeparate enroll the first patient separately from the rest of
##' the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
##' in this patient.
##' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
##' giving the MCMC options for each evaluation in the trial. By default,
##' the standard options are used
##' @param parallel should the simulation runs be parallelized across the
##' clusters of the computer? (not default)
##' @param nCores how many cores should be used for parallel computing?
##' Defaults to the number of cores on the machine, maximum 5.
##' @param \dots not used
##' @param derive a named list of functions which derives statistics, based on the
##' vector of posterior MTD samples. Each list element must therefore accept
##' one and only one argument, which is a numeric vector, and return a number.
##'
##' @return an object of class \code{\linkS4class{DualSimulations}}
##'
##' @example examples/design-method-simulate-DualDesign.R
##' @importFrom mvtnorm rmvnorm
##' @export
##' @keywords methods
setMethod("simulate",
  signature =
    signature(object = "DualDesign"),
  def =
    function(object, nsim = 1L, seed = NULL,
             trueTox, trueBiomarker, args = NULL,
             sigma2W, rho = 0,
             firstSeparate = FALSE,
             mcmcOptions = McmcOptions(),
             parallel = FALSE,
             nCores =
               min(parallel::detectCores(), 5), derive = list(),
             ...) {
      ## checks and extracts
      assert_function(trueTox)
      assert_function(trueBiomarker)
      assert_number(sigma2W, lower = 0)
      assert_number(rho, lower = -1, upper = 1)
      assert_flag(firstSeparate)
      assert_count(nsim, positive = TRUE)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      args <- as.data.frame(args)
      nArgs <- max(nrow(args), 1L)

      ## get names of arguments (excluding the first one which is the dose)
      trueToxArgnames <- names(formals(trueTox))[-1]
      trueBiomarkerArgnames <- names(formals(trueBiomarker))[-1]

      ## this is the covariance matrix we assume:
      trueCov <- matrix(
        c(
          sigma2W, sqrt(sigma2W) * rho,
          sqrt(sigma2W) * rho, 1
        ),
        nrow = 2, byrow = TRUE
      )

      ## seed handling
      RNGstate <- set_seed(seed)

      ## from this,
      ## generate the individual seeds for the simulation runs
      simSeeds <- sample(x = seq_len(1e5), size = as.integer(nsim))

      ## the function to produce the run a single simulation
      ## with index "iterSim"
      runSim <- function(iterSim) {
        ## set the seed for this run
        set.seed(simSeeds[iterSim])

        ## what is now the argument for the true functions?
        ## (appropriately recycled)
        thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

        ## so the true tox function is:
        thisTrueTox <- function(dose) {
          do.call(
            trueTox,
            ## First argument: the dose
            c(
              dose,
              ## Following arguments: take only those that
              ## are required by the tox function
              as.list(thisArgs)[trueToxArgnames]
            )
          )
        }

        ## and the true biomarker function is:
        thisTrueBiomarker <- function(dose) {
          do.call(
            trueBiomarker,
            ## First argument: the dose
            c(
              dose,
              ## Following arguments: take only those that
              ## are required by the biomarker function
              as.list(thisArgs)[trueBiomarkerArgnames]
            )
          )
        }

        ## start the simulated data with the provided one
        thisData <- object@data

        ## shall we stop the trial?
        ## First, we want to continue with the starting dose.
        ## This variable is updated after each cohort in the loop.
        stopit <- FALSE

        ## what is the next dose to be used?
        ## initialize with starting dose
        thisDose <- object@startingDose

        if (thisData@placebo) {
          ## what is the probability for tox. at placebo?
          thisProb.PL <- thisTrueTox(object@data@doseGrid[1])
          thisMeanZ.PL <- qlogis(thisProb.PL)

          ## what is the biomarker mean at placebo?
          thisMeanBiomarker.PL <- thisTrueBiomarker(object@data@doseGrid[1])
        }

        # In case there are placebo, extract true Toxicity and Efficacy for placebo

        ## inside this loop we simulate the whole trial, until stopping
        while (!stopit) {
          ## what is the probability for tox. at this dose?
          thisProb <- thisTrueTox(thisDose)
          ## and the transformation to the z scale is:
          thisMeanZ <- qlogis(thisProb)

          ## what is the biomarker mean at this dose?
          thisMeanBiomarker <- thisTrueBiomarker(thisDose)

          ## what is the cohort size at this dose?
          thisSize <- size(object@cohort_size,
            dose = thisDose,
            data = thisData
          )

          ## In case there are placebo
          ## what is the cohort size at this dose for Placebo?
          if (thisData@placebo) {
            thisSize.PL <- size(object@pl_cohort_size,
              dose = thisDose,
              data = thisData
            )
          }

          ## simulate tox and biomarker response: depends on whether we
          ## separate the first patient or not.
          tmp <-
            if (firstSeparate && (thisSize > 1L)) {
              ## dose the first patient
              tmpStart <- mvtnorm::rmvnorm(
                n = 1,
                mean =
                  c(
                    thisMeanBiomarker,
                    thisMeanZ
                  ),
                sigma = trueCov
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                tmpStart.PL <- mvtnorm::rmvnorm(
                  n = 1,
                  mean =
                    c(
                      thisMeanBiomarker.PL,
                      thisMeanZ.PL
                    ),
                  sigma = trueCov
                )
              }


              ## if there is no DLT:
              if (tmpStart[, 2] < 0) {
                ## enroll the remaining patients
                tmpStart <-
                  rbind(
                    tmpStart,
                    mvtnorm::rmvnorm(
                      n = thisSize - 1,
                      mean =
                        c(
                          thisMeanBiomarker,
                          thisMeanZ
                        ),
                      sigma = trueCov
                    )
                  )

                if (thisData@placebo && (thisSize.PL > 0L)) {
                  tmpStart.PL <-
                    rbind(
                      tmpStart.PL,
                      mvtnorm::rmvnorm(
                        n = thisSize.PL,
                        mean =
                          c(
                            thisMeanBiomarker.PL,
                            thisMeanZ.PL
                          ),
                        sigma = trueCov
                      )
                    )
                }
              }

              if (thisData@placebo && (thisSize.PL > 0L)) {
                list(tmpStart = tmpStart, tmpStart.PL = tmpStart.PL)
              } else {
                list(tmpStart = tmpStart)
              }
            } else {
              ## we can directly dose all patients
              tmpStart <- mvtnorm::rmvnorm(
                n = thisSize,
                mean =
                  c(
                    thisMeanBiomarker,
                    thisMeanZ
                  ),
                sigma = trueCov
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                tmpStart.PL <- mvtnorm::rmvnorm(
                  n = thisSize.PL,
                  mean =
                    c(
                      thisMeanBiomarker.PL,
                      thisMeanZ.PL
                    ),
                  sigma = trueCov
                )
              }

              if (thisData@placebo && (thisSize.PL > 0L)) {
                list(tmpStart = tmpStart, tmpStart.PL = tmpStart.PL)
              } else {
                list(tmpStart = tmpStart)
              }
            }

          ## extract biomarker and DLT samples
          thisBiomarkers <- tmp$tmpStart[, 1]
          thisDLTs <- as.integer(tmp$tmpStart[, 2] > 0)

          # in case there are placebo
          if (thisData@placebo && (thisSize.PL > 0L)) {
            thisBiomarkers.PL <- tmp$tmpStart.PL[, 1]
            thisDLTs.PL <- as.integer(tmp$tmpStart.PL[, 2] > 0)

            ## update the data first with placebo...
            thisData <- update(
              object = thisData,
              x = object@data@doseGrid[1],
              y = thisDLTs.PL,
              w = thisBiomarkers.PL,
              check = FALSE
            )

            ### ... and then with active dose
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs,
              w = thisBiomarkers,
              new_cohort = FALSE
            )
          } else {
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs,
              w = thisBiomarkers
            )
          }


          ## what is the dose limit?
          doselimit <- maxDose(object@increments,
            data = thisData
          )

          ## generate samples from the model
          thisSamples <- mcmc(
            data = thisData,
            model = object@model,
            options = mcmcOptions
          )

          ## => what is the next best dose?
          thisDose <- nextBest(object@nextBest,
            doselimit = doselimit,
            samples = thisSamples,
            model = object@model,
            data = thisData
          )$value

          ## evaluate stopping rules
          stopit <- stopTrial(object@stopping,
            dose = thisDose,
            samples = thisSamples,
            model = object@model,
            data = thisData
          )
          stopit_results <- h_unpack_stopit(stopit)
        }

        ## get the fit
        thisFit <- fit(
          object = thisSamples,
          model = object@model,
          data = thisData
        )

        # Get the MTD estimate from the samples.

        target_dose_samples <- dose(
          mean(object@nextBest@target),
          model = object@model,
          samples = thisSamples
        )

        # Create a function for additional statistical summary.

        additional_stats <- lapply(derive, function(f) f(target_dose_samples))


        ## return the results
        thisResult <-
          list(
            data = thisData,
            dose = thisDose,
            fitTox =
              subset(thisFit,
                select =
                  c(middle, lower, upper)
              ),
            fit_biomarker =
              subset(thisFit,
                select =
                  c(
                    middleBiomarker, lowerBiomarker,
                    upperBiomarker
                  )
              ),
            rho_est = median(thisSamples@data$rho),
            sigma2w_est = median(1 / thisSamples@data$precW),
            stop =
              attr(
                stopit,
                "message"
              ),
            additional_stats = additional_stats,
            report_results = stopit_results
          )

        return(thisResult)
      }

      resultList <- get_result_list(
        fun = runSim,
        nsim = nsim,
        vars =
          c(
            "simSeeds",
            "args",
            "nArgs",
            "firstSeparate",
            "trueTox",
            "trueBiomarker",
            "trueCov",
            "object",
            "mcmcOptions"
          ),
        parallel = parallel,
        n_cores = nCores
      )

      ## put everything in the Simulations format:

      ## setup the list for the simulated data objects
      dataList <- lapply(resultList, "[[", "data")

      ## the vector of the final dose recommendations
      recommendedDoses <- as.numeric(sapply(resultList, "[[", "dose"))

      ## vector of rho estimates
      rhoEstimates <- as.numeric(sapply(resultList, "[[", "rho_est"))

      ## vector of sigma2W estimates
      sigma2Westimates <- as.numeric(sapply(resultList, "[[", "sigma2w_est"))

      ## setup the list for the final tox fits
      fitToxList <- lapply(resultList, "[[", "fitTox")

      ## setup the list for the final biomarker fits
      fitBiomarkerList <- lapply(resultList, "[[", "fit_biomarker")

      ## the reasons for stopping
      stopReasons <- lapply(resultList, "[[", "stop")

      # individual stopping rule results as matrix, labels as column names
      stop_results <- lapply(resultList, "[[", "report_results")
      stop_report <- as.matrix(do.call(rbind, stop_results))

      ## For dual simulations summary of additional statistics.
      additional_stats <- lapply(resultList, "[[", "additional_stats")

      ## return the results in the DualSimulations class object
      ret <- DualSimulations(
        data = dataList,
        doses = recommendedDoses,
        rho_est = rhoEstimates,
        sigma2w_est = sigma2Westimates,
        fit = fitToxList,
        fit_biomarker = fitBiomarkerList,
        stop_report = stop_report,
        stop_reasons = stopReasons,
        additional_stats = additional_stats,
        seed = RNGstate
      )

      return(ret)
    }
)


## ============================================================

##' Obtain hypothetical trial course table for a design
##'
##' This generic function takes a design and generates a dataframe
##' showing the beginning of several hypothetical trial courses under
##' the design. This means, from the generated dataframe one can read off:
##' - how many cohorts are required in the optimal case (no DLTs observed) in
##'   order to reach the highest dose of the specified dose grid (or until
##'   the stopping rule is fulfilled)
##' - assuming no DLTs are observed until a certain dose level, what the next
##'   recommended dose is for all possible number of DLTs observed
##' - the actual relative increments that will be used in these cases
##' - whether the trial would stop at a certain cohort
##' Examining the "single trial" behavior of a dose escalation design is
##' the first important step in evaluating a design, and cannot be replaced by
##' studying solely the operating characteristics in "many trials". The cohort
##' sizes are also taken from the design, assuming no DLTs occur until the dose
##' listed.
##'
##' @param object the design (\code{\linkS4class{Design}} or
##' \code{\linkS4class{RuleDesign}} object) we want to examine
##' @param \dots additional arguments (see methods)
##' @param maxNoIncrement maximum number of contiguous next doses at 0
##' DLTs that are the same as before, i.e. no increment (default to 100)
##'
##' @return The data frame
##'
##' @export
##' @keywords methods regression
setGeneric("examine",
  def =
    function(object, ..., maxNoIncrement = 100L) {
      ## check maxNoIncrement argument
      assert_count(maxNoIncrement, positive = TRUE)

      ## there should be no default method,
      ## therefore just forward to next method!
      standardGeneric("examine")
    },
  valueClass = "data.frame"
)


##' @describeIn examine Examine a model-based CRM
##'
##' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
##' giving the MCMC options for each evaluation in the trial. By default,
##' the standard options are used
##'
##' @example examples/design-method-examine-Design.R
setMethod("examine",
  signature =
    signature(object = "Design"),
  def =
    function(object,
             mcmcOptions = McmcOptions(),
             ...,
             maxNoIncrement) {
      ## start with the empty table
      ret <- data.frame(
        dose = numeric(),
        DLTs = integer(),
        nextDose = numeric(),
        stop = logical(),
        increment = integer()
      )

      ## start the base data with the provided one
      baseData <- object@data

      ## are we finished and can stop?
      stopit <- FALSE

      ## counter how many contiguous doses at 0 DLTs with
      ## no increment
      noIncrementCounter <- 0L

      ## what is the next dose to be used?
      ## initialize with starting dose
      thisDose <- object@startingDose

      ## inside this loop we continue filling up the table, until
      ## stopping
      while (!stopit) {
        ## what is the cohort size at this dose?
        thisSize <- size(object@cohort_size,
          dose = thisDose,
          data = baseData
        )

        if (baseData@placebo) {
          thisSize.PL <- size(object@pl_cohort_size,
            dose = thisDose,
            data = baseData
          )
        }

        ## for all possible number of DLTs:
        for (numDLTs in 0:thisSize)
        {
          ## update data with corresponding DLT vector
          if (baseData@placebo && (thisSize.PL > 0L)) {
            thisData <- update(
              object = baseData,
              x = baseData@doseGrid[1],
              y = rep(0, thisSize.PL),
              check = FALSE
            )

            thisData <-
              update(
                object = thisData,
                x = thisDose,
                y =
                  rep(
                    x = c(0, 1),
                    times =
                      c(
                        thisSize - numDLTs,
                        numDLTs
                      )
                  ),
                new_cohort = FALSE
              )
          } else {
            thisData <-
              update(
                object = baseData,
                x = thisDose,
                y =
                  rep(
                    x = c(0, 1),
                    times =
                      c(
                        thisSize - numDLTs,
                        numDLTs
                      )
                  )
              )
          }

          ## what is the dose limit?
          doselimit <- maxDose(object@increments,
            data = thisData
          )

          ## generate samples from the model
          thisSamples <- mcmc(
            data = thisData,
            model = object@model,
            options = mcmcOptions
          )

          ## => what is the next best dose?
          nextDose <- nextBest(object@nextBest,
            doselimit = doselimit,
            samples = thisSamples,
            model = object@model,
            data = thisData
          )$value

          ## compute relative increment in percent
          thisIncrement <-
            round((nextDose - thisDose) / thisDose * 100)

          ## evaluate stopping rules
          stopThisTrial <- stopTrial(object@stopping,
            dose = nextDose,
            samples = thisSamples,
            model = object@model,
            data = thisData
          )

          ## append information to the data frame
          ret <- rbind(
            ret,
            list(
              dose = thisDose,
              DLTs = numDLTs,
              nextDose = nextDose,
              stop = stopThisTrial,
              increment = as.integer(thisIncrement)
            )
          )
        }

        ## change base data
        if (baseData@placebo && (thisSize.PL > 0L)) {
          baseData <-
            update(
              object = baseData,
              x = baseData@doseGrid[1],
              y = rep(0, thisSize.PL),
              check = FALSE
            )

          baseData <-
            update(
              object = baseData,
              x = thisDose,
              y = rep(0, thisSize),
              new_cohort = FALSE
            )
        } else {
          baseData <-
            update(
              object = baseData,
              x = thisDose,
              y = rep(0, thisSize)
            )
        }

        ## what are the results if 0 DLTs?
        resultsNoDLTs <- subset(
          tail(ret, thisSize + 1),
          dose == thisDose & DLTs == 0
        )

        ## what is the new dose then accordingly?
        newDose <- as.numeric(resultsNoDLTs$nextDose)

        ## what is the difference to the previous dose?
        doseDiff <- newDose - thisDose

        ## update the counter for no increments of the dose
        if (doseDiff == 0) {
          noIncrementCounter <- noIncrementCounter + 1L
        } else {
          noIncrementCounter <- 0L
        }

        ## would stopping rule be fulfilled already?
        stopAlready <- resultsNoDLTs$stop

        ## update dose
        thisDose <- newDose

        ## too many times no increment?
        stopNoIncrement <- (noIncrementCounter >= maxNoIncrement)
        if (stopNoIncrement) {
          warning(paste(
            "Stopping because",
            noIncrementCounter,
            "times no increment vs. previous dose"
          ))
        }

        ## check if we can stop:
        ## either when we have reached the highest dose in the
        ## next cohort, or when the stopping rule is already
        ## fulfilled, or when too many times no increment
        stopit <- (thisDose >= max(object@data@doseGrid)) ||
          stopAlready || stopNoIncrement
      }

      return(ret)
    }
)



##' @describeIn examine Examine a rule-based design
##' @example examples/design-method-examine-RuleDesign.R
setMethod("examine",
  signature =
    signature(object = "RuleDesign"),
  def =
    function(object,
             ...,
             maxNoIncrement) {
      ## start with the empty table
      ret <- data.frame(
        dose = numeric(),
        DLTs = integer(),
        nextDose = numeric(),
        stop = logical(),
        increment = integer()
      )

      ## start the base data with the provided one
      baseData <- object@data

      ## are we finished and can stop?
      stopit <- FALSE

      ## counter how many contiguous doses at 0 DLTs with
      ## no increment
      noIncrementCounter <- 0L

      ## what is the next dose to be used?
      ## initialize with starting dose
      thisDose <- object@startingDose

      ## inside this loop we continue filling up the table, until
      ## stopping
      while (!stopit) {
        ## what is the cohort size at this dose?
        thisSize <- size(object@cohort_size,
          dose = thisDose,
          data = baseData
        )

        ## for all possible number of DLTs:
        for (numDLTs in 0:thisSize)
        {
          ## update data with corresponding DLT vector
          thisData <-
            update(
              object = baseData,
              x = thisDose,
              y =
                rep(
                  x = c(0, 1),
                  times =
                    c(
                      thisSize - numDLTs,
                      numDLTs
                    )
                )
            )

          ## evaluate the rule
          thisOutcome <- nextBest(object@nextBest,
            data = thisData
          )

          ## next dose
          nextDose <- thisOutcome$value

          ## do we stop here?
          stopThisTrial <- thisOutcome$stopHere

          ## compute relative increment in percent
          thisIncrement <-
            round((nextDose - thisDose) / thisDose * 100)

          ## append information to the data frame
          ret <- rbind(
            ret,
            list(
              dose = thisDose,
              DLTs = numDLTs,
              nextDose = nextDose,
              stop = stopThisTrial,
              increment = as.integer(thisIncrement)
            )
          )
        }

        ## change base data
        baseData <-
          update(
            object = baseData,
            x = thisDose,
            y = rep(0, thisSize)
          )

        ## what are the results if 0 DLTs?
        resultsNoDLTs <- subset(
          tail(ret, thisSize + 1),
          dose == thisDose & DLTs == 0
        )

        ## what is the new dose then accordingly?
        newDose <- as.numeric(resultsNoDLTs$nextDose)

        ## what is the difference to the previous dose?
        doseDiff <- newDose - thisDose

        ## update the counter for no increments of the dose
        if (doseDiff == 0) {
          noIncrementCounter <- noIncrementCounter + 1L
        } else {
          noIncrementCounter <- 0L
        }

        ## would stopping rule be fulfilled already?
        stopAlready <- resultsNoDLTs$stop

        ## update dose
        thisDose <- newDose

        ## too many times no increment?
        stopNoIncrement <- (noIncrementCounter >= maxNoIncrement)
        if (stopNoIncrement) {
          warning(paste(
            "Stopping because",
            noIncrementCounter,
            "times no increment vs. previous dose"
          ))
        }

        ## check if we can stop:
        ## either when we have reached the highest dose in the
        ## next cohort, or when the stopping rule is already
        ## fulfilled, or when too many times no increment
        stopit <- (thisDose >= max(object@data@doseGrid)) ||
          stopAlready || stopNoIncrement
      }

      return(ret)
    }
)

##' @describeIn examine Examine a model-based CRM
##'
##' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
##' giving the MCMC options for each evaluation in the trial. By default,
##' the standard options are used
##'
##' @example examples/design-method-examine-DADesign.R
setMethod("examine",
  signature =
    signature(object = "DADesign"),
  def =
    function(object, mcmcOptions = McmcOptions(), ...,
             maxNoIncrement) {
      # A function to return follow up fulfull yes (TRUE) vs no (FALSE);
      ready_to_open <- function(day, window, thisSurv) {
        size <- length(thisSurv)
        # the date the patient starts;
        start_time <- apply(rbind(thisSurv[-size], window$patientGap[-1]), 2, min)
        # the relative time for each patient on the specified "date";
        individule_check <- day - cumsum(c(0, start_time))
        # the minial number should be 0;
        individule_check[individule_check < 0] <- 0
        follow_up <- apply(rbind(thisSurv, individule_check), 2, min)
        return(all((follow_up - apply(rbind(window$patientFollow, thisSurv), 2, min)) >= 0) & (max(follow_up) >= min(window$patientFollowMin, max(thisSurv))))
      }

      ## assume we have surfficient patients, i.e. patient can be immediately enrolled
      ## once the trial accumulation is open. This function will tell you when to open
      ## the next cohort;
      # this function applys to all trials;
      nextOpen <- function(window, thisSurv) {
        size <- length(thisSurv)

        window$patientGap <- window$patientGap[1:size] ## if length(window$pt)>length(thisSurv), assume the first length(thisSurv) patients were enrolled;
        ## if the DLT happens before the end of DLT window, then the next
        ## cohort/enrollment of the next patient would happened earlier;
        start_time <- apply(rbind(thisSurv[-size], window$patientGap[-1]), 2, min)
        # duration of the cohort (all DLT windows finished);
        maxT <- max(thisSurv + cumsum(c(0, start_time)))

        meetrequire <- sapply(1:maxT, function(i) {
          ready_to_open(i, window, thisSurv)
        })
        if (sum(meetrequire) > 0) {
          # the earliest time that the require is met;
          time <- min(c(1:maxT)[meetrequire])
        } else {
          time <- maxT
        }

        return(time)
      }

      ## start with the empty table
      ret <- data.frame(
        DLTsearly_1 = integer(), ## JZ: add a cohort index;
        dose = numeric(),
        DLTs = integer(),
        nextDose = numeric(),
        stop = logical(),
        increment = integer()
      )

      ## start the base data with the provided one
      baseData <- object@data

      ## are we finished and can stop?
      stopit <- FALSE

      ## what is the next dose to be used?
      ## initialize with starting dose
      thisDose <- object@startingDose

      ## initial {fact} variables;
      factDLTs <- baseData@y
      factSurv <- baseData@u
      factT0 <- baseData@t0

      ## Initiate "trialtime" which is zero. This is the global time for studies;
      trialtime <- 0

      ## when the current cohort open?
      pretime <- 0

      ## the duration of DLT window
      Tmax <- baseData@Tmax

      ## number of patients with un-completed DLT window;
      ## assume no patient is under DLT observation period at the beginning;
      preSize <- 0

      ## inside this loop we continue filling up the table, until
      ## stopping
      while (!stopit) {
        ## what is the cohort size at this dose?
        thisSize <- size(object@cohort_size,
          dose = thisDose,
          data = baseData
        )

        ## what's the safetywindow
        thisSafetywindow <- windowLength(object@safetyWindow, thisSize)


        # initial parameters
        thisT0 <- trialtime + cumsum(thisSafetywindow$patientGap)

        factDLTs <- c(factDLTs, rep(0, thisSize))

        factSurv <- c(factSurv, rep(Tmax, thisSize))

        factT0 <- c(factT0, thisT0)

        ## The time that the next cohort open
        trialtime <- trialtime + nextOpen(
          window = thisSafetywindow,
          thisSurv = rep(Tmax, thisSize)
        )

        ## In the DA-CRM, we should count the number of patients who is still within the DLT window;
        ## Thus the loop for numDLTs should be 0:nFollow;
        nFollow <- thisSize + preSize

        ## Identify the censored patients;
        ## "thiscensored" will be used in the cases that numDLTs>0;
        npt <- length(baseData@x) # total number of patients

        thiscensored <- c(c(1:npt)[(trialtime - baseData@t0) < baseData@Tmax & baseData@y == 0], (npt + 1):(npt + thisSize))



        ## for all possible number of DLTs:
        for (numDLTs in 0:nFollow)
        {
          ## If numDLTs>0, two extreme cases will be examinated;
          ## (1) DLTs occur on patients with the longer follow ups;
          ## (2) DLTs occur on patients with the shorter follow ups;






          if (numDLTs == 0) {
            baseData <- update(
              object = baseData,
              y = factDLTs, #### the x will be constantly updated according to u
              u = factSurv,
              t0 = factT0,
              x = thisDose,
              trialtime = trialtime
            ) #### the u will be constantly updated



            ## what is the dose limit?
            doselimit <- maxDose(object@increments,
              data = baseData
            )

            ## generate samples from the model
            thisSamples <- mcmc(
              data = baseData,
              model = object@model,
              options = mcmcOptions
            )

            ## => what is the next best dose?
            nextDose <- nextBest(object@nextBest,
              doselimit = doselimit,
              samples = thisSamples,
              model = object@model,
              data = baseData
            )$value

            # ##remove savePlot
            #
            #                                                   savePlot(plot(baseData),name=paste("Dose",thisDose,0,"DLT",nextDose,sep="_"))
            #
            ## compute relative increment in percent
            thisIncrement <-
              round((nextDose - thisDose) / thisDose * 100)

            ## evaluate stopping rules
            stopThisTrial <- stopTrial(object@stopping,
              dose = nextDose,
              samples = thisSamples,
              model = object@model,
              data = baseData
            )

            ## append information to the data frame
            ret <- rbind(
              ret,
              list(
                DLTsearly_1 = 0,
                dose = thisDose,
                DLTs = numDLTs,
                nextDose = nextDose,
                stop = stopThisTrial,
                increment = as.integer(thisIncrement)
              )
            )
            ### comment here to show only no DLTs;
            #                                            }
          } else {
            for (DLTsearly in 1:numDLTs) {
              # Update current {fact} variables
              thisDLTs <- factDLTs
              thisSurv <- factSurv

              if (DLTsearly == 1) {
                # scenario 1: The patients with longest follow up have DLTs

                thisDLTs[thiscensored][1:numDLTs] <- rep(1, rep(numDLTs))

                thisSurv[thiscensored][1:numDLTs] <- apply(rbind(rep(Tmax, numDLTs), c(trialtime - factT0[thiscensored][1:numDLTs])), 2, min)


                thisData <- update(
                  object = baseData,
                  y = thisDLTs, #### the y will be updated according to u
                  u = thisSurv,
                  t0 = factT0,
                  x = thisDose,
                  trialtime = trialtime
                ) #### the u will be updated
              } else {
                # scenario 2: The patients with shortest follow up have DLTs

                thisDLTs[rev(thiscensored)][1:numDLTs] <- rep(1, rep(numDLTs))

                thisSurv[rev(thiscensored)][1:numDLTs] <- c(apply(rbind(rep(1, numDLTs), pretime + 1 - factT0[rev(thiscensored)][1:numDLTs]), 2, max))

                if (numDLTs >= thisSize) {
                  thistime <- 1 + max(thisT0)
                } else {
                  thistime <- trialtime
                }

                thisData <- update(
                  object = baseData,
                  y = thisDLTs, #### the y will be updated according to u
                  u = thisSurv,
                  t0 = factT0,
                  x = thisDose,
                  trialtime = thistime
                ) #### the u will be updated
              }


              ## what is the dose limit?
              doselimit <- maxDose(object@increments,
                data = thisData
              )

              ## generate samples from the model
              thisSamples <- mcmc(
                data = thisData,
                model = object@model,
                options = mcmcOptions
              )

              ## => what is the next best dose?
              nextDose <- nextBest(object@nextBest,
                doselimit = doselimit,
                samples = thisSamples,
                model = object@model,
                data = thisData
              )$value

              # ##remove savePlot
              #                                                   savePlot(plot(thisData),name=paste("Dose",thisDose,numDLTs,"DLT",DLTsearly,nextDose,sep="_"))
              #
              ## compute relative increment in percent
              thisIncrement <-
                round((nextDose - thisDose) / thisDose * 100)

              ## evaluate stopping rules
              stopThisTrial <- stopTrial(object@stopping,
                dose = nextDose,
                samples = thisSamples,
                model = object@model,
                data = thisData
              )

              ## append information to the data frame
              ret <- rbind(
                ret,
                list(
                  DLTsearly_1 = DLTsearly,
                  dose = thisDose,
                  DLTs = numDLTs,
                  nextDose = nextDose,
                  stop = stopThisTrial,
                  increment = as.integer(thisIncrement)
                )
              )
            }
          }
        }

        ## update pretime
        pretime <- trialtime

        ## what are the results if 0 DLTs?
        resultsNoDLTs <- subset(
          ret,
          dose == thisDose & DLTs == 0
        )

        ## what is the new dose according to table?
        newDose <- as.numeric(resultsNoDLTs$nextDose)

        ## what is the difference to the previous dose?
        doseDiff <- newDose - thisDose

        ## would stopping rule be fulfilled already?
        stopAlready <- any(resultsNoDLTs$stop)

        ## update dose
        thisDose <- max(newDose)

        ## number of patients with un-completed DLT window;
        preSize <- sum(baseData@u[baseData@y == 0] < baseData@Tmax)

        ## update the counter for no increments of the dose
        if (all(doseDiff == 0)) {
          noIncrementCounter <- noIncrementCounter + 1L
        } else {
          noIncrementCounter <- 0L
        }

        ## too many times no increment?
        stopNoIncrement <- (noIncrementCounter >= maxNoIncrement)
        if (stopNoIncrement) {
          warning(paste(
            "Stopping because",
            noIncrementCounter,
            "times no increment vs. previous dose"
          ))
        }

        ## check if we can stop:
        ## either when we have reached the highest dose in the
        ## next cohort, or when the stopping rule is already
        ## fulfilled, or when too many times no increment
        stopit <- (thisDose >= max(object@data@doseGrid)) ||
          stopAlready || stopNoIncrement
      }

      return(ret)
    }
)

## ===================================================================================
## ----------------------------------------------------------------------------------------
##  Simulate design using DLE responses only with DLE samples (pseudo DLE model)
## ------------------------------------------------------------------------------------
##' This is a methods to simulate dose escalation procedure only using the DLE responses.
##' This is a method based on the \code{\linkS4class{TDsamplesDesign}} where model used are of
##' \code{\linkS4class{ModelTox}} class object DLE samples are also used
##'
##' @param object the \code{\linkS4class{TDsamplesDesign}} object we want to simulate the data from
##' @param nsim the number of simulations (default :1)
##' @param seed see \code{\link{set_seed}}
##' @param truth a function which takes as input a dose (vector) and returns the true probability
##' (vector) of the occurrence of a DLE. Additional arguments can be supplied in \code{args}.
##' @param args data frame with arguments for the \code{truth} function. The
##' column names correspond to the argument names, the rows to the values of the
##' arguments. The rows are appropriately recycled in the \code{nsim}
##' simulations. In order to produce outcomes from the posterior predictive
##' distribution, e.g, pass an \code{object} that contains the data observed so
##' far, \code{truth} contains the \code{prob} function from the model in
##' \code{object}, and \code{args} contains posterior samples from the model.
##' @param firstSeparate enroll the first patient separately from the rest of
##' the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
##' in this patient.
##' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
##' giving the MCMC options for each evaluation in the trial. By default,
##' the standard options are used
##' @param parallel should the simulation runs be parallelized across the
##' clusters of the computer? (not default)
##' @param nCores how many cores should be used for parallel computing?
##' Defaults to the number of cores on the machine, maximum 5.
##' @param \dots not used
##'
##' @example examples/design-method-simulateTDsamplesDesign.R
##'
##' @return an object of class \code{\linkS4class{PseudoSimulations}}
##'
##'  @export
##'  @keywords methods
setMethod("simulate",
  signature =
    signature(
      object = "TDsamplesDesign",
      nsim = "ANY",
      seed = "ANY"
    ),
  def =
    function(object, nsim = 1L, seed = NULL,
             truth, args = NULL, firstSeparate = FALSE,
             mcmcOptions = McmcOptions(),
             parallel = FALSE, nCores =
               min(parallel::detectCores(), 5L),
             ...) {
      ## checks and extracts
      assert_function(truth)
      assert_flag(firstSeparate)
      assert_count(nsim, positive = TRUE)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      args <- as.data.frame(args)
      nArgs <- max(nrow(args), 1L)


      ## seed handling
      RNGstate <- set_seed(seed)

      ## from this,
      ## generate the individual seeds for the simulation runs
      simSeeds <- sample(x = seq_len(1e5), size = as.integer(nsim))

      ## the function to produce the run a single simulation
      ## with index "iterSim"
      runSim <- function(iterSim) {
        ## set the seed for this run
        set.seed(simSeeds[iterSim])

        ## what is now the argument for the truth?
        ## (appropriately recycled)
        thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

        ## so this truth is...
        thisTruth <- function(dose) {
          do.call(
            truth,
            ## First argument: the dose
            c(
              dose,
              ## Following arguments
              thisArgs
            )
          )
        }

        ## start the simulated data with the provided one
        thisData <- object@data

        # In case there are placebo
        if (thisData@placebo) {
          ## what is the probability for tox. at placebo?
          thisProb.PL <- thisTruth(object@data@doseGrid[1])
        }


        ## shall we stop the trial?
        ## First, we want to continue with the starting dose.
        ## This variable is updated after each cohort in the loop.
        stopit <- FALSE

        ## what is the next dose to be used?
        ## initialize with starting dose
        thisDose <- object@startingDose

        ## inside this loop we simulate the whole trial, until stopping
        while (!stopit) {
          ## what is the probability for tox. at this dose?
          thisProb <- thisTruth(thisDose)


          ## what is the cohort size at this dose?
          thisSize <- size(object@cohort_size,
            dose = thisDose,
            data = thisData
          )

          ## In case there are placebo
          if (thisData@placebo) {
            thisSize.PL <- size(object@pl_cohort_size,
              dose = thisDose,
              data = thisData
            )
          }



          ## simulate DLTs: depends on whether we
          ## separate the first patient or not.
          if (firstSeparate && (thisSize > 1L)) {
            ## dose the first patient
            thisDLTs <- rbinom(
              n = 1L,
              size = 1L,
              prob = thisProb
            )

            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisDLTs.PL <- rbinom(
                n = 1L,
                size = 1L,
                prob = thisProb.PL
              )
            }

            ## if there is no DLT:
            if (thisDLTs == 0) {
              ## enroll the remaining patients
              thisDLTs <- c(
                thisDLTs,
                rbinom(
                  n = thisSize - 1L,
                  size = 1L,
                  prob = thisProb
                )
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisDLTs.PL <- c(
                  thisDLTs.PL,
                  rbinom(
                    n = thisSize.PL,
                    size = 1L,
                    prob = thisProb.PL
                  )
                )
              }
            }
          } else {
            ## we can directly dose all patients
            thisDLTs <- rbinom(
              n = thisSize,
              size = 1L,
              prob = thisProb
            )
            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisDLTs.PL <- rbinom(
                n = thisSize.PL,
                size = 1L,
                prob = thisProb.PL
              )
            }
          }

          ## update the data with this placebo (if any) cohort and then with active dose
          if (thisData@placebo && (thisSize.PL > 0L)) {
            thisData <- update(
              object = thisData,
              x = object@data@doseGrid[1],
              y = thisDLTs.PL
            )

            ## update the data with active dose
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs,
              new_cohort = FALSE
            )
          } else {
            ## update the data with this cohort
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs
            )
          }

          ## Update the model with thisData
          thisModel <- update(object@model,
            data = thisData
          )

          ## what is the dose limit?
          doselimit <- maxDose(object@increments,
            data = thisData
          )

          ## generate samples from the model
          thisSamples <- mcmc(
            data = thisData,
            model = thisModel,
            options = mcmcOptions
          )

          ## => what is the next best dose?

          next_bd <- nextBest(object@nextBest,
            doselimit = doselimit,
            samples = thisSamples,
            model = thisModel,
            data = thisData,
            in_sim = TRUE
          )

          thisDose <- next_bd$next_dose_drt
          thisTDtargetDuringTrial <- next_bd$dose_target_drt
          thisTDtargetEndOfTrial <- next_bd$dose_target_eot
          thisTDtargetEndOfTrialatdoseGrid <- next_bd$next_dose_eot
          thisCITDEOT <- list(lower = next_bd$ci_dose_target_eot[1], upper = next_bd$ci_dose_target_eot[2])
          thisratioTDEOT <- next_bd$ci_ratio_dose_target_eot

          ## evaluate stopping rules
          stopit <- stopTrial(object@stopping,
            dose = thisDose,
            samples = thisSamples,
            model = thisModel,
            data = thisData
          )
          stopit_results <- h_unpack_stopit(stopit)
        }

        ## get the fit
        thisFit <- fit(
          object = thisSamples,
          model = thisModel,
          data = thisData
        )


        ## return the results
        thisResult <-
          list(
            data = thisData,
            dose = thisDose,
            TDtargetDuringTrial = thisTDtargetDuringTrial,
            TDtargetEndOfTrial = thisTDtargetEndOfTrial,
            TDtargetEndOfTrialatdoseGrid = thisTDtargetEndOfTrialatdoseGrid,
            TDtargetDuringTrialatdoseGrid = thisDose,
            CITDEOT = thisCITDEOT,
            ratioTDEOT = thisratioTDEOT,
            fit =
              subset(thisFit,
                select = c(middle, lower, upper)
              ),
            stop =
              attr(
                stopit,
                "message"
              ),
            report_results = stopit_results
          )
        return(thisResult)
      }

      resultList <- get_result_list(
        fun = runSim,
        nsim = nsim,
        vars =
          c(
            "simSeeds",
            "args",
            "nArgs",
            "firstSeparate",
            "truth",
            "object",
            "mcmcOptions"
          ),
        parallel = parallel,
        n_cores = nCores
      )

      ## put everything in the Simulations format:

      ## setup the list for the simulated data objects
      dataList <- lapply(resultList, "[[", "data")

      ## set up list for the final TD during Trial Estimate
      TDtargetDuringTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrial"))

      ## set up list for the final TD End of Trial Estimate
      TDtargetEndOfTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrial"))

      ## set up list for the final TD during Trial estimate at dose Grid
      TDtargetDuringTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrialatdoseGrid"))

      ## set up list for the final TD End Of Trial estimate at dose Grid
      TDtargetEndOfTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrialatdoseGrid"))


      ## the vector of the final dose recommendations
      recommendedDoses <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrialatdoseGrid"))

      ## Set up the list for the final 95% CI obtained
      CIList <- lapply(resultList, "[[", "CITDEOT")

      ## Set up the list for the final ratios obtained
      ratioList <- as.numeric(sapply(resultList, "[[", "ratioTDEOT"))

      ## Set up the list for the final TDEOT 95% CI obtained
      CITDEOTList <- lapply(resultList, "[[", "CITDEOT")

      ## Set up the list for the final TDEOT ratios obtained
      ratioTDEOTList <- as.numeric(sapply(resultList, "[[", "ratioTDEOT"))

      ## setup the list for the final fits
      fitList <- lapply(resultList, "[[", "fit")

      ## the reasons for stopping
      stopReasons <- lapply(resultList, "[[", "stop")

      # individual stopping rule results as matrix, labels as column names
      stop_results <- lapply(resultList, "[[", "report_results")
      stop_report <- as.matrix(do.call(rbind, stop_results))


      ## return the results in the Simulations class object
      ret <- PseudoSimulations(
        data = dataList,
        doses = recommendedDoses,
        fit = fitList,
        final_td_target_during_trial_estimates = TDtargetDuringTrialList,
        final_td_target_end_of_trial_estimates = TDtargetEndOfTrialList,
        final_td_target_during_trial_at_dose_grid = TDtargetDuringTrialDoseGridList,
        final_td_target_end_of_trial_at_dose_grid = TDtargetEndOfTrialDoseGridList,
        final_cis = CIList,
        final_ratios = ratioList,
        final_tdeot_cis = CITDEOTList,
        final_tdeot_ratios = ratioTDEOTList,
        stop_reasons = stopReasons,
        stop_report = stop_report,
        seed = RNGstate
      )

      return(ret)
    }
)
## -------------------------------------------------------------------------------------
## Simulate design using DLE responses only without samples (pseudo DLE model)
## --------------------------------------------------------------------------------
###
##' This is a methods to simulate dose escalation procedure only using the DLE responses.
##' This is a method based on the \code{\linkS4class{TDDesign}} where model used are of
##' \code{\linkS4class{ModelTox}} class object and no samples are involved.
##'
##' @param object the \code{\linkS4class{TDDesign}} object we want to simulate the data from
##' @param nsim the number of simulations (default :1)
##' @param seed see \code{\link{set_seed}}
##' @param truth a function which takes as input a dose (vector) and returns the true probability
##' (vector) of the occurrence of a DLE. Additional arguments can be supplied in \code{args}.
##' @param args data frame with arguments for the \code{truth} function. The
##' column names correspond to the argument names, the rows to the values of the
##' arguments. The rows are appropriately recycled in the \code{nsim}
##' simulations. In order to produce outcomes from the posterior predictive
##' distribution, e.g, pass an \code{object} that contains the data observed so
##' far, \code{truth} contains the \code{prob} function from the model in
##' \code{object}, and \code{args} contains posterior samples from the model.
##' @param firstSeparate enroll the first patient separately from the rest of
##' the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
##' in this patient.
##' @param parallel should the simulation runs be parallelized across the
##' clusters of the computer? (not default)
##' @param nCores how many cores should be used for parallel computing?
##' Defaults to the number of cores on the machine, maximum 5.
##' @param \dots not used
##'
##' @example examples/design-method-simulateTDDesign.R
##'
##' @return an object of class \code{\linkS4class{PseudoSimulations}}
##'
##'  @export
##'  @keywords methods
setMethod("simulate",
  signature =
    signature(
      object = "TDDesign",
      nsim = "ANY",
      seed = "ANY"
    ),
  def =
    function(object, nsim = 1L, seed = NULL,
             truth, args = NULL, firstSeparate = FALSE,
             parallel = FALSE, nCores =
               min(parallel::detectCores(), 5L),
             ...) {
      ## checks and extracts
      assert_function(truth)
      assert_flag(firstSeparate)
      assert_count(nsim, positive = TRUE)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      args <- as.data.frame(args)
      nArgs <- max(nrow(args), 1L)

      ## seed handling
      RNGstate <- set_seed(seed)

      ## from this,
      ## generate the individual seeds for the simulation runs
      simSeeds <- sample(x = seq_len(1e5), size = as.integer(nsim))

      ## the function to produce the run a single simulation
      ## with index "iterSim"
      runSim <- function(iterSim) {
        ## set the seed for this run
        set.seed(simSeeds[iterSim])

        ## what is now the argument for the truth?
        ## (appropriately recycled)
        thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

        ## so this truth is...
        thisTruth <- function(dose) {
          do.call(
            truth,
            ## First argument: the dose
            c(
              dose,
              ## Following arguments
              thisArgs
            )
          )
        }

        ## start the simulated data with the provided one
        thisData <- object@data

        # In case there are placebo
        if (thisData@placebo) {
          ## what is the probability for tox. at placebo?
          thisProb.PL <- thisTruth(object@data@doseGrid[1])
        }

        ## shall we stop the trial?
        ## First, we want to continue with the starting dose.
        ## This variable is updated after each cohort in the loop.
        stopit <- FALSE

        ## what is the next dose to be used?
        ## initialize with starting dose
        thisDose <- object@startingDose

        ## inside this loop we simulate the whole trial, until stopping
        while (!stopit) {
          ## what is the probability for tox. at this dose?
          thisProb <- thisTruth(thisDose)

          ## what is the cohort size at this dose?
          thisSize <- size(object@cohort_size,
            dose = thisDose,
            data = thisData
          )

          ## In case there are placebo
          if (thisData@placebo) {
            thisSize.PL <- size(object@pl_cohort_size,
              dose = thisDose,
              data = thisData
            )
          }

          ## simulate DLTs: depends on whether we
          ## separate the first patient or not.
          if (firstSeparate && (thisSize > 1L)) {
            ## dose the first patient
            thisDLTs <- rbinom(
              n = 1L,
              size = 1L,
              prob = thisProb
            )

            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisDLTs.PL <- rbinom(
                n = 1L,
                size = 1L,
                prob = thisProb.PL
              )
            }

            ## if there is no DLT:
            if (thisDLTs == 0) {
              ## enroll the remaining patients
              thisDLTs <- c(
                thisDLTs,
                rbinom(
                  n = thisSize - 1L,
                  size = 1L,
                  prob = thisProb
                )
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisDLTs.PL <- c(
                  thisDLTs.PL,
                  rbinom(
                    n = thisSize.PL,
                    size = 1L,
                    prob = thisProb.PL
                  )
                )
              }
            }
          } else {
            ## we can directly dose all patients
            thisDLTs <- rbinom(
              n = thisSize,
              size = 1L,
              prob = thisProb
            )

            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisDLTs.PL <- rbinom(
                n = thisSize.PL,
                size = 1L,
                prob = thisProb.PL
              )
            }
          }

          ## update the data with this placebo (if any) cohort and then with active dose
          if (thisData@placebo && (thisSize.PL > 0L)) {
            thisData <- update(
              object = thisData,
              x = object@data@doseGrid[1],
              y = thisDLTs.PL
            )

            ## update the data with active dose
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs,
              new_cohort = FALSE
            )
          } else {
            ## update the data with this cohort
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs
            )
          }

          ## Update model estimates with thisData
          thisModel <- update(object@model,
            data = thisData
          )

          ## what is the dose limit?
          doselimit <- maxDose(object@increments, data = thisData)


          ## => what is the next best dose?
          next_bd <- nextBest(object@nextBest,
            doselimit = doselimit,
            model = thisModel,
            data = thisData,
            in_sim = TRUE
          )

          thisDose <- next_bd$next_dose_drt
          thisTDtargetDuringTrial <- next_bd$dose_target_drt
          thisTDtargetEndOfTrial <- next_bd$dose_target_eot
          thisTDtargetEndOfTrialatdoseGrid <- next_bd$next_dose_eot
          thisCITDEOT <- list(lower = next_bd$ci_dose_target_eot[1], upper = next_bd$ci_dose_target_eot[2])
          thisratioTDEOT <- next_bd$ci_ratio_dose_target_eot

          ## evaluate stopping rules
          stopit <- stopTrial(object@stopping,
            dose = thisDose,
            model = thisModel,
            data = thisData
          )
          stopit_results <- h_unpack_stopit(stopit)
        }
        ## get the fit
        prob_fun <- probFunction(thisModel, phi1 = thisModel@phi1, phi2 = thisModel@phi2)
        thisFit <- list(
          phi1 = thisModel@phi1,
          phi2 = thisModel@phi2,
          probDLE = prob_fun(object@data@doseGrid)
        )



        ## return the results
        thisResult <-
          list(
            data = thisData,
            dose = thisDose,
            TDtargetDuringTrial = thisTDtargetDuringTrial,
            TDtargetEndOfTrial = thisTDtargetEndOfTrial,
            TDtargetEndOfTrialatdoseGrid = thisTDtargetEndOfTrialatdoseGrid,
            TDtargetDuringTrialatdoseGrid = thisDose,
            CITDEOT = thisCITDEOT,
            ratioTDEOT = thisratioTDEOT,
            fit = thisFit,
            stop =
              attr(
                stopit,
                "message"
              ),
            report_results = stopit_results
          )
        return(thisResult)
      }


      resultList <- get_result_list(
        fun = runSim,
        nsim = nsim,
        vars =
          c(
            "simSeeds",
            "args",
            "nArgs",
            "firstSeparate",
            "truth",
            "object"
          ),
        parallel = parallel,
        n_cores = nCores
      )

      ## put everything in the Simulations format:

      ## setup the list for the simulated data objects
      dataList <- lapply(resultList, "[[", "data")


      ## set up list for the final TD during Trial Estimate
      TDtargetDuringTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrial"))

      ## set up list for the final TD End of Trial Estimate
      TDtargetEndOfTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrial"))

      ## set up list for the final TD during Trial estimate at dose Grid
      TDtargetDuringTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrialatdoseGrid"))

      ## set up list for the final TD End Of Trial estimate at dose Grid
      TDtargetEndOfTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrialatdoseGrid"))


      ## the vector of the final dose recommendations
      recommendedDoses <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrialatdoseGrid"))

      ## Set up the list for the final 95% CI obtained
      CIList <- lapply(resultList, "[[", "CITDEOT")

      ## Set up the list for the final ratios obtained
      ratioList <- as.numeric(sapply(resultList, "[[", "ratioTDEOT"))

      ## Set up the list for the final TDEOT 95% CI obtained
      CITDEOTList <- lapply(resultList, "[[", "CITDEOT")

      ## Set up the list for the final TDEOT ratios obtained
      ratioTDEOTList <- as.numeric(sapply(resultList, "[[", "ratioTDEOT"))
      ## set up the list for the final fits

      fitList <- lapply(resultList, "[[", "fit")

      ## the reasons for stopping
      stopReasons <- lapply(resultList, "[[", "stop")

      # individual stopping rule results as matrix, labels as column names
      stop_results <- lapply(resultList, "[[", "report_results")
      stop_report <- as.matrix(do.call(rbind, stop_results))


      ## return the results in the Simulations class object
      ret <- PseudoSimulations(
        data = dataList,
        doses = recommendedDoses,
        fit = fitList,
        final_td_target_during_trial_estimates = TDtargetDuringTrialList,
        final_td_target_end_of_trial_estimates = TDtargetEndOfTrialList,
        final_td_target_during_trial_at_dose_grid = TDtargetDuringTrialDoseGridList,
        final_td_target_end_of_trial_at_dose_grid = TDtargetEndOfTrialDoseGridList,
        final_cis = CIList,
        final_ratios = ratioList,
        final_tdeot_cis = CITDEOTList,
        final_tdeot_ratios = ratioTDEOTList,
        stop_reasons = stopReasons,
        stop_report = stop_report,
        seed = RNGstate
      )

      return(ret)
    }
)
## -----------------------------------------------------------------------------------------------
## Simulate design using DLE and efficacy responses without DLE and efficacy samples (pseudo models)
## --------------------------------------------------------------------------------------------
###
##' This is a methods to simulate dose escalation procedure using both DLE and efficacy responses.
##' This is a method based on the \code{\linkS4class{DualResponsesDesign}} where DLEmodel used are of
##' \code{\linkS4class{ModelTox}} class object and efficacy model used are of \code{\linkS4class{ModelEff}}
##' class object. In addition, no DLE and efficacy samples are involved or generated in the simulation
##' process
##'
##' @param object the \code{\linkS4class{DualResponsesDesign}} object we want to simulate the data from
##' @param nsim the number of simulations (default :1)
##' @param seed see \code{\link{set_seed}}
##' @param trueDLE a function which takes as input a dose (vector) and returns the true probability
##' (vector) of the occurrence of a DLE. Additional arguments can be supplied in \code{args}.
##' @param trueEff a function which takes as input a dose (vector) and returns the expected efficacy
##' responses (vector). Additional arguments can be supplied in \code{args}.
##' @param trueNu the precision, the inverse of the variance of the efficacy responses
##' @param args data frame with arguments for the \code{trueDLE} and
##' \code{trueEff} function. The column names correspond to the argument
##' names, the rows to the values of the arguments. The rows are appropriately
##' recycled in the \code{nsim} simulations.
##' @param firstSeparate enroll the first patient separately from the rest of
##' the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
##' in this patient.
##' @param parallel should the simulation runs be parallelized across the
##' clusters of the computer? (not default)
##' @param nCores how many cores should be used for parallel computing?
##' Defaults to the number of cores on the machine, maximum 5.
##' @param \dots not used
##'
##' @example examples/design-method-simulateDualResponsesDesign.R
##'
##' @return an object of class \code{\linkS4class{PseudoDualSimulations}}
##'
##' @export
##' @keywords methods

setMethod("simulate",
  signature =
    signature(
      object = "DualResponsesDesign",
      nsim = "ANY",
      seed = "ANY"
    ),
  def =
    function(object, nsim = 1L, seed = NULL,
             trueDLE, trueEff, trueNu,
             args = NULL, firstSeparate = FALSE,
             parallel = FALSE, nCores =
               min(parallel::detectCores(), 5L),
             ...) {
      ## checks and extracts
      assert_function(trueDLE)
      assert_function(trueEff)
      assert_true(trueNu > 0)
      assert_flag(firstSeparate)
      assert_count(nsim, positive = TRUE)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      args <- as.data.frame(args)
      nArgs <- max(nrow(args), 1L)

      ## get names of arguments (excluding the first one which is the dose)
      trueDLEArgnames <- names(formals(trueDLE))[-1]
      trueEffArgnames <- names(formals(trueEff))[-1]



      ## seed handling
      RNGstate <- set_seed(seed)

      ## from this,
      ## generate the individual seeds for the simulation runs
      simSeeds <- sample(x = seq_len(1e5), size = as.integer(nsim))

      ## the function to produce the run a single simulation
      ## with index "iterSim"
      runSim <- function(iterSim) {
        ## set the seed for this run
        set.seed(simSeeds[iterSim])

        ## what is now the argument for the truth?
        ## (appropriately recycled)
        thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

        ## so this truth DLE function is...
        thisTruthDLE <- function(dose) {
          do.call(
            trueDLE,
            ## First argument: the dose
            c(
              dose,
              ## Following arguments: take only those that
              ## are required by the DLE function
              as.list(thisArgs)[trueDLEArgnames]
            )
          )
        }

        ## and the truth Eff function is:
        thisTruthEff <- function(dose) {
          do.call(
            trueEff,
            ## First argument: the dose
            c(
              dose,
              ## Following arguments: take only those that
              ## are required by the Eff function
              as.list(thisArgs)[trueEffArgnames]
            )
          )
        }

        ## start the simulated data with the provided one
        thisData <- object@data

        ## find true sigma2 to generate responses

        trueSigma2 <- 1 / trueNu

        ## start the simulated data with the provided one
        thisData <- object@data

        if (thisData@placebo) {
          ## what is the probability for tox. at placebo?
          thisProb.PL <- thisTruthDLE(object@data@doseGrid[1])

          ## what is the mean efficacy at placebo?
          thisMeanEff.PL <- thisTruthEff(object@data@doseGrid[1])
        }

        ## shall we stop the trial?
        ## First, we want to continue with the starting dose.
        ## This variable is updated after each cohort in the loop.
        stopit <- FALSE

        ## what is the next dose to be used?
        ## initialize with starting dose
        thisDose <- object@startingDose

        ## inside this loop we simulate the whole trial, until stopping
        while (!stopit) {
          ## what is the probability for tox. at this dose?
          thisDLEProb <- thisTruthDLE(thisDose)
          thisMeanEff <- thisTruthEff(thisDose)

          ## what is the cohort size at this dose?
          thisSize <- size(object@cohort_size,
            dose = thisDose,
            data = thisData
          )


          ## In case there are placebo
          ## what is the cohort size at this dose for Placebo?
          if (thisData@placebo) {
            thisSize.PL <- size(object@pl_cohort_size,
              dose = thisDose,
              data = thisData
            )
          }


          ## simulate DLTs: depends on whether we
          ## separate the first patient or not.
          if (firstSeparate && (thisSize > 1L)) {
            ## dose the first patient
            thisDLTs <- rbinom(
              n = 1L,
              size = 1L,
              prob = thisDLEProb
            )


            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisDLTs.PL <- rbinom(
                n = 1L,
                size = 1L,
                prob = thisProb.PL
              )
            }

            thisEff <- rnorm(
              n = 1L,
              mean = thisMeanEff,
              sd = sqrt(trueSigma2)
            )

            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisEff.PL <- rnorm(
                n = 1L,
                mean = thisMeanEff.PL,
                sd = sqrt(trueSigma2)
              )
            }

            ## if there is no DLT:
            if (thisDLTs == 0) {
              ## enroll the remaining patients
              thisDLTs <- c(
                thisDLTs,
                rbinom(
                  n = thisSize - 1L,
                  size = 1L,
                  prob = thisDLEProb
                )
              )

              thisEff <- c(
                thisEff,
                rnorm(
                  n = thisSize - 1L,
                  mean = thisMeanEff,
                  sd = sqrt(trueSigma2)
                )
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisDLTs.PL <- c(
                  thisDLTs.PL,
                  rbinom(
                    n = thisSize.PL,
                    size = 1L,
                    prob = thisProb.PL
                  )
                )
                thisEff.PL <- c(
                  thisMeanEff.PL,
                  rnorm(
                    n = thisSize.PL,
                    mean = thisMeanEff.PL,
                    sd = sqrt(trueSigma2)
                  )
                )
              }
            }
          } else {
            ## we can directly dose all patients
            thisDLTs <- rbinom(
              n = thisSize,
              size = 1L,
              prob = thisDLEProb
            )
            thisEff <- rnorm(
              n = thisSize,
              mean = thisMeanEff,
              sd = sqrt(trueSigma2)
            )

            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisDLTs.PL <- rbinom(
                n = thisSize.PL,
                size = 1L,
                prob = thisProb.PL
              )
              thisEff.PL <- rnorm(
                n = thisSize.PL,
                mean = thisMeanEff.PL,
                sd = sqrt(trueSigma2)
              )
            }
          }

          ## update the data with this placebo (if any) cohort and then with active dose
          if (thisData@placebo && (thisSize.PL > 0L)) {
            thisData <- update(
              object = thisData,
              x = object@data@doseGrid[1],
              y = thisDLTs.PL,
              w = thisEff.PL,
              check = FALSE
            )

            ## update the data with active dose
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs,
              w = thisEff,
              new_cohort = FALSE
            )
          } else {
            ## update the data with this cohort
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs,
              w = thisEff
            )
          }
          ## Update model estimate in DLE and Eff models
          thisDLEModel <- update(
            object = object@model,
            data = thisData
          )

          thisEffModel <- update(
            object = object@eff_model,
            data = thisData
          )

          thisNu <- thisEffModel@nu


          thisSigma2 <- if (thisEffModel@use_fixed) {
            1 / thisNu
          } else {
            1 / (as.numeric(thisNu["a"] / thisNu["b"]))
          }


          ## what is the dose limit?
          doselimit <- maxDose(object@increments, data = thisData)



          ## => what is the next best dose?
          next_bd <- nextBest(object@nextBest,
            doselimit = doselimit,
            model = thisDLEModel,
            data = thisData,
            model_eff = thisEffModel,
            in_sim = TRUE
          )

          thisDose <- next_bd$next_dose
          thisTDtargetDuringTrial <- next_bd$dose_target_drt
          thisTDtargetDuringTrialAtDoseGrid <- next_bd$next_dose_drt
          thisTDtargetEndOfTrial <- next_bd$dose_target_eot
          thisTDtargetEndOfTrialAtDoseGrid <- next_bd$next_dose_eot
          thisGstar <- next_bd$dose_max_gain
          thisGstarAtDoseGrid <- next_bd$next_dose_max_gain

          Recommend <- min(thisTDtargetEndOfTrialAtDoseGrid, thisGstarAtDoseGrid)

          ## Find the 95 % CI and its ratio (upper to the lower of this 95% CI of each of the estimates)
          thisCITDEOT <- list(lower = next_bd$ci_dose_target_eot[1], upper = next_bd$ci_dose_target_eot[2])
          thisratioTDEOT <- next_bd$ci_ratio_dose_target_eot

          thisCIGstar <- list(lower = next_bd$ci_dose_max_gain[1], upper = next_bd$ci_dose_max_gain[2])
          thisratioGstar <- next_bd$ci_ratio_dose_max_gain

          ## Find the optimal dose
          OptimalDose <- min(thisGstar, thisTDtargetEndOfTrial)

          if (OptimalDose == thisGstar) {
            thisratio <- thisratioGstar
            thisCI <- thisCIGstar
          } else {
            thisratio <- thisratioTDEOT
            thisCI <- thisCITDEOT
          }

          ## evaluate stopping rules
          stopit <- stopTrial(object@stopping,
            dose = thisDose,
            model = thisDLEModel,
            data = thisData,
            Effmodel = thisEffModel
          )
          stopit_results <- h_unpack_stopit(stopit)
        }

        ## get the fits
        prob_fun <- probFunction(thisDLEModel, phi1 = thisDLEModel@phi1, phi2 = thisDLEModel@phi2)
        thisDLEFit <- list(
          phi1 = thisDLEModel@phi1,
          phi2 = thisDLEModel@phi2,
          probDLE = prob_fun(object@data@doseGrid)
        )

        eff_fun <- efficacyFunction(thisEffModel, theta1 = thisEffModel@theta1, theta2 = thisEffModel@theta2)
        thisEffFit <- list(
          theta1 = thisEffModel@theta1,
          theta2 = thisEffModel@theta2,
          ExpEff = eff_fun(object@data@doseGrid)
        )


        ## return the results
        thisResult <- list(
          data = thisData,
          dose = thisDose,
          TDtargetDuringTrial = thisTDtargetDuringTrial,
          TDtargetDuringTrialAtDoseGrid = thisTDtargetDuringTrialAtDoseGrid,
          TDtargetEndOfTrial = thisTDtargetEndOfTrial,
          TDtargetEndOfTrialAtDoseGrid = thisTDtargetEndOfTrialAtDoseGrid,
          Gstar = thisGstar,
          GstarAtDoseGrid = thisGstarAtDoseGrid,
          Recommend = Recommend,
          OptimalDose = OptimalDose,
          OptimalDoseAtDoseGrid = Recommend,
          ratio = thisratio,
          CI = thisCI,
          ratioGstar = thisratioGstar,
          CIGstar = thisCIGstar,
          ratioTDEOT = thisratioTDEOT,
          CITDEOT = thisCITDEOT,
          fitDLE = thisDLEFit,
          fitEff = thisEffFit,
          sigma2est = thisSigma2,
          stop = attr(
            stopit,
            "message"
          ),
          report_results = stopit_results
        )

        return(thisResult)
      }


      resultList <- get_result_list(
        fun = runSim,
        nsim = nsim,
        vars =
          c(
            "simSeeds",
            "args",
            "nArgs",
            "firstSeparate",
            "trueDLE",
            "trueEff",
            "trueNu",
            "object"
          ),
        parallel = parallel,
        n_cores = nCores
      )


      ## put everything in the Simulations format:

      ## setup the list for the simulated data objects
      dataList <- lapply(resultList, "[[", "data")

      ## the vector of the final dose recommendations
      recommendedDoses <- as.numeric(sapply(resultList, "[[", "Recommend"))


      ## set up list for the final TD during Trial Estimate
      TDtargetDuringTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrial"))

      ## set up list for the final TD End of Trial Estimate
      TDtargetEndOfTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrial"))

      ## set up list for the final TD during Trial estimate at dose Grid
      TDtargetDuringTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrialAtDoseGrid"))

      ## set up list for the final TD End Of Trial estimate at dose Grid
      TDtargetEndOfTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrialAtDoseGrid"))

      ## set up list for the final Gstar estimates
      GstarList <- as.numeric(sapply(resultList, "[[", "Gstar"))

      ## set up list for the final Gstar estimates at dose grid
      GstarAtDoseGridList <- as.numeric(sapply(resultList, "[[", "GstarAtDoseGrid"))

      ## set up list for final optimal dose estimates
      OptimalDoseList <- as.numeric(sapply(resultList, "[[", "OptimalDose"))

      ## set up list for final optimal dose estimates at dose Grid
      OptimalDoseAtDoseGridList <- as.numeric(sapply(resultList, "[[", "Recommend"))

      ## Set up the list for the final 95% CI obtained
      CIList <- lapply(resultList, "[[", "CI")

      ## Set up the list for the final ratios obtained
      ratioList <- as.numeric(sapply(resultList, "[[", "ratio"))

      ## Set up the list for the final 95% CI of the TDtarget End Of Trial obtained
      CITDEOTList <- lapply(resultList, "[[", "CITDEOT")

      ## Set up the list for the final ratios of the TDtarget End Of Trial obtained
      ratioTDEOTList <- as.numeric(sapply(resultList, "[[", "ratioTDEOT"))

      ## Set up the list for the final 95% CI of the Gstar obtained
      CIGstarList <- lapply(resultList, "[[", "CIGstar")

      ## Set up the list for the final ratios of the Gstar obtained
      ratioGstarList <- as.numeric(sapply(resultList, "[[", "ratioGstar"))



      ## set up the list for the final fits
      fitDLEList <- lapply(resultList, "[[", "fitDLE")
      fitEffList <- lapply(resultList, "[[", "fitEff")


      ## the vector of the sigma2
      sigma2Estimates <- as.numeric(sapply(resultList, "[[", "sigma2est"))

      ## the reasons for stopping
      stopReasons <- lapply(resultList, "[[", "stop")

      # individual stopping rule results as matrix, labels as column names
      stop_results <- lapply(resultList, "[[", "report_results")
      stop_report <- as.matrix(do.call(rbind, stop_results))


      ## return the results in the Simulations class object
      ret <- PseudoDualSimulations(
        data = dataList,
        doses = recommendedDoses,
        final_td_target_during_trial_estimates = TDtargetDuringTrialList,
        final_td_target_end_of_trial_estimates = TDtargetEndOfTrialList,
        final_td_target_during_trial_at_dose_grid = TDtargetDuringTrialDoseGridList,
        final_td_target_end_of_trial_at_dose_grid = TDtargetEndOfTrialDoseGridList,
        final_cis = CIList,
        final_ratios = ratioList,
        final_gstar_estimates = GstarList,
        final_gstar_at_dose_grid = GstarAtDoseGridList,
        final_gstar_cis = CIGstarList,
        final_gstar_ratios = ratioGstarList,
        final_tdeot_cis = CITDEOTList,
        final_tdeot_ratios = ratioTDEOTList,
        final_optimal_dose = OptimalDoseList,
        final_optimal_dose_at_dose_grid = OptimalDoseAtDoseGridList,
        fit = fitDLEList,
        fit_eff = fitEffList,
        sigma2_est = sigma2Estimates,
        stop_reasons = stopReasons,
        stop_report = stop_report,
        seed = RNGstate
      )
      return(ret)
    }
)

## =========================================================================
## -----------------------------------------------------------------------------------------------
## Simulate design using DLE and efficacy responses with DLE and efficacy samples (pseudo models)
## --------------------------------------------------------------------------------------------
###
##' This is a methods to simulate dose escalation procedure using both DLE and efficacy responses.
##' This is a method based on the \code{\linkS4class{DualResponsesSamplesDesign}} where DLEmodel
##' used are of
##' \code{\linkS4class{ModelTox}} class object and efficacy model used are of
##' \code{\linkS4class{ModelEff}}
##' class object (special case is \code{\linkS4class{EffFlexi}} class model object).
##' In addition, DLE and efficacy samples are involved or generated in the simulation
##' process
##'
##' @param object the \code{\linkS4class{DualResponsesSamplesDesign}} object we want to
##' simulate the data from
##' @param nsim the number of simulations (default :1)
##' @param seed see \code{\link{set_seed}}
##' @param trueDLE a function which takes as input a dose (vector) and returns the true probability
##' (vector) of the occurrence of a DLE. Additional arguments can be supplied in \code{args}.
##' @param trueEff a function which takes as input a dose (vector) and returns the expected
##' efficacy responses (vector). Additional arguments can be supplied in \code{args}.
##' @param trueNu (not with code{\linkS4class{EffFlexi}}) the precision, the inverse of the
##' variance of the efficacy responses
##' @param trueSigma2 (only with code{\linkS4class{EffFlexi}}) the true variance of the efficacy
##' responses which must be a single positive scalar.
##' @param trueSigma2betaW (only with code{\linkS4class{EffFlexi}}) the true variance for the
##' random walk model used for smoothing. This must be a single positive scalar.
##' @param args data frame with arguments for the \code{trueDLE} and
##' \code{trueEff} function. The column names correspond to the argument
##' names, the rows to the values of the arguments. The rows are appropriately
##' recycled in the \code{nsim} simulations.
##' @param firstSeparate enroll the first patient separately from the rest of
##' the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
##' in this patient.
##' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
##' giving the MCMC options for each evaluation in the trial. By default,
##' the standard options are used
##' @param parallel should the simulation runs be parallelized across the
##' clusters of the computer? (not default)
##' @param nCores how many cores should be used for parallel computing?
##' Defaults to the number of cores on the machine, maximum 5.
##'
##' @param ... not used.
##'
##' @example examples/design-method-simulateDualResponsesSamplesDesign.R
##'
##' @return an object of class \code{\linkS4class{PseudoDualSimulations}} or
##' \code{\linkS4class{PseudoDualFlexiSimulations}}
##'
##' @export
##' @keywords methods
setMethod("simulate",
  signature =
    signature(
      object = "DualResponsesSamplesDesign",
      nsim = "ANY",
      seed = "ANY"
    ),
  def =
    function(object, nsim = 1L, seed = NULL,
             trueDLE, trueEff, trueNu = NULL,
             trueSigma2 = NULL, trueSigma2betaW = NULL,
             args = NULL, firstSeparate = FALSE,
             mcmcOptions = McmcOptions(),
             parallel = FALSE, nCores =
               min(parallel::detectCores(), 5L),
             ...) {
      ## common checks and extracts
      assert_function(trueDLE)
      assert_flag(firstSeparate)
      assert_count(nsim, positive = TRUE)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      ## check if special case applies
      isFlexi <- is(object@eff_model, "EffFlexi")

      ## conditional code from here on:
      if (isFlexi) {
        ## special checks and extracts
        stopifnot(
          trueSigma2 > 0,
          trueSigma2betaW > 0,
          is.numeric(trueEff),
          length(trueEff) == length(object@data@doseGrid)
        )

        args <- as.data.frame(args)
        nArgs <- max(nrow(args), 1L)

        ## get names of arguments (excluding the first one which is the dose)
        trueDLEArgnames <- names(formals(trueDLE))[-1]

        ## seed handling
        RNGstate <- set_seed(seed)

        ## from this,
        ## generate the individual seeds for the simulation runs
        simSeeds <- sample(x = seq_len(1e5), size = as.integer(nsim))

        ## the function to produce the run a single simulation
        ## with index "iterSim"
        runSim <- function(iterSim) {
          ## set the seed for this run
          set.seed(simSeeds[iterSim])

          ## what is now the argument for the truth?
          ## (appropriately recycled)
          thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

          ## so this truth is...
          thisTruthDLE <- function(dose) {
            do.call(
              trueDLE,
              ## First argument: the dose
              c(
                dose,
                ## Following arguments
                thisArgs
              )
            )
          }

          ## get the true Eff
          thisTruthEff <- trueEff

          ## start the simulated data with the provided one
          thisData <- object@data

          ## shall we stop the trial?
          ## First, we want to continue with the starting dose.
          ## This variable is updated after each cohort in the loop.
          stopit <- FALSE

          ## what is the next dose to be used?
          ## initialize with starting dose
          thisDose <- object@startingDose

          ## Start with specified sigma2 and sigma2betaW
          thisSigma2 <- trueSigma2
          thisSigma2betaW <- trueSigma2betaW


          ## inside this loop we simulate the whole trial, until stopping
          while (!stopit) {
            ## what is the probability for tox. at this dose?
            thisDLEProb <- thisTruthDLE(thisDose)
            thisDoseIndex <- which(thisDose == thisData@doseGrid)
            thisMeanEff <- thisTruthEff[thisDoseIndex]



            ## what is the cohort size at this dose?
            thisSize <- size(object@cohort_size,
              dose = thisDose,
              data = thisData
            )

            if (thisData@placebo) {
              thisSize.PL <- size(object@pl_cohort_size,
                dose = thisDose,
                data = thisData
              )
            }

            ## simulate DLTs: depends on whether we
            ## separate the first patient or not.
            if (firstSeparate && (thisSize > 1L)) {
              ## dose the first patient
              thisDLTs <- rbinom(
                n = 1L,
                size = 1L,
                prob = thisDLEProb
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisDLTs.PL <- rbinom(
                  n = 1L,
                  size = 1L,
                  prob = thisProb.PL
                )
              }

              thisEff <- rnorm(
                n = 1L,
                mean = thisMeanEff,
                sd = sqrt(trueSigma2)
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisEff.PL <- rnorm(
                  n = 1L,
                  mean = thisMeanEff.PL,
                  sd = sqrt(trueSigma2)
                )
              }

              ## if there is no DLT:
              if (thisDLTs == 0) {
                ## enroll the remaining patients
                thisDLTs <- c(
                  thisDLTs,
                  rbinom(
                    n = thisSize - 1L,
                    size = 1L,
                    prob = thisDLEProb
                  )
                )
                thisEff <- c(
                  thisEff,
                  rnorm(
                    n = thisSize - 1L,
                    mean = thisMeanEff,
                    sd = sqrt(trueSigma2)
                  )
                )
                if (thisData@placebo && (thisSize.PL > 0L)) {
                  thisDLTs.PL <- c(
                    thisDLTs.PL,
                    rbinom(
                      n = thisSize.PL,
                      size = 1L,
                      prob = thisProb.PL
                    )
                  )
                  thisEff.PL <- c(
                    thisMeanEff.PL,
                    rnorm(
                      n = thisSize.PL,
                      mean = thisMeanEff.PL,
                      sd = sqrt(trueSigma2)
                    )
                  )
                }
              }
            } else {
              ## we can directly dose all patients
              thisDLTs <- rbinom(
                n = thisSize,
                size = 1L,
                prob = thisDLEProb
              )

              thisEff <- rnorm(
                n = thisSize,
                mean = thisMeanEff,
                sd = sqrt(trueSigma2)
              )
              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisDLTs.PL <- rbinom(
                  n = thisSize.PL,
                  size = 1L,
                  prob = thisProb.PL
                )
                thisEff.PL <- rnorm(
                  n = thisSize.PL,
                  mean = thisMeanEff.PL,
                  sd = sqrt(trueSigma2)
                )
              }
            }

            ## update the data with this placebo (if any) cohort and then with active dose
            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisData <- update(
                object = thisData,
                x = object@data@doseGrid[1],
                y = thisDLTs.PL,
                w = thisEff.PL,
                check = FALSE
              )

              ## update the data with active dose
              thisData <- update(
                object = thisData,
                x = thisDose,
                y = thisDLTs,
                w = thisEff,
                new_cohort = FALSE
              )
            } else {
              ## update the data with this cohort
              thisData <- update(
                object = thisData,
                x = thisDose,
                y = thisDLTs,
                w = thisEff
              )
            }

            ## Update model estimate in DLE model
            thisDLEModel <- update(
              object = object@model,
              data = thisData
            )

            thisEffModel <- update(
              object = object@eff_model,
              data = thisData
            )


            ## what is the dose limit?
            doselimit <- maxDose(object@increments,
              data = thisData
            )

            ## generate DLE and Eff samples from the DLE and Eff model
            thisDLEsamples <- mcmc(
              data = thisData,
              model = thisDLEModel,
              options = mcmcOptions
            )

            thisEffsamples <- mcmc(
              data = thisData,
              model = thisEffModel,
              options = mcmcOptions
            )

            thisSigma2 <- mean(thisEffsamples@data$sigma2W)

            thisSigma2betaW <- mean(thisEffsamples@data$sigma2betaW)

            ## => what is the next best dose?

            next_bd <- nextBest(object@nextBest,
              doselimit = doselimit,
              samples = thisDLEsamples,
              model = thisDLEModel,
              model_eff = thisEffModel,
              samples_eff = thisEffsamples,
              data = thisData,
              in_sim = TRUE
            )

            thisDose <- next_bd$next_dose
            thisTDtargetDuringTrial <- next_bd$dose_target_drt
            thisTDtargetDuringTrialAtDoseGrid <- next_bd$next_dose_drt
            thisTDtargetEndOfTrial <- next_bd$dose_target_eot
            thisTDtargetEndOfTrialAtDoseGrid <- next_bd$next_dose_eot
            thisGstar <- next_bd$dose_max_gain
            thisGstarAtDoseGrid <- next_bd$next_dose_max_gain

            Recommend <- min(thisTDtargetEndOfTrialAtDoseGrid, thisGstarAtDoseGrid)

            ## Find the 95 % CI and its ratio (upper to the lower of this 95% CI of each of the estimates)
            thisCITDEOT <- list(lower = next_bd$ci_dose_target_eot[1], upper = next_bd$ci_dose_target_eot[2])
            thisratioTDEOT <- next_bd$ci_ratio_dose_target_eot

            thisCIGstar <- list(lower = next_bd$ci_dose_max_gain[1], upper = next_bd$ci_dose_max_gain[2])
            thisratioGstar <- next_bd$ci_ratio_dose_max_gain

            ## Find the optimal dose
            OptimalDose <- min(thisGstar, thisTDtargetEndOfTrial)

            if (OptimalDose == thisGstar) {
              thisratio <- thisratioGstar
              thisCI <- thisCIGstar
            } else {
              thisratio <- thisratioTDEOT
              thisCI <- thisCITDEOT
            }

            ## evaluate stopping rules
            stopit <- stopTrial(object@stopping,
              dose = thisDose,
              samples = thisDLEsamples,
              model = thisDLEModel,
              data = thisData,
              TDderive = object@nextBest@derive,
              Effmodel = thisEffModel,
              Effsamples = thisEffsamples,
              Gstarderive = object@nextBest@mg_derive
            )
            stopit_results <- h_unpack_stopit(stopit)
          }

          ## get the fits

          thisDLEFit <- fit(
            object = thisDLEsamples,
            model = thisDLEModel,
            data = thisData
          )

          thisEffFit <- fit(
            object = thisEffsamples,
            model = thisEffModel,
            data = thisData
          )


          ## return the results
          thisResult <-
            list(
              data = thisData,
              dose = thisDose,
              TDtargetDuringTrial = thisTDtargetDuringTrial,
              TDtargetDuringTrialAtDoseGrid = thisTDtargetDuringTrialAtDoseGrid,
              TDtargetEndOfTrial = thisTDtargetEndOfTrial,
              TDtargetEndOfTrialAtDoseGrid = thisTDtargetEndOfTrialAtDoseGrid,
              Gstar = thisGstar,
              GstarAtDoseGrid = thisGstarAtDoseGrid,
              Recommend = Recommend,
              OptimalDose = OptimalDose,
              OptimalDoseAtDoseGrid = Recommend,
              ratio = thisratio,
              CI = thisCI,
              ratioGstar = thisratioGstar,
              CIGstar = thisCIGstar,
              ratioTDEOT = thisratioTDEOT,
              CITDEOT = thisCITDEOT,
              fitDLE = subset(thisDLEFit,
                select =
                  c(middle, lower, upper)
              ),
              fitEff = subset(thisEffFit,
                select =
                  c(middle, lower, upper)
              ),
              sigma2est = thisSigma2,
              sigma2betaWest = thisSigma2betaW,
              stop =
                attr(
                  stopit,
                  "message"
                ),
              report_results = stopit_results
            )

          return(thisResult)
        }

        resultList <- get_result_list(
          fun = runSim,
          nsim = nsim,
          vars =
            c(
              "simSeeds",
              "args",
              "nArgs",
              "firstSeparate",
              "trueDLE",
              "trueEff",
              "trueSigma2",
              "trueSigma2betaW",
              "object",
              "mcmcOptions"
            ),
          parallel = parallel,
          n_cores = nCores
        )

        ## put everything in the Simulations format:

        ## setup the list for the simulated data objects
        dataList <- lapply(resultList, "[[", "data")

        ## the vector of the final dose recommendations
        recommendedDoses <- as.numeric(sapply(resultList, "[[", "Recommend"))


        ## set up the list for the final fits for both DLE and efficacy
        fitDLEList <- lapply(resultList, "[[", "fitDLE")
        fitEffList <- lapply(resultList, "[[", "fitEff")

        ## the vector of sigma2 estimates
        sigma2Estimates <- as.numeric(sapply(resultList, "[[", "sigma2est"))

        ## the vector of sigma2betaW estimates
        sigma2betaWEstimates <- as.numeric(sapply(resultList, "[[", "sigma2betaWest"))

        ## set up list for the final TD during Trial Estimate
        TDtargetDuringTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrial"))

        ## set up list for the final TD End of Trial Estimate
        TDtargetEndOfTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrial"))

        ## set up list for the final TD during Trial estimate at dose Grid
        TDtargetDuringTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrialAtDoseGrid"))

        ## set up list for the final TD End Of Trial estimate at dose Grid
        TDtargetEndOfTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrialAtDoseGrid"))

        ## set up list for the final Gstar estimates
        GstarList <- as.numeric(sapply(resultList, "[[", "Gstar"))

        ## set up list for the final Gstar estimates at dose grid
        GstarAtDoseGridList <- as.numeric(sapply(resultList, "[[", "GstarAtDoseGrid"))

        ## set up list for final optimal dose estimates
        OptimalDoseList <- as.numeric(sapply(resultList, "[[", "OptimalDose"))

        ## set up list for final optimal dose estimates at dose Grid
        OptimalDoseAtDoseGridList <- as.numeric(sapply(resultList, "[[", "Recommend"))

        ## Set up the list for the final 95% CI obtained
        CIList <- lapply(resultList, "[[", "CI")

        ## Set up the list for the final ratios obtained
        ratioList <- as.numeric(sapply(resultList, "[[", "ratio"))

        ## Set up the list for the final 95% CI of the TDtarget End Of Trial obtained
        CITDEOTList <- lapply(resultList, "[[", "CITDEOT")

        ## Set up the list for the final ratios of the TDtarget End Of Trial obtained
        ratioTDEOTList <- as.numeric(sapply(resultList, "[[", "ratioTDEOT"))

        ## Set up the list for the final 95% CI of the Gstar obtained
        CIGstarList <- lapply(resultList, "[[", "CIGstar")

        ## Set up the list for the final ratios of the Gstar obtained
        ratioGstarList <- as.numeric(sapply(resultList, "[[", "ratioGstar"))

        ## the reasons for stopping
        stopReasons <- lapply(resultList, "[[", "stop")

        # individual stopping rule results as matrix, labels as column names
        stop_results <- lapply(resultList, "[[", "report_results")
        stop_report <- as.matrix(do.call(rbind, stop_results))


        ## return the results in the Simulations class object
        ret <- PseudoDualFlexiSimulations(
          data = dataList,
          doses = recommendedDoses,
          final_td_target_during_trial_estimates = TDtargetDuringTrialList,
          final_td_target_end_of_trial_estimates = TDtargetEndOfTrialList,
          final_td_target_during_trial_at_dose_grid = TDtargetDuringTrialDoseGridList,
          final_td_target_end_of_trial_at_dose_grid = TDtargetEndOfTrialDoseGridList,
          final_cis = CIList,
          final_ratios = ratioList,
          final_gstar_estimates = GstarList,
          final_gstar_at_dose_grid = GstarAtDoseGridList,
          final_gstar_cis = CIGstarList,
          final_gstar_ratios = ratioGstarList,
          final_tdeot_cis = CITDEOTList,
          final_tdeot_ratios = ratioTDEOTList,
          final_optimal_dose = OptimalDoseList,
          final_optimal_dose_at_dose_grid = OptimalDoseAtDoseGridList,
          fit = fitDLEList,
          fit_eff = fitEffList,
          sigma2_est = sigma2Estimates,
          sigma2_beta_west = sigma2betaWEstimates,
          stop_reasons = stopReasons,
          stop_report = stop_report,
          seed = RNGstate
        )

        return(ret)
      } else {
        stopifnot(
          trueNu > 0,
          is.function(trueEff)
        )


        args <- as.data.frame(args)
        nArgs <- max(nrow(args), 1L)

        ## get names of arguments (excluding the first one which is the dose)
        trueDLEArgnames <- names(formals(trueDLE))[-1]
        trueEffArgnames <- names(formals(trueEff))[-1]



        ## seed handling
        RNGstate <- set_seed(seed)

        ## from this,
        ## generate the individual seeds for the simulation runs
        simSeeds <- sample(x = seq_len(1e5), size = nsim)

        ## the function to produce the run a single simulation
        ## with index "iterSim"
        runSim <- function(iterSim) {
          ## set the seed for this run
          set.seed(simSeeds[iterSim])

          ## what is now the argument for the truth?
          ## (appropriately recycled)
          thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

          ## so this truth DLE function is...
          thisTruthDLE <- function(dose) {
            do.call(
              trueDLE,
              ## First argument: the dose
              c(
                dose,
                ## Following arguments: take only those that
                ## are required by the DLE function
                as.list(thisArgs)[trueDLEArgnames]
              )
            )
          }

          ## and the truth Eff function is:
          thisTruthEff <- function(dose) {
            do.call(
              trueEff,
              ## First argument: the dose
              c(
                dose,
                ## Following arguments: take only those that
                ## are required by the Eff function
                as.list(thisArgs)[trueEffArgnames]
              )
            )
          }

          ## find true sigma2 to generate responses

          trueSigma2 <- 1 / trueNu

          ## start the simulated data with the provided one
          thisData <- object@data


          ## shall we stop the trial?
          ## First, we want to continue with the starting dose.
          ## This variable is updated after each cohort in the loop.
          stopit <- FALSE

          ## what is the next dose to be used?
          ## initialize with starting dose
          thisDose <- object@startingDose

          ## inside this loop we simulate the whole trial, until stopping
          while (!stopit) {
            ## what is the probability for tox. at this dose?
            thisDLEProb <- thisTruthDLE(thisDose)
            thisMeanEff <- thisTruthEff(thisDose)

            ## what is the cohort size at this dose?
            thisSize <- size(object@cohort_size,
              dose = thisDose,
              data = thisData
            )

            ## simulate DLTs: depends on whether we
            ## separate the first patient or not.
            if (firstSeparate && (thisSize > 1L)) {
              ## dose the first patient
              thisDLTs <- rbinom(
                n = 1L,
                size = 1L,
                prob = thisDLEProb
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisDLTs.PL <- rbinom(
                  n = 1L,
                  size = 1L,
                  prob = thisProb.PL
                )
              }

              thisEff <- rnorm(
                n = 1L,
                mean = thisMeanEff,
                sd = sqrt(trueSigma2)
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisEff.PL <- rnorm(
                  n = 1L,
                  mean = thisMeanEff.PL,
                  sd = sqrt(trueSigma2)
                )
              }

              ## if there is no DLT:
              if (thisDLTs == 0) {
                ## enroll the remaining patients
                thisDLTs <- c(
                  thisDLTs,
                  rbinom(
                    n = thisSize - 1L,
                    size = 1L,
                    prob = thisDLEProb
                  )
                )
                thisEff <- c(
                  thisEff,
                  rnorm(
                    n = thisSize - 1L,
                    mean = thisMeanEff,
                    sd = sqrt(trueSigma2)
                  )
                )


                if (thisData@placebo && (thisSize.PL > 0L)) {
                  thisDLTs.PL <- c(
                    thisDLTs.PL,
                    rbinom(
                      n = thisSize.PL,
                      size = 1L,
                      prob = thisProb.PL
                    )
                  )
                  thisEff.PL <- c(
                    thisMeanEff.PL,
                    rnorm(
                      n = thisSize.PL,
                      mean = thisMeanEff.PL,
                      sd = sqrt(trueSigma2)
                    )
                  )
                }
              }
            } else {
              ## we can directly dose all patients
              thisDLTs <- rbinom(
                n = thisSize,
                size = 1L,
                prob = thisDLEProb
              )
              thisEff <- rnorm(
                n = thisSize,
                mean = thisMeanEff,
                sd = sqrt(trueSigma2)
              )

              if (thisData@placebo && (thisSize.PL > 0L)) {
                thisDLTs.PL <- rbinom(
                  n = thisSize.PL,
                  size = 1L,
                  prob = thisProb.PL
                )
                thisEff.PL <- rnorm(
                  n = thisSize.PL,
                  mean = thisMeanEff.PL,
                  sd = sqrt(trueSigma2)
                )
              }
            }



            ## update the data with this placebo (if any) cohort and then with active dose
            if (thisData@placebo && (thisSize.PL > 0L)) {
              thisData <- update(
                object = thisData,
                x = object@data@doseGrid[1],
                y = thisDLTs.PL,
                w = thisEff.PL,
                check = FALSE
              )

              ## update the data with active dose
              thisData <- update(
                object = thisData,
                x = thisDose,
                y = thisDLTs,
                w = thisEff,
                new_cohort = FALSE
              )
            } else {
              ## update the data with this cohort
              thisData <- update(
                object = thisData,
                x = thisDose,
                y = thisDLTs,
                w = thisEff
              )
            }


            ## Update model estimate in DLE and Eff models
            thisDLEModel <- update(
              object = object@model,
              data = thisData
            )

            thisEffModel <- update(
              object = object@eff_model,
              data = thisData
            )

            thisNu <- thisEffModel@nu

            thisDLEsamples <- mcmc(
              data = thisData,
              model = thisDLEModel,
              options = mcmcOptions
            )

            thisEffsamples <- mcmc(
              data = thisData,
              model = thisEffModel,
              options = mcmcOptions
            )

            thisSigma2 <- if (thisEffModel@use_fixed) {
              1 / thisNu
            } else {
              1 / (as.numeric(thisNu["a"] / thisNu["b"]))
            }

            ## what is the dose limit?
            doselimit <- maxDose(object@increments, data = thisData)



            ## => what is the next best dose?
            next_bd <- nextBest(object@nextBest,
              doselimit = doselimit,
              samples = thisDLEsamples,
              model = thisDLEModel,
              data = thisData,
              model_eff = thisEffModel,
              samples_eff = thisEffsamples,
              in_sim = TRUE
            )

            thisDose <- next_bd$next_dose
            thisTDtargetDuringTrial <- next_bd$dose_target_drt
            thisTDtargetDuringTrialAtDoseGrid <- next_bd$next_dose_drt
            thisTDtargetEndOfTrial <- next_bd$dose_target_eot
            thisTDtargetEndOfTrialAtDoseGrid <- next_bd$next_dose_eot
            thisGstar <- next_bd$dose_max_gain
            thisGstarAtDoseGrid <- next_bd$next_dose_max_gain

            Recommend <- min(thisTDtargetEndOfTrialAtDoseGrid, thisGstarAtDoseGrid)

            ## Find the 95 % CI and its ratio (upper to the lower of this 95% CI of each of the estimates)
            thisCITDEOT <- list(lower = next_bd$ci_dose_target_eot[1], upper = next_bd$ci_dose_target_eot[2])
            thisratioTDEOT <- next_bd$ci_ratio_dose_target_eot

            thisCIGstar <- list(lower = next_bd$ci_dose_max_gain[1], upper = next_bd$ci_dose_max_gain[2])
            thisratioGstar <- next_bd$ci_ratio_dose_max_gain

            ## Find the optimal dose
            OptimalDose <- min(thisGstar, thisTDtargetEndOfTrial)

            if (OptimalDose == thisGstar) {
              thisratio <- thisratioGstar
              thisCI <- thisCIGstar
            } else {
              thisratio <- thisratioTDEOT
              thisCI <- thisCITDEOT
            }


            ## evaluate stopping rules
            stopit <- stopTrial(object@stopping,
              dose = thisDose,
              samples = thisDLEsamples,
              model = thisDLEModel,
              data = thisData,
              TDderive = object@nextBest@derive,
              Effmodel = thisEffModel,
              Effsamples = thisEffsamples,
              Gstarderive = object@nextBest@mg_derive
            )
            stopit_results <- h_unpack_stopit(stopit)
          }
          ## get the fit
          thisDLEFit <- fit(
            object = thisDLEsamples,
            model = thisDLEModel,
            data = thisData
          )

          thisEffFit <- fit(
            object = thisEffsamples,
            model = thisEffModel,
            data = thisData
          )


          ## return the results
          thisResult <- list(
            data = thisData,
            dose = thisDose,
            TDtargetDuringTrial = thisTDtargetDuringTrial,
            TDtargetDuringTrialAtDoseGrid = thisTDtargetDuringTrialAtDoseGrid,
            TDtargetEndOfTrial = thisTDtargetEndOfTrial,
            TDtargetEndOfTrialAtDoseGrid = thisTDtargetEndOfTrialAtDoseGrid,
            Gstar = thisGstar,
            GstarAtDoseGrid = thisGstarAtDoseGrid,
            Recommend = Recommend,
            OptimalDose = OptimalDose,
            OptimalDoseAtDoseGrid = Recommend,
            ratio = thisratio,
            CI = thisCI,
            ratioGstar = thisratioGstar,
            CIGstar = thisCIGstar,
            ratioTDEOT = thisratioTDEOT,
            CITDEOT = thisCITDEOT,
            fitDLE = subset(thisDLEFit,
              select =
                c(middle, lower, upper)
            ),
            fitEff = subset(thisEffFit,
              select =
                c(middle, lower, upper)
            ),
            sigma2est = thisSigma2,
            stop = attr(
              stopit,
              "message"
            ),
            report_results = stopit_results
          )

          return(thisResult)
        }


        resultList <- get_result_list(
          fun = runSim,
          nsim = nsim,
          vars =
            c(
              "simSeeds",
              "args",
              "nArgs",
              "firstSeparate",
              "trueDLE",
              "trueEff",
              "trueNu",
              "object"
            ),
          parallel = parallel,
          n_cores = nCores
        )


        ## put everything in the Simulations format:

        ## setup the list for the simulated data objects
        dataList <- lapply(resultList, "[[", "data")

        ## the vector of the final dose recommendations
        recommendedDoses <- as.numeric(sapply(resultList, "[[", "Recommend"))

        ## set up list for the final TD during Trial Estimate
        TDtargetDuringTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrial"))

        ## set up list for the final TD End of Trial Estimate
        TDtargetEndOfTrialList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrial"))

        ## set up list for the final TD during Trial estimate at dose Grid
        TDtargetDuringTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetDuringTrialAtDoseGrid"))

        ## set up list for the final TD End Of Trial estimate at dose Grid
        TDtargetEndOfTrialDoseGridList <- as.numeric(sapply(resultList, "[[", "TDtargetEndOfTrialAtDoseGrid"))

        ## set up list for the final Gstar estimates
        GstarList <- as.numeric(sapply(resultList, "[[", "Gstar"))

        ## set up list for the final Gstar estimates at dose grid
        GstarAtDoseGridList <- as.numeric(sapply(resultList, "[[", "GstarAtDoseGrid"))

        ## set up list for final optimal dose estimates
        OptimalDoseList <- as.numeric(sapply(resultList, "[[", "OptimalDose"))

        ## set up list for final optimal dose estimates at dose Grid
        OptimalDoseAtDoseGridList <- as.numeric(sapply(resultList, "[[", "Recommend"))

        ## Set up the list for the final 95% CI obtained
        CIList <- lapply(resultList, "[[", "CI")

        ## Set up the list for the final ratios obtained
        ratioList <- as.numeric(sapply(resultList, "[[", "ratio"))

        ## Set up the list for the final 95% CI of the TDtarget End Of Trial obtained
        CITDEOTList <- lapply(resultList, "[[", "CITDEOT")

        ## Set up the list for the final ratios of the TDtarget End Of Trial obtained
        ratioTDEOTList <- as.numeric(sapply(resultList, "[[", "ratioTDEOT"))

        ## Set up the list for the final 95% CI of the Gstar obtained
        CIGstarList <- lapply(resultList, "[[", "CIGstar")

        ## Set up the list for the final ratios of the Gstar obtained
        ratioGstarList <- as.numeric(sapply(resultList, "[[", "ratioGstar"))

        ## set up the list for the final fits for both DLE and efficacy
        fitDLEList <- lapply(resultList, "[[", "fitDLE")
        fitEffList <- lapply(resultList, "[[", "fitEff")
        ## the vector of the sigma2
        sigma2Estimates <- as.numeric(sapply(resultList, "[[", "sigma2est"))

        ## the reasons for stopping
        stopReasons <- lapply(resultList, "[[", "stop")

        # individual stopping rule results as matrix, labels as column names
        stop_results <- lapply(resultList, "[[", "report_results")
        stop_report <- as.matrix(do.call(rbind, stop_results))


        ## return the results in the Simulations class object
        ret <- PseudoDualSimulations(
          data = dataList,
          doses = recommendedDoses,
          final_td_target_during_trial_estimates = TDtargetDuringTrialList,
          final_td_target_end_of_trial_estimates = TDtargetEndOfTrialList,
          final_td_target_during_trial_at_dose_grid = TDtargetDuringTrialDoseGridList,
          final_td_target_end_of_trial_at_dose_grid = TDtargetEndOfTrialDoseGridList,
          final_cis = CIList,
          final_ratios = ratioList,
          final_gstar_estimates = GstarList,
          final_gstar_at_dose_grid = GstarAtDoseGridList,
          final_gstar_cis = CIGstarList,
          final_gstar_ratios = ratioGstarList,
          final_tdeot_cis = CITDEOTList,
          final_tdeot_ratios = ratioTDEOTList,
          final_optimal_dose = OptimalDoseList,
          final_optimal_dose_at_dose_grid = OptimalDoseAtDoseGridList,
          fit = fitDLEList,
          fit_eff = fitEffList,
          sigma2_est = sigma2Estimates,
          stop_reasons = stopReasons,
          stop_report = stop_report,
          seed = RNGstate
        )
        return(ret)
      }
    }
)

## --------------------------------------------------------------------------

##' Simulate outcomes from a time-to-DLT augmented CRM design (`DADesign`)
##'
##' @param object the \code{\linkS4class{DADesign}} object we want to simulate
##'   data from
##' @param nsim the number of simulations (default: 1)
##' @param seed see \code{\link{set_seed}}
##' @param truthTox a function which takes as input a dose (vector) and returns the
##'   true probability (vector) for toxicity and the time DLT occurs. Additional
##'   arguments can be supplied in \code{args}.
##' @param truthSurv a CDF which takes as input a time (vector) and returns
##'   the true cumulative probability (vector) that the DLT would occur conditioning on the patient
##'   has DLTs.
##' @param trueTmax (`number` or `NULL`)\cr the true maximum time at which DLTs can occur. Note that this must be larger thank `Tmax` from the `object`'s base data, which is the length of the DLT window, i.e. until which time DLTs are officially declared as such and used in the trial.
##' @param args data frame with arguments for the \code{truth} function. The
##'   column names correspond to the argument names, the rows to the values of the
##'   arguments. The rows are appropriately recycled in the \code{nsim}
##'   simulations. In order to produce outcomes from the posterior predictive
##'   distribution, e.g, pass an \code{object} that contains the data observed so
##'   far, \code{truth} contains the \code{prob} function from the model in
##'   \code{object}, and \code{args} contains posterior samples from the model.
##' @param firstSeparate enroll the first patient separately from the rest of
##'   the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
##'   in this patient.
##' @param deescalate deescalation when a DLT occurs in cohorts with lower dose
##'   level.
##' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
##'   giving the MCMC options for each evaluation in the trial. By default,
##'   the standard options are used.
##' @param DA document or rename this parameter to make it more meaningful
##' @param parallel should the simulation runs be parallelized across the
##'   clusters of the computer? (not default)
##' @param nCores how many cores should be used for parallel computing?
##' Defaults to the number of cores on the machine (maximum 5)
##' @param \dots not used
##' @param derive a named list of functions which derives statistics, based on the
##' vector of posterior MTD samples. Each list element must therefore accept
##' one and only one argument, which is a numeric vector, and return a number.
##'
##' @return an object of class \code{\linkS4class{Simulations}}
##'
##' @example examples/design-method-simulate-DADesign.R
##' @export
##' @keywords methods
setMethod("simulate",
  signature =
    signature(
      object = "DADesign",
      nsim = "ANY",
      seed = "ANY"
    ),
  def =
    function(object, nsim = 1L, seed = NULL,
             truthTox, truthSurv, trueTmax = NULL, args = NULL, firstSeparate = FALSE,
             deescalate = TRUE,
             mcmcOptions = McmcOptions(),
             DA = TRUE,
             parallel = FALSE, nCores = min(parallel::detectCores(), 5),
             derive = list(),
             ...) {
      ## checks and extracts
      assert_function(truthTox)
      assert_function(truthSurv)
      assert_flag(firstSeparate)
      assert_count(nsim, positive = TRUE)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      args <- as.data.frame(args)
      nArgs <- max(nrow(args), 1L)

      ## seed handling
      RNGstate <- set_seed(seed)

      ## from this,
      ## generate the individual seeds for the simulation runs
      simSeeds <- sample(x = seq_len(1e5), size = as.integer(nsim))

      ## Define functions which are useful in DLT Surv generation
      inverse <- function(f, lower = -100, upper = 100) {
        function(y) {
          uniroot((function(x) f(x) - y),
            lower = lower, upper = upper
          )[1]$root
        }
      }

      ## The DLT window length
      thisData <- object@data
      Tmax <- thisData@Tmax

      if (is.null(trueTmax)) {
        trueTmax <- Tmax
      } else if (trueTmax < Tmax) {
        warning("trueTmax < Tmax! trueTmax is set to Tmax")
        trueTmax <- Tmax
      }

      ## Calculate the inverse function of Surv to DLT CDF
      itruthSurv <- inverse(truthSurv, 0, trueTmax)

      ## generate random variable of Surv to DLT data; return Tmax when no
      ## DLT
      rtruthSurv <- function(DLT, Tmax_, itruthSurv = itruthSurv) {
        u_i <- rep(-100, length(DLT)) # remember to check this

        if (sum(DLT == 0) > 0) {
          u_i[DLT == 0] <- Tmax_
        }

        if (sum(DLT == 1) > 0) {
          u_i[DLT == 1] <- unlist(lapply(runif(sum(DLT == 1), 0, 1), itruthSurv))
        }

        # make sure that results are always positive, otherwise we get
        # problems below.
        u_i[u_i == 0] <- 0.5
        return(u_i)
      }

      # A function to return follow up fulfull yes (TRUE) vs no (FALSE);
      ready_to_open <- function(day, window, thisSurv) {
        size <- length(thisSurv)
        # the date the patient starts;
        start_time <- apply(rbind(thisSurv[-size], window$patientGap[-1]), 2, min)
        # the relative time for each patient on the specified "date";
        individule_check <- day - cumsum(c(0, start_time))
        # the minial number should be 0;
        individule_check[individule_check < 0] <- 0
        follow_up <- apply(rbind(thisSurv, individule_check), 2, min)
        return(all((follow_up - apply(rbind(window$patientFollow, thisSurv), 2, min)) >= 0) & (max(follow_up) >= min(window$patientFollowMin, max(thisSurv))))
      }

      ## assume we have surfficient patients, i.e. patient can be immediately enrolled
      ## once the trial accumulation is open. This function will tell you when to open
      ## the next cohort;
      # this function applys to all trials;
      nextOpen <- function(window, thisSurv) {
        size <- length(thisSurv)

        window$patientGap <- window$patientGap[1:size] ## if length(window$pt)>length(thisSurv), assume the first length(thisSurv) patients were enrolled;
        ## if the DLT happens before the end of DLT window, then the next
        ## cohort/enrollment of the next patient would happened earlier;
        start_time <- apply(rbind(thisSurv[-size], window$patientGap[-1]), 2, min)
        # duration of the cohort (all DLT windows finished);
        maxT <- max(thisSurv + cumsum(c(0, start_time)))

        meetrequire <- sapply(1:maxT, function(i) {
          ready_to_open(i, window, thisSurv)
        })
        if (sum(meetrequire) > 0) {
          # the earliest time that the require is met;
          time <- min(c(1:maxT)[meetrequire])
        } else {
          time <- maxT
        }

        return(time)
      }

      ## the function to produce the run a single simulation
      ## with index "iterSim"
      runSim <- function(iterSim) {
        ## set the seed for this run
        set.seed(simSeeds[iterSim])

        # check<<-simSeeds[iterSim]
        ## what is now the argument for the truth?
        ## (appropriately recycled)
        thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop = FALSE]

        ## so this truth is...
        thisTruth <- function(dose) {
          do.call(
            truthTox,
            ## First argument: the dose
            c(
              dose,
              ## Following arguments
              thisArgs
            )
          )
        }

        ## start the simulated data with the provided one
        thisData <- object@data

        # In case there are placebo
        if (thisData@placebo) {
          ## what is the probability for tox. at placebo?
          thisProb.PL <- thisTruth(object@data@doseGrid[1])
        }

        ## shall we stop the trial?
        ## First, we want to continue with the starting dose.
        ## This variable is updated after each cohort in the loop.
        stopit <- FALSE

        ## the time clock for the study, set to 0 when the first simulated
        ## patient initially dosed;
        trialtime <- 0

        # initiate the true DLT, true DLT Surv and the C1/D1 for each patients
        factDLTs <- thisData@y
        factSurv <- thisData@u
        factT0 <- thisData@t0

        ## what is the next dose to be used?
        ## initialize with starting dose
        thisDose <- object@startingDose

        ## inside this loop we simulate the whole trial, until stopping
        while (!stopit) {
          ## what is the probability for tox. at this dose?
          thisProb <- thisTruth(thisDose)

          ## what is the cohort size at this dose?
          thisSize <- size(object@cohort_size,
            dose = thisDose,
            data = thisData
          )

          thisSafetywindow <- windowLength(object@safetyWindow, thisSize)
          # better: add a checkpoint in safetywindow--dim(safetywindow$pt)==thisSize;

          ## In case there are placebo
          if (thisData@placebo) {
            thisSize.PL <- size(object@pl_cohort_size,
              dose = thisDose,
              data = thisData
            )
          }


          ## simulate DLTs: depends on whether we
          ## separate the first patient or not.
          ## amended on May 24: if any patient had DLT before the
          ## first patient finished a staggered window
          ## further enrollment will be stopped;

          if (firstSeparate && (thisSize > 1L)) {
            ## dose the first patient
            thisDLTs <- rbinom(
              n = 1L,
              size = 1L,
              prob = thisProb
            )

            if (thisData@placebo) {
              thisDLTs.PL <- rbinom(
                n = 1L,
                size = 1L,
                prob = thisProb.PL
              )
            }

            thisSurv <- ceiling(rtruthSurv(DLT = thisDLTs, Tmax_ = trueTmax, itruthSurv = itruthSurv))

            if (Tmax < trueTmax) {
              thisDLTs[thisDLTs == 1 & thisSurv > Tmax] <- 0

              thisSurv <- apply(rbind(thisSurv, rep(Tmax, length(thisSurv))), 2, min)
            }

            thisT0 <- trialtime

            ## if there is no DLT during Safety window:
            ## and no DLTs of previous patients-->


            # need to update the DataDA object
            tempData <- update(
              object = thisData,
              y = c(factDLTs, thisDLTs), #### the y will be updated according to u
              u = c(factSurv, thisSurv),
              t0 = c(factT0, thisT0),
              x = thisDose,
              trialtime = trialtime + thisSafetywindow$patientGap[2]
            ) #### the u will be updated over time

            temptime <- (tempData@u + tempData@t0)[tempData@y == 1 & tempData@x <= thisDose]


            # identify number of DLTs occurs during the thisSafetywindow$pt[2]
            # if(thisSurv>thisSafetywindow$pt[2])
            if (sum(temptime > trialtime) == 0) {
              ## enroll the remaining patients
              thisDLTs <- c(
                thisDLTs,
                rbinom(
                  n = thisSize - 1L,
                  size = 1L,
                  prob = thisProb
                )
              )

              thisSurv <- c(
                thisSurv,
                ceiling(rtruthSurv(thisDLTs[-1], trueTmax, itruthSurv = itruthSurv))
              )

              if (Tmax < trueTmax) {
                thisDLTs[thisDLTs == 1 & thisSurv > Tmax] <- 0

                thisSurv <- apply(rbind(thisSurv, rep(Tmax, length(thisSurv))), 2, min)
              }

              # in case any DLT happens before the end of the safety window;
              real_window <- apply(rbind(thisSurv[-thisSize], thisSafetywindow$patientGap[-1]), 2, min)


              thisT0 <- trialtime + c(0, cumsum(real_window))

              if (thisData@placebo && (thisSize.PL > 1L)) {
                thisDLTs.PL <- c(
                  thisDLTs.PL,
                  rbinom(
                    n = thisSize.PL - 1L,
                    size = 1L,
                    prob = thisProb.PL
                  )
                )
              }
            }

            rm(tempData)
            rm(temptime)
          } else {
            ## we can directly dose all patients
            thisDLTs <- rbinom(
              n = thisSize,
              size = 1L,
              prob = thisProb
            )

            thisSurv <- ceiling(rtruthSurv(thisDLTs, trueTmax, itruthSurv = itruthSurv))
            ## should return a vector with a same dimention as thisDLTs
            if (Tmax < trueTmax) {
              thisDLTs[thisDLTs == 1 & thisSurv > Tmax] <- 0

              thisSurv <- apply(rbind(thisSurv, rep(Tmax, length(thisSurv))), 2, min)
            }
            # in case any DLT happens before the end of the safety window;
            real_window <- apply(rbind(thisSurv[-thisSize], thisSafetywindow$patientGap[-1]), 2, min)

            thisT0 <- trialtime + c(0, cumsum(real_window))
            ## should return a vector with a same dimention as thisDLTs

            if (thisData@placebo) {
              thisDLTs.PL <- rbinom(
                n = thisSize.PL,
                size = 1L,
                prob = thisProb.PL
              )
            }
          }


          ## update the data with this placebo (if any)
          ## cohort and then with active dose
          if (thisData@placebo) {
            thisData <- update(
              object = thisData,
              x = object@data@doseGrid[1],
              y = thisDLTs.PL
            )

            ## update the data with active dose
            thisData <- update(
              object = thisData,
              x = thisDose,
              y = thisDLTs,
              new_cohort = FALSE
            )

            ## JZ: additional part for DADesign--when to start the next cohort
            trialtime <- trialtime + nextOpen(
              window = thisSafetywindow,
              thisSurv = thisSurv
            )
          } else {
            ## JZ: since the whole y and u column need update.
            ## factDLTs and factSuev get update and then calculate the y and u value in
            ## thisData object
            #
            #                                                   ## update the data with this cohort
            #                                                   thisData <- update(object=thisData,
            #                                                                      x=thisDose,  ####the x will be constantly updated according to u
            #                                                                      y=thisDLTs,
            #                                                                      u=thisSurv)  ####the u will be constantly updated

            factDLTs <- c(factDLTs, thisDLTs)

            factSurv <- c(factSurv, thisSurv) # better: check the data type of factSurv and thisSurv;

            factT0 <- c(factT0, thisT0)


            tempnext <- nextOpen(
              window = thisSafetywindow,
              thisSurv = thisSurv
            )

            ##### if there are DLTs, patients in the higher cohorts will be dosed a lower dose or discontinue.
            if (deescalate == TRUE) {
              newDLTid <- ((factSurv + factT0) > trialtime & (factSurv + factT0 - trialtime) <= tempnext & factDLTs == 1)

              newDLTnum <- c(1:length(factDLTs))[newDLTid]

              newDLTnum <- newDLTnum[newDLTnum <= (length(factDLTs) - length(thisDLTs))]

              # if(ifelse(sum(newDLTnum)==0,Inf,min(newDLTnum))<=(length(factDLTs)-length(thisDLTs))){
              if (length(newDLTnum) > 0) {
                for (DLT_loop in newDLTnum) {
                  newDLTtime <- (factSurv + factT0)[DLT_loop]

                  # identify higher dose--impacted patients:
                  deescalateID <- c(DLT_loop:length(factDLTs))[c(thisData@x, rep(thisDose, length(thisDLTs)))[DLT_loop:length(factDLTs)] > thisData@x[DLT_loop]]


                  ## DLT will be observed once the followup time >= the time to DLT
                  factDLTs[deescalateID] <- as.integer(factDLTs * (newDLTtime >= factT0 + factSurv))[deescalateID]

                  ## update DLT free survival time
                  factSurv[deescalateID] <- apply(rbind(factSurv, newDLTtime - factT0), 2, min)[deescalateID]
                }

                tempnext <- min(tempnext, max((factSurv + factT0)[(length(factDLTs) - length(thisDLTs) + 1):length(factDLTs)]) - trialtime)
              }
            }






            ## JZ: future work: additional part for DADesign--when to start the next cohort
            ## nextOpen can be modified to incorporate different patient enrollment rate;
            ## currently assume we have sufficient patients;
            ## If there is a gap between cohorts for cohort manager meeting, it can be
            ## added to here;

            trialtime <- trialtime + tempnext



            ## Update thisData
            ## according to what can be observed by the time when the next cohort open;


            thisData <- update(
              object = thisData,
              y = factDLTs, #### the y will be updated according to u
              u = factSurv,
              t0 = factT0,
              x = thisDose,
              trialtime = trialtime
            ) #### the u will be updated over time

            try(if (length(thisData@x) != length(thisData@u) || length(thisData@u) != length(thisData@y)) {
              stop("x,y,u dimention error")
            })
          }


          # testthisdata<<-thisData

          ## what is the dose limit? #JZ should
          ## still work for the DataDA object
          doselimit <- maxDose(object@increments,
            data = thisData
          )



          ## generate samples from the model
          if (DA == TRUE) {
            thisSamples <- mcmc(
              data = thisData,
              model = object@model,
              options = mcmcOptions
            )
          } else if (DA == FALSE) {
            temp_model <- LogisticLogNormal(
              mean = object@model@params@mean,
              cov = object@model@params@cov,
              ref_dose = object@model@refDose
            )

            trunk_Data <- Data(
              x = thisData@x, y = thisData@y,
              doseGrid = thisData@doseGrid,
              cohort = thisData@cohort,
              ID = thisData@ID
            )

            thisSamples <- mcmc(
              data = trunk_Data,
              model = temp_model,
              options = mcmcOptions
            )
          }

          ## => what is the next best dose?
          thisDose <- nextBest(object@nextBest,
            doselimit = doselimit,
            samples = thisSamples,
            model = object@model,
            data = thisData
          )$value

          ## evaluate stopping rules
          stopit <- stopTrial(object@stopping,
            dose = thisDose,
            samples = thisSamples,
            model = object@model,
            data = thisData
          )
          stopit_results <- h_unpack_stopit(stopit)
        }

        ## get the fit
        thisFit <- fit(
          object = thisSamples,
          model = object@model,
          data = thisData
        )

        # Get the MTD estimate from the samples.

        target_dose_samples <- dose(
          mean(object@nextBest@target),
          model = object@model,
          samples = thisSamples
        )

        # Create a function for additional statistical summary.
        additional_stats <- lapply(derive, function(f) f(target_dose_samples))

        ## return the results
        thisResult <-
          list(
            data = thisData,
            dose = thisDose,
            duration = trialtime,
            fit =
              subset(thisFit,
                select = c(middle, lower, upper)
              ),
            stop =
              attr(
                stopit,
                "message"
              ),
            report_results = stopit_results,
            additional_stats = additional_stats
          )
        return(thisResult)
      }

      resultList <- get_result_list(
        fun = runSim, ## remove
        nsim = nsim,
        vars =
          c(
            "simSeeds",
            "args",
            "nArgs",
            "firstSeparate",
            "truthTox",
            "truthSurv",
            "object",
            "mcmcOptions",
            "nextOpen",
            "ready_to_open"
          ),
        parallel = parallel,
        n_cores = nCores
      )

      ## put everything in the Simulations format:

      ## setup the list for the simulated data objects
      dataList <- lapply(resultList, "[[", "data")


      ## the vector of the final dose recommendations
      recommendedDoses <- as.numeric(sapply(resultList, "[[", "dose"))

      ## the vector of the final trial duration;
      trialduration <- as.numeric(sapply(resultList, "[[", "duration"))

      ## setup the list for the final fits
      fitList <- lapply(resultList, "[[", "fit")

      ## the reasons for stopping
      stopReasons <- lapply(resultList, "[[", "stop")

      # individual stopping rule results as matrix, labels as column names
      stop_results <- lapply(resultList, "[[", "report_results")
      stop_report <- as.matrix(do.call(rbind, stop_results))

      additional_stats <- lapply(resultList, "[[", "additional_stats")

      ## return the results in the Simulations class object
      ret <- DASimulations(
        data = dataList,
        doses = recommendedDoses,
        fit = fitList,
        trialduration = trialduration,
        stop_report = stop_report,
        stop_reasons = stopReasons,
        additional_stats = additional_stats,
        seed = RNGstate
      )


      return(ret)
    }
)

## --------------------------------------------------------------------------
# nolint end

# simulate ----

## DesignGrouped ----

#' Simulate Method for the [`DesignGrouped`] Class
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' A simulate method for [`DesignGrouped`] designs.
#'
#' @param object (`DesignGrouped`)\cr the design we want to simulate trials from.
#' @param nsim (`number`)\cr how many trials should be simulated.
#' @param seed (`RNGstate`)\cr generated with [set_seed()].
#' @param truth (`function`)\cr a function which takes as input a dose (vector) and
#'   returns the true probability (vector) for toxicity for the mono arm.
#'   Additional arguments can be supplied in `args`.
#' @param combo_truth (`function`)\cr same as `truth` but for the combo arm.
#' @param args (`data.frame`)\cr optional `data.frame` with arguments that work
#'   for both the `truth` and `combo_truth` functions. The column names correspond to
#'   the argument names, the rows to the values of the arguments. The rows are
#'   appropriately recycled in the `nsim` simulations.
#' @param firstSeparate (`flag`)\cr whether to enroll the first patient separately
#'   from the rest of the cohort and close the cohort in case a DLT occurs in this
#'   first patient.
#' @param mcmcOptions (`McmcOptions`)\cr MCMC options for each evaluation in the trial.
#' @param parallel (`flag`)\cr whether the simulation runs are parallelized across the
#'   cores of the computer.
#' @param nCores (`number`)\cr how many cores should be used for parallel computing.
#' @param ... not used.
#'
#' @return A list of `mono` and `combo` simulation results as [`Simulations`] objects.
#'
#' @aliases simulate-DesignGrouped
#' @export
#' @example examples/Design-method-simulate-DesignGrouped.R
#'
setMethod(
  "simulate",
  signature =
    signature(
      object = "DesignGrouped",
      nsim = "ANY",
      seed = "ANY"
    ),
  def =
    function(object,
             nsim = 1L,
             seed = NULL,
             truth,
             combo_truth,
             args = data.frame(),
             firstSeparate = FALSE,
             mcmcOptions = McmcOptions(),
             parallel = FALSE,
             nCores = min(parallelly::availableCores(), 5),
             ...) {
      nsim <- as.integer(nsim)
      assert_function(truth)
      assert_function(combo_truth)
      assert_data_frame(args)
      assert_count(nsim, positive = TRUE)
      assert_flag(firstSeparate)
      assert_flag(parallel)
      assert_count(nCores, positive = TRUE)

      n_args <- max(nrow(args), 1L)
      rng_state <- set_seed(seed)
      sim_seeds <- sample.int(n = 2147483647, size = nsim)

      run_sim <- function(iter_sim) {
        set.seed(sim_seeds[iter_sim])
        current <- list(mono = list(), combo = list())
        # Define true toxicity functions.
        current$args <- args[(iter_sim - 1) %% n_args + 1, , drop = FALSE]
        current$mono$truth <- function(dose) do.call(truth, c(dose, current$args))
        current$combo$truth <- function(dose) do.call(combo_truth, c(dose, current$args))
        # Start the simulated data with the provided one.
        current$mono$data <- object@mono@data
        current$combo$data <- object@combo@data
        # We are in the first cohort and continue for mono and combo.
        current$first <- TRUE
        current$mono$stop <- current$combo$stop <- FALSE

        # What are the next doses to be used? Initialize with starting doses.
        if (object@same_dose_for_all || (!object@first_cohort_mono_only && object@same_dose_for_start)) {
          current$mono$dose <- current$combo$dose <- min(object@mono@startingDose, object@combo@startingDose)
        } else {
          current$mono$dose <- object@mono@startingDose
          current$combo$dose <- object@combo@startingDose
        }

        # Inside this loop we simulate the whole trial, until stopping.
        while (!(current$mono$stop && current$combo$stop)) {
          if (!current$mono$stop) {
            cohort_size_mono <- size(
              object@mono@cohort_size,
              dose = current$mono$dose,
              data = current$mono$data
            )
            this_prob_mono <- current$mono$truth(current$mono$dose)
            current$mono$data <- current$mono$data %>%
              h_determine_dlts(
                dose = current$mono$dose,
                prob = this_prob_mono,
                cohort_size = cohort_size_mono,
                first_separate = firstSeparate
              )
          }
          if (!current$combo$stop && (!current$first || !object@first_cohort_mono_only)) {
            cohort_size_combo <- size(
              object@combo@cohort_size,
              dose = current$combo$dose,
              data = current$combo$data
            )
            this_prob_combo <- current$combo$truth(current$combo$dose)
            current$combo$data <- current$combo$data %>%
              h_determine_dlts(
                dose = current$combo$dose,
                prob = this_prob_combo,
                cohort_size = cohort_size_combo,
                first_separate = firstSeparate
              )
          }

          current$grouped <- h_group_data(current$mono$data, current$combo$data)
          current$samples <- mcmc(current$grouped, object@model, mcmcOptions)
          if (!current$mono$stop) {
            current$mono$limit <- maxDose(object@mono@increments, data = current$mono$data)
            current$mono$dose <- object@mono@nextBest %>%
              nextBest(current$mono$limit, current$samples, object@model, current$grouped, group = "mono")
            current$mono$dose <- current$mono$dose$value
          }
          if (!current$combo$stop && (!current$first || !object@first_cohort_mono_only)) {
            current$combo$limit <- if (is.na(current$mono$dose)) {
              0
            } else {
              maxDose(object@combo@increments, current$combo$data) %>%
                min(current$mono$dose, na.rm = TRUE)
            }
            current$combo$dose <- object@combo@nextBest %>%
              nextBest(current$combo$limit, current$samples, object@model, current$grouped, group = "combo")
            current$combo$dose <- current$combo$dose$value
            current$combo$stop <- object@combo@stopping %>%
              stopTrial(current$combo$dose, current$samples, object@model, current$combo$data, group = "combo")
            current$combo$results <- h_unpack_stopit(current$combo$stop)
          }
          if (!current$mono$stop) {
            current$mono$stop <- object@mono@stopping %>%
              stopTrial(
                current$mono$dose, current$samples, object@model, current$mono$data,
                group = "mono", external = current$combo$stop
              )
            current$mono$results <- h_unpack_stopit(current$mono$stop)
          }
          if (object@same_dose_for_all && !current$mono$stop && !current$combo$stop) {
            current$mono$dose <- current$combo$dose <- min(current$mono$dose, current$combo$dose)
          }
          if (current$first) {
            current$first <- FALSE
            if (object@first_cohort_mono_only && object@same_dose_for_start) {
              current$mono$dose <- current$combo$dose <- min(current$mono$dose, current$combo$dose)
            }
          }
        }
        current$mono$fit <- fit(current$samples, object@model, current$grouped, group = "mono")
        current$combo$fit <- fit(current$samples, object@model, current$grouped, group = "combo")
        lapply(
          X = current[c("mono", "combo")], FUN = with,
          list(
            data = data, dose = dose, fit = subset(fit, select = -dose),
            stop = attr(stop, "message"), results = results
          )
        )
      }
      vars_needed <- c("simSeeds", "args", "nArgs", "truth", "combo_truth", "firstSeparate", "object", "mcmcOptions")

      result_list <- get_result_list(run_sim, nsim, vars_needed, parallel, nCores)
      # Now we have a list with each element containing mono and combo. Reorder this a bit:
      result_list <- list(
        mono = lapply(result_list, "[[", "mono"),
        combo = lapply(result_list, "[[", "combo")
      )
      # Put everything in a list with both mono and combo Simulations:
      lapply(result_list, function(this_list) {
        data_list <- lapply(this_list, "[[", "data")
        recommended_doses <- as.numeric(sapply(this_list, "[[", "dose"))
        fit_list <- lapply(this_list, "[[", "fit")
        stop_reasons <- lapply(this_list, "[[", "stop")
        report_results <- lapply(this_list, "[[", "results")
        stop_report <- as.matrix(do.call(rbind, report_results))
        additional_stats <- lapply(this_list, "[[", "additional_stats")


        Simulations(
          data = data_list,
          doses = recommended_doses,
          fit = fit_list,
          stop_reasons = stop_reasons,
          stop_report = stop_report,
          additional_stats = additional_stats,
          seed = rng_state
        )
      })
    }
)

# tidy ----

## tidy-DualDesign ----

#' @rdname tidy
#' @aliases tidy-DualDesign
#' @example examples/Design-method-tidyDualDesign.R
#'
#' @export
setMethod(
  f = "tidy",
  signature = signature(x = "DualDesign"),
  definition = function(x, ...) {
    # Some Design objects have complex attributes whose structure is not supported.
    rv <- h_tidy_all_slots(x, attributes = FALSE) %>% h_tidy_class(x)
    if (length(rv) == 1) {
      rv[[names(rv)[1]]] %>% h_tidy_class(x)
    } else {
      rv
    }
  }
)
