#' @include helpers.R
#' @include Data-class.R
#' @include Simulations-validity.R
#' @include CrmPackClass-class.R
NULL

# GeneralSimulations ----

## class ----

#' `GeneralSimulations` @description `r lifecycle::badge("stable")`
#' This class captures trial simulations.
#' Here also the random generator state before starting the simulation is
#' saved, in order to be able to reproduce the outcome. For this just use
#' [`set.seed`] with the `seed` as argument before running
#' [`simulate,Design-method`].
#'
#' @slot data (`list`)\cr produced [`Data`] objects.
#' @slot doses (`numeric`)\cr final dose recommendations.
#' @slot seed (`integer`)\cr random generator state before starting the simulation.
#'
#' @aliases GeneralSimulations
#' @export
.GeneralSimulations <-
  setClass(
    Class = "GeneralSimulations",
    slots = c(
      data = "list",
      doses = "numeric",
      seed = "integer"
    ),
    prototype = prototype(
      data =
        list(
          Data(
            x = 1:2,
            y = 0:1,
            doseGrid = 1:2,
            ID = 1L:2L,
            cohort = 1L:2L
          ),
          Data(
            x = 3:4,
            y = 0:1,
            doseGrid = 3:4,
            ID = 1L:2L,
            cohort = 1L:2L
          )
        ),
      doses = c(1, 2),
      seed = 1L
    ),
    contains = "CrmPackClass",
    validity = v_general_simulations
  )

## constructor ----

#' @rdname GeneralSimulations-class
#'
#' @param data (`list`)\cr see slot definition.
#' @param doses (`numeric`)\cr see slot definition.
#' @param seed (`integer`)\cr see slot definition.
#'
#' @example examples/Simulations-class-GeneralSimulations.R
#' @export
GeneralSimulations <- function(data,
                               doses,
                               seed) {
  assert_integerish(seed)
  .GeneralSimulations(
    data = data,
    doses = doses,
    seed = as.integer(seed)
  )
}


## default constructor

#' @rdname GeneralSimulations-class
#' @note Typically, end users will not use the `.DefaultGeneralSimulations()` function.
#' @export
.DefaultGeneralSimulations <- function() {
  GeneralSimulations(
    data = list(
      Data(x = 1:3, y = c(0, 1, 0), doseGrid = 1:3, ID = 1L:3L, cohort = 1L:3L),
      Data(x = 4:6, y = c(0, 1, 0), doseGrid = 4:6, ID = 1L:3L, cohort = 1L:3L)
    ),
    doses = c(1, 2),
    seed = 123
  )
}


# Simulations ----

## class ----

#' `Simulations`
#'
#' @description `r lifecycle::badge("stable")`
#'
#' This class captures the trial simulations from model based designs.
#' Additional slots `fit`, `stop_reasons`, `stop_report`,`additional_stats` compared to
#' the general class [`GeneralSimulations`].
#'
#' @slot fit (`list`)\cr final fits
#' @slot stop_reasons (`list`)\cr stopping reasons for each simulation run
#' @slot stop_report matrix of stopping rule outcomes
#' @slot additional_stats list of additional statistical summary
#' @aliases Simulations
#' @export
.Simulations <-
  setClass(
    Class = "Simulations",
    slots = c(
      fit = "list",
      stop_report = "matrix",
      stop_reasons = "list",
      additional_stats = "list"
    ),
    prototype = prototype(
      fit =
        list(
          c(0.1, 0.2),
          c(0.1, 0.2)
        ),
      stop_report = matrix(TRUE, nrow = 2),
      stop_reasons =
        list("A", "A"),
      additional_stats =
        list(a = 1, b = 1)
    ),
    contains = "GeneralSimulations",
    validity = v_simulations
  )

## constructor ----

#' @rdname Simulations-class
#'
#' @param fit (`list`)\cr see slot definition.
#' @param stop_reasons (`list`)\cr see slot definition.
#' @param stop_report see [`Simulations`]
#' @param additional_stats (`list`)\cr see slot definition.
#' @param \dots additional parameters from [`GeneralSimulations`]
#'
#' @example examples/Simulations-class-Simulations.R
#' @export
Simulations <- function(fit,
                        stop_reasons,
                        stop_report,
                        additional_stats,
                        ...) {
  start <- GeneralSimulations(...)
  .Simulations(start,
    fit = fit,
    stop_report = stop_report,
    stop_reasons = stop_reasons,
    additional_stats = additional_stats
  )
}

## default constructor ----

#' @rdname Simulations-class
#' @note Typically, end users will not use the `.DefaultSimulations()` function.
#' @export
.DefaultSimulations <- function() {
  design <- .DefaultDesign()
  myTruth <- probFunction(design@model, alpha0 = 7, alpha1 = 8)

  simulate(
    design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    mcmcOptions = .DefaultMcmcOptions(),
    parallel = FALSE
  )
}

# DualSimulations ----

## class ----

#' `DualSimulations`
#'
#' @description `r lifecycle::badge("stable")`
#'
#' This class captures the trial simulations from dual-endpoint model based
#' designs. In comparison to the parent class [`Simulations`],
#' it contains additional slots to capture the dose-biomarker `fits`, and the
#' `sigma2W` and `rho` estimates.
#'
#' @slot rho_est (`numeric`)\cr vector of final posterior median rho estimates
#' @slot sigma2w_est (`numeric`)\cr vector of final posterior median sigma2W estimates
#' @slot fit_biomarker (`list`)\cr with the final dose-biomarker curve fits
#' @aliases DualSimulations
#' @export
.DualSimulations <-
  setClass(
    Class = "DualSimulations",
    slots = c(
      rho_est = "numeric",
      sigma2w_est = "numeric",
      fit_biomarker = "list"
    ),
    prototype = prototype(
      rho_est = c(0.2, 0.3),
      sigma2w_est = c(0.2, 0.3),
      fit_biomarker =
        list(
          c(0.1, 0.2),
          c(0.1, 0.2)
        )
    ),
    contains = "Simulations",
    validity = v_dual_simulations
  )


## constructor ----

#' @rdname DualSimulations-class
#'
#' @param rho_est (`numeric`)\cr see [`DualSimulations`]
#' @param sigma2w_est (`numeric`)\cr [`DualSimulations`]
#' @param fit_biomarker (`list`)\cr see [`DualSimulations`]
#' @param \dots additional parameters from [`Simulations`]
#'
#' @example examples/Simulations-class-DualSimulations.R
#' @export
DualSimulations <- function(rho_est,
                            sigma2w_est,
                            fit_biomarker,
                            ...) {
  start <- Simulations(...)
  .DualSimulations(start,
    rho_est = rho_est,
    sigma2w_est = sigma2w_est,
    fit_biomarker = fit_biomarker
  )
}

## default constructor ----

#' @rdname DualSimulations-class
#' @note Typically, end users will not use the `.DefaultDualSimulations()` function.
#' @export
.DefaultDualSimulations <- function() {
  DualSimulations(
    rho_est = c(0.25, 0.35),
    sigma2w_est = c(0.15, 0.25),
    fit_biomarker = list(c(0.3, 0.4), c(0.4, 0.5)),
    fit = list(
      c(0.1, 0.2),
      c(0.3, 0.4)
    ),
    stop_report = matrix(c(TRUE, FALSE), nrow = 2),
    stop_reasons = list("A", "B"),
    additional_stats = list(a = 1, b = 1),
    data = list(
      Data(
        x = 1:2,
        y = 0:1,
        doseGrid = 1:2,
        ID = 1L:2L,
        cohort = 1L:2L
      ),
      Data(
        x = 3:4,
        y = 0:1,
        doseGrid = 3:4,
        ID = 1L:2L,
        cohort = 1L:2L
      )
    ),
    doses = c(1, 2),
    seed = 123L
  )
}

#' `GeneralSimulationsSummary`
#'
#' @description `r lifecycle::badge("stable")`
#'
#' This class captures the summary of general simulations output. Note that objects
#' should not be created by users, therefore no initialization
#' function is provided for this class.
#'
#' @slot target (`numeric`)\cr target toxicity interval
#' @slot target_dose_interval (`numeric`)\cr corresponding target dose interval
#' @slot nsim (`integer`)\cr number of simulations
#' @slot prop_dlts (`ANY`)\cr A numeric array (multi-dimensional) or list representing proportions of DLTs in the trials
#' @slot mean_tox_risk (`numeric`)\cr mean toxicity risks for the patients
#' @slot dose_selected (`numeric`)\cr doses selected as MTD
#' @slot tox_at_doses_selected (`numeric`)\cr true toxicity at doses selected
#' @slot prop_at_target (`numeric`)\cr Proportion of trials selecting target MTD
#' @slot dose_most_selected (`numeric`)\cr dose most often selected as MTD
#' @slot obs_tox_rate_at_dose_most_selected (`numeric`)\cr observed toxicity rate at dose most often selected
#' @slot n_obs (`ANY`)\cr A numeric array (multi-dimensional) or list representing number of patients overall.
#' @slot n_above_target (`integer`)\cr number of patients treated above target tox interval
#' @slot dose_grid (`numeric`)\cr the dose grid that has been used
#' @slot placebo (`logical`)\cr set to TRUE (default is FALSE) for a design with placebo
#' @aliases GeneralSimulationsSummary
#' @export
.GeneralSimulationsSummary <-
  setClass(
    Class = "GeneralSimulationsSummary",
    slots = c(
      target = "numeric",
      target_dose_interval = "numeric",
      nsim = "integer",
      prop_dlts = "ANY",
      mean_tox_risk = "numeric",
      dose_selected = "numeric",
      tox_at_doses_selected = "numeric",
      prop_at_target = "numeric",
      dose_most_selected = "numeric",
      obs_tox_rate_at_dose_most_selected = "numeric",
      n_obs = "ANY",
      n_above_target = "integer",
      dose_grid = "numeric",
      placebo = "logical"
    )
  )

## default constructor ----

#' @rdname GeneralSimulationsSummary-class
#' @note Typically, end users will not use the `.DefaultGeneralSimulationsSummary()` function.
#' @export
.DefaultGeneralSimulationsSummary <- function() {
  stop(
    paste(
      "Class GeneralSimulationsSummary cannot be instantiated directly.",
      "Please use one of its subclasses instead."
    )
  )
}

## SimulationsSummary ----

## class ----

#' `SimulationsSummary`
#'
#' @description `r lifecycle::badge("stable")`
#'
#' In addition to the slots in the parent class [`GeneralSimulationsSummary`],
#' it contains two slots with model fit information.
#'
#' @slot stop_report (`matrix`)\cr matrix of stopping rule outcomes
#' @slot fit_at_dose_most_selected (`numeric`)\cr fitted toxicity rate at dose most often selected
#' @slot additional_stats (`list`)\cr list of additional statistical summary
#' @slot mean_fit (`list`)\cr list with the average, lower (2.5%) and upper (97.5%)
#' quantiles of the mean fitted toxicity at each dose level
#'
#' @aliases SimulationsSummary
#' @export
.SimulationsSummary <-
  setClass(
    Class = "SimulationsSummary",
    slots = c(
      stop_report = "matrix",
      fit_at_dose_most_selected = "numeric",
      additional_stats = "list",
      mean_fit = "list"
    ),
    contains = "GeneralSimulationsSummary"
  )

## default constructor ----

#' @rdname SimulationsSummary-class
#' @note Typically, end users will not use the `.DefaultSimulationsSummary()` function.
#' @export
.DefaultSimulationsSummary <- function() {
  stop(paste(
    "Class SimulationsSummary cannot be instantiated directly.",
    "Please use one of its subclasses instead."
  ))
}

# DualSimulationsSummary ----

# class ----

#' `DualSimulationsSummary`
#'
#' @description `r lifecycle::badge("stable")`
#' This class captures the summary of dual-endpoint simulations output.
#' In comparison to its parent class [`SimulationsSummary`], it has additional slots.
#'
#' @slot biomarker_fit_at_dose_most_selected (`numeric`)\cr fitted biomarker level at most often selected dose.
#' @slot mean_biomarker_fit (`list`)\cr list with average, lower (2.5%) and upper (97.5%) quantiles of
#' mean fitted biomarker level at each dose
#' @aliases DualSimulationsSummary
#' @export
.DualSimulationsSummary <-
  setClass(
    Class = "DualSimulationsSummary",
    slots = c(
      biomarker_fit_at_dose_most_selected = "numeric",
      mean_biomarker_fit = "list"
    ),
    contains = "SimulationsSummary"
  )

# default constructor

#' @rdname DualSimulationsSummary-class
#' @note Typically, end users will not use the `.DefaultDualSimulationsSummary()` function.
#' @export
.DefaultDualSimulationsSummary <- function() {
  empty_data <- DataDual(doseGrid = c(1, 3, 5, 10, 15, 20, 25, 30))

  my_model <- DualEndpointRW(
    mean = c(0, 1),
    cov = matrix(c(1, 0, 0, 1), nrow = 2),
    sigma2betaW = 0.01,
    sigma2W = c(a = 0.1, b = 0.1),
    rho = c(a = 1, b = 1),
    rw1 = TRUE
  )

  my_next_best <- NextBestDualEndpoint(
    target = c(0.9, 1),
    overdose = c(0.35, 1),
    max_overdose_prob = 0.25
  )

  my_size1 <- CohortSizeRange(
    intervals = c(0, 30),
    cohort_size = c(1, 3)
  )
  my_size2 <- CohortSizeDLT(
    intervals = c(0, 1),
    cohort_size = c(1, 3)
  )
  my_size <- maxSize(my_size1, my_size2)

  my_stopping1 <- StoppingTargetBiomarker(
    target = c(0.9, 1),
    prob = 0.5
  )

  my_stopping <- my_stopping1 | StoppingMinPatients(10) | StoppingMissingDose()

  my_increments <- IncrementsRelative(
    intervals = c(0, 20),
    increments = c(1, 0.33)
  )

  my_design <- DualDesign(
    model = my_model,
    data = empty_data,
    nextBest = my_next_best,
    stopping = my_stopping,
    increments = my_increments,
    cohort_size = CohortSizeConst(3),
    startingDose = 3
  )

  beta_mod <- function(dose, e0, eMax, delta1, delta2, scal) {
    maxDens <- (delta1^delta1) * (delta2^delta2) / ((delta1 + delta2)^(delta1 + delta2))
    dose <- dose / scal
    e0 + eMax / maxDens * (dose^delta1) * (1 - dose)^delta2
  }

  true_biomarker <- function(dose) {
    beta_mod(dose, e0 = 0.2, eMax = 0.6, delta1 = 5, delta2 = 5 * 0.5 / 0.5, scal = 100)
  }

  true_tox <- function(dose) {
    pnorm((dose - 60) / 10)
  }

  x <- simulate(
    object = my_design,
    trueTox = true_tox,
    trueBiomarker = true_biomarker,
    sigma2W = 0.01,
    rho = 0,
    nsim = 1,
    parallel = FALSE,
    seed = 3,
    startingDose = 6,
    mcmcOptions = .DefaultMcmcOptions()
  )
}

# PseudoSimulations ----

## class ----

#' `PseudoSimulations`
#'
#' @description `r lifecycle::badge("stable")`
#' This class captures trial simulations from designs using pseudo model.
#' It has additional slots `fit` and `stop_reasons` compared to the
#' general class [`GeneralSimulations`].
#'
#' @slot fit (`list`)\cr final fit values.
#' @slot final_td_target_during_trial_estimates (`numeric`)\cr final estimates of the `td_target_during_trial`.
#' @slot final_td_target_end_of_trial_estimates (`numeric`)\cr final estimates of the `td_target_end_of_trial`.
#' @slot final_td_target_during_trial_at_dose_grid (`numeric`)
#'        \cr dose levels at dose grid closest below the final `td_target_during_trial` estimates.
#' @slot final_td_target_end_of_trial_at_dose_grid (`numeric`)
#'        \cr dose levels at dose grid closest below the final `td_target_end_of_trial` estimates.
#' @slot final_tdeot_cis (`list`)\cr 95% credibility intervals of the final estimates for `td_target_end_of_trial`.
#' @slot final_tdeot_ratios (`numeric`)\cr ratio of the upper to the lower 95%
#'        credibility intervals for `td_target_end_of_trial`.
#' @slot final_cis (`list`)\cr final 95% credibility intervals for `td_target_end_of_trial` estimates.
#' @slot final_ratios (`numeric`)\cr final ratios of the upper to the lower 95%
#'        credibility interval for `td_target_end_of_trial`.
#' @slot stop_report (`matrix`)\cr outcomes of stopping rules.
#' @slot stop_reasons (`list`)\cr reasons for stopping each simulation run.
#'
#' @aliases PseudoSimulations
#' @export
.PseudoSimulations <-
  setClass(
    Class = "PseudoSimulations",
    slots = c(
      fit = "list",
      final_td_target_during_trial_estimates = "numeric",
      final_td_target_end_of_trial_estimates = "numeric",
      final_td_target_during_trial_at_dose_grid = "numeric",
      final_td_target_end_of_trial_at_dose_grid = "numeric",
      final_tdeot_cis = "list",
      final_tdeot_ratios = "numeric",
      final_cis = "list",
      final_ratios = "numeric",
      stop_report = "matrix",
      stop_reasons = "list"
    ),
    prototype = prototype(
      final_td_target_during_trial_estimates = c(0.1, 0.1),
      final_td_target_end_of_trial_estimates = c(0.1, 0.1),
      final_td_target_during_trial_at_dose_grid = c(0.1, 0.1),
      final_td_target_end_of_trial_at_dose_grid = c(0.1, 0.1),
      final_tdeot_cis = list(c(0.1, 0.2), c(0.1, 0.2)),
      final_tdeot_ratios = c(0.1, 0.1),
      final_cis = list(c(0.1, 0.2), c(0.1, 0.2)),
      final_ratios = c(0.1, 0.1),
      stop_report = matrix(TRUE, nrow = 2),
      stop_reasons = list("A", "A")
    ),
    contains = "GeneralSimulations",
    validity = v_pseudo_simulations
  )

## constructor ----

#' @rdname PseudoSimulations-class
#'
#' @param fit (`list`)\cr see slot definition.
#' @param final_td_target_during_trial_estimates (`numeric`)\cr see slot definition.
#' @param final_td_target_end_of_trial_estimates (`numeric`)\cr see slot definition.
#' @param final_td_target_during_trial_at_dose_grid (`numeric`)\cr see slot definition.
#' @param final_td_target_end_of_trial_at_dose_grid (`numeric`)\cr see slot definition.
#' @param final_tdeot_cis (`list`)\cr see slot definition.
#' @param final_tdeot_ratios (`numeric`)\cr see slot definition.
#' @param final_cis (`list`)\cr see slot definition.
#' @param final_ratios (`numeric`)\cr see slot definition.
#' @param stop_report see [`PseudoSimulations`]
#' @param stop_reasons (`list`)\cr see slot definition.
#' @param \dots additional parameters from [`GeneralSimulations`]
#'
#' @export
PseudoSimulations <- function(fit,
                              final_td_target_during_trial_estimates,
                              final_td_target_end_of_trial_estimates,
                              final_td_target_during_trial_at_dose_grid,
                              final_td_target_end_of_trial_at_dose_grid,
                              final_tdeot_cis,
                              final_tdeot_ratios,
                              final_cis,
                              final_ratios,
                              stop_report,
                              stop_reasons,
                              ...) {
  start <- GeneralSimulations(...)
  .PseudoSimulations(start,
    fit = fit,
    final_td_target_during_trial_estimates = final_td_target_during_trial_estimates,
    final_td_target_end_of_trial_estimates = final_td_target_end_of_trial_estimates,
    final_td_target_during_trial_at_dose_grid = final_td_target_during_trial_at_dose_grid,
    final_td_target_end_of_trial_at_dose_grid = final_td_target_end_of_trial_at_dose_grid,
    final_tdeot_cis = final_tdeot_cis,
    final_tdeot_ratios = final_tdeot_ratios,
    final_cis = final_cis,
    final_ratios = final_ratios,
    stop_report = stop_report,
    stop_reasons = stop_reasons
  )
}

## default constructor ----

#' @rdname PseudoSimulations-class
#' @note Typically, end users will not use the `.DefaultPseudoSimulations()` function.
#' @export
.DefaultPseudoSimulations <- function() {
  stop("Class PseudoSimulations cannot be instantiated directly. Please use one of its subclasses instead.")
}

# PseudoDualSimulations ----

## class ----

#' `PseudoDualSimulations`
#'
#' @description `r lifecycle::badge("stable")`
#' This class conducts trial simulations for designs using both the
#' DLE and efficacy responses. It defines final values for
#' efficacy fit and DLE, estimates of Gstar, optimal dose and sigma2.
#'
#' @slot fit_eff (`list`)\cr final values of efficacy fit.
#' @slot final_gstar_estimates (`numeric`)\cr final Gstar estimates.
#' @slot final_gstar_at_dose_grid (`numeric`)\cr final Gstar estimates at dose grid.
#' @slot final_gstar_cis (`list`)\cr list of 95% confidence interval for Gstar estimates.
#' @slot final_gstar_ratios (`numeric`)\cr ratios of confidence intervals for Gstar estimates.
#' @slot final_optimal_dose (`numeric`)\cr final optimal dose.
#' @slot final_optimal_dose_at_dose_grid (`numeric`)\cr final optimal dose at dose grid.
#' @slot sigma2_est (`numeric`)\cr final sigma2 estimates.
#'
#' @aliases PseudoDualSimulations
#' @export
.PseudoDualSimulations <-
  setClass(
    Class = "PseudoDualSimulations",
    slots = c(
      fit_eff = "list",
      final_gstar_estimates = "numeric",
      final_gstar_at_dose_grid = "numeric",
      final_gstar_cis = "list",
      final_gstar_ratios = "numeric",
      final_optimal_dose = "numeric",
      final_optimal_dose_at_dose_grid = "numeric",
      sigma2_est = "numeric"
    ),
    prototype = prototype(
      final_gstar_estimates = c(0.1, 0.1),
      final_gstar_at_dose_grid = c(0.1, 0.1),
      final_gstar_cis = list(c(0.1, 0.2), c(0.1, 0.2)),
      final_gstar_ratios = c(0.01, 0.01),
      final_optimal_dose = c(0.01, 0.01),
      final_optimal_dose_at_dose_grid = c(0.01, 0.01),
      sigma2_est = c(0.001, 0.002)
    ),
    contains = "PseudoSimulations",
    validity = v_pseudo_dual_simulations
  )

## constructor ----

#' @rdname PseudoDualSimulations-class
#'
#' @param fit_eff (`list`)\cr see slot definition.
#' @param final_gstar_estimates (`numeric`)\cr see slot definition.
#' @param final_gstar_at_dose_grid (`numeric`)\cr see slot definition.
#' @param final_gstar_cis (`list`)\cr see slot definition.
#' @param final_gstar_ratios (`numeric`)\cr see slot definition.
#' @param final_optimal_dose (`numeric`)\cr see slot definition.
#' @param final_optimal_dose_at_dose_grid (`numeric`)\cr see slot definition.
#' @param sigma2_est (`numeric`)\cr see slot definition.
#' @param \dots additional parameters from [`PseudoSimulations`]
#' @export
PseudoDualSimulations <- function(fit_eff,
                                  final_gstar_estimates,
                                  final_gstar_at_dose_grid,
                                  final_gstar_cis,
                                  final_gstar_ratios,
                                  final_optimal_dose,
                                  final_optimal_dose_at_dose_grid,
                                  sigma2_est,
                                  ...) {
  start <- PseudoSimulations(...)
  .PseudoDualSimulations(start,
    fit_eff = fit_eff,
    final_gstar_estimates = final_gstar_estimates,
    final_gstar_at_dose_grid = final_gstar_at_dose_grid,
    final_gstar_cis = final_gstar_cis,
    final_gstar_ratios = final_gstar_ratios,
    final_optimal_dose = final_optimal_dose,
    final_optimal_dose_at_dose_grid = final_optimal_dose_at_dose_grid,
    sigma2_est = sigma2_est
  )
}

## default constructor ----

#' @rdname PseudoDualSimulations-class
#' @note Do not use the `.DefaultPseudoDualSimulations()` function.
#' @export
.DefaultPseudoDualSimulations <- function() {
  stop("Class PseudoDualSimulations cannot be instantiated directly. Please use a subclass.")
}

# PseudoDualFlexiSimulations ----

## class ----

#' `PseudoDualFlexiSimulations`
#'
#' @description `r lifecycle::badge("stable")`
#' This class captures the trial simulations design using both the DLE and
#' efficacy responses using [`EffFlexi`] efficacy model.
#' It extends [`PseudoDualSimulations`] by adding the capability to capture the sigma2betaW estimates.
#'
#' @slot sigma2_beta_west (`numeric`)\cr the vector of the final posterior mean sigma2betaW estimates
#' @aliases PseudoDualFlexiSimulations
#' @export
.PseudoDualFlexiSimulations <-
  setClass(
    Class = "PseudoDualFlexiSimulations",
    slots = c(sigma2_beta_west = "numeric"),
    prototype = prototype(sigma2_beta_west = c(0.001, 0.002)),
    contains = "PseudoDualSimulations"
  )

## constructor ----

#' @rdname PseudoDualFlexiSimulations-class
#'
#' @param sigma2_beta_west (`numeric`)\cr the vector of the final posterior mean sigma2betaW estimates
#' @param \dots additional parameters from [`PseudoDualSimulations`]
#'
#' @export
PseudoDualFlexiSimulations <- function(sigma2_beta_west,
                                       ...) {
  start <- PseudoDualSimulations(...)
  .PseudoDualFlexiSimulations(start,
    sigma2_beta_west = sigma2_beta_west
  )
}

## default constructor ----

#' @rdname PseudoDualFlexiSimulations-class
#' @note Typically, end users will not use the `.DefaultPseudoFlexiSimulations()` function.
#' @export
.DefaultPseudoDualFlexiSimulations <- function() {
  stop("Class PseudoFlexiSimulations cannot be instantiated directly. Please use one of its subclasses instead.")
}

# nolint start
## -------------------------------------------------------------------------------------------------------
## ================================================================================================

##' Class for the summary of pseudo-models simulations output
##'
##' Note that objects should not be created by users, therefore no
##' initialization function is provided for this class.
##'
##' @slot targetEndOfTrial the target probability of DLE wanted at the end of a trial
##' @slot targetDoseEndOfTrial the dose level corresponds to the target probability
##' of DLE wanted at the end of a trial, TDEOT
##' @slot targetDoseEndOfTrialAtDoseGrid the dose level at dose grid corresponds to the target probability
##' of DLE wanted at the end of a trial
##' @slot targetDuringTrial the target probability of DLE wanted during a trial
##' @slot targetDoseDuringTrial the dose level corresponds to the target probability of DLE
##' wanted during the trial. TDDT
##' @slot targetDoseDuringTrialAtDoseGrid the dose level at dose grid corresponds to the target probability
##' of DLE wanted during a trial
##' @slot TDEOTSummary the six-number table summary, include the lowest, the 25th precentile (lower quartile),
##' the 50th percentile (median), the mean, the 27th percentile and the highest values of the
##' final dose levels obtained corresponds to the target probability of DLE
##' want at the end of a trial across all simulations
##' @slot TDDTSummary the six-number table summary, include the lowest, the 25th precentile (lower quartile),
##' the 50th percentile (median), the mean, the 27th percentile and the highest values of the
##' final dose levels obtained corresponds to the target probability of DLE
##' want during a trial across all simulations
##' @slot FinalDoseRecSummary the six-number table summary, include the lowest, the 25th precentile (lower quartile),
##' the 50th percentile (median), the mean, the 27th percentile and the highest values of the
##' final optimal doses, which is either the TDEOT when only DLE response are incorporated into
##' the escalation procedure or the minimum of the TDEOT and Gstar when DLE and efficacy responses are
##' incorporated, across all simulations
##' @slot ratioTDEOTSummary the six-number summary table of the final ratios of the upper to the lower 95%
##' credibility intervals of the final TDEOTs across all simulations
##' @slot FinalRatioSummary the six-number summary table of the final ratios of the upper to the lower 95%
##' credibility intervals of the final optimal doses across all simulations
##' #@slot doseRec the dose level that will be recommend for subsequent study
##' @slot nsim number of simulations
##' @slot propDLE proportions of DLE in the trials
##' @slot meanToxRisk mean toxicity risks for the patients
##' @slot doseSelected doses selected as MTD (targetDoseEndOfTrial)
##' @slot toxAtDosesSelected true toxicity at doses selected
##' @slot propAtTargetEndOfTrial Proportion of trials selecting at the doseGrid closest below the MTD, the
##' targetDoseEndOfTrial
##' @slot propAtTargetDuringTrial Proportion of trials selecting at the doseGrid closest below the
##' targetDoseDuringTrial
##' @slot doseMostSelected dose most often selected as MTD
##' @slot obsToxRateAtDoseMostSelected observed toxicity rate at dose most often
##' selected
##' @slot nObs number of patients overall
##' @slot nAboveTargetEndOfTrial number of patients treated above targetDoseEndOfTrial
##' @slot nAboveTargetDuringTrial number of patients treated above targetDoseDuringTrial
##' @slot doseGrid the dose grid that has been used
##' @slot fitAtDoseMostSelected fitted toxicity rate at dose most often selected
##' @slot meanFit list with the average, lower (2.5%) and upper (97.5%)
##' quantiles of the mean fitted toxicity at each dose level
##' @slot stop_report matrix of stopping rule outcomes
##'
##' @export
##' @keywords classes
.PseudoSimulationsSummary <-
  setClass(
    Class = "PseudoSimulationsSummary",
    representation(
      targetEndOfTrial = "numeric",
      targetDoseEndOfTrial = "numeric",
      targetDoseEndOfTrialAtDoseGrid = "numeric",
      targetDuringTrial = "numeric",
      targetDoseDuringTrial = "numeric",
      targetDoseDuringTrialAtDoseGrid = "numeric",
      TDEOTSummary = "table",
      TDDTSummary = "table",
      FinalDoseRecSummary = "table",
      ratioTDEOTSummary = "table",
      FinalRatioSummary = "table",
      # doseRec="numeric",
      nsim = "integer",
      propDLE = "numeric",
      meanToxRisk = "numeric",
      doseSelected = "numeric",
      toxAtDosesSelected = "numeric",
      propAtTargetEndOfTrial = "numeric",
      propAtTargetDuringTrial = "numeric",
      doseMostSelected = "numeric",
      obsToxRateAtDoseMostSelected = "numeric",
      nObs = "integer",
      nAboveTargetEndOfTrial = "integer",
      nAboveTargetDuringTrial = "integer",
      doseGrid = "numeric",
      fitAtDoseMostSelected = "numeric",
      meanFit = "list",
      stop_report = "matrix"
    )
  )

## default constructor ----

#' @rdname GeneralSimulationsSummary-class
#' @note Typically, end users will not use the `.DefaultPseudoSimulationsSummary()` function.
#' @export
.DefaultPseudoSimulationsSummary <- function() {
  stop(paste0("Class PseudoSimulationsSummary cannot be instantiated directly.  Please use one of its subclasses instead."))
}

## ---------------------------------------------------------------------------------------------
##' Class for the summary of the dual responses simulations using pseudo models
##'
##' It contains all slots from \code{\linkS4class{PseudoSimulationsSummary}} object. In addition to
##' the slots in the parent class \code{\linkS4class{PseudoSimulationsSummary}}, it contains four
##' more slots for the efficacy model fit information.
##'
##' Note that objects should not be created by users, therefore no initialization function
##' is provided for this class.
##'
##' @slot targetGstar the target dose level such that its gain value is at maximum
##' @slot targetGstarAtDoseGrid the dose level at dose Grid closest and below Gstar
##' @slot GstarSummary the six-number table summary (lowest, 25th, 50th (median), 75th percentile, mean
##' and highest value) of the final Gstar values obtained across all simulations
##' @slot ratioGstarSummary the six-number summary table of the ratios of the upper to the lower 95%
##' credibility intervals of the final Gstar across all simulations
##' @slot EffFitAtDoseMostSelected fitted expected mean efficacy value at dose most often
##' selected
##' @slot meanEffFit list with mean, lower (2.5%) and upper (97.5%) quantiles of the fitted expected
##' efficacy value at each dose level.
##'
##' @export
##' @keywords class
.PseudoDualSimulationsSummary <-
  setClass(
    Class = "PseudoDualSimulationsSummary",
    contains = "PseudoSimulationsSummary",
    representation =
      representation(
        targetGstar = "numeric",
        targetGstarAtDoseGrid = "numeric",
        GstarSummary = "table",
        ratioGstarSummary = "table",
        EffFitAtDoseMostSelected = "numeric",
        meanEffFit = "list"
      )
  )

## default constructor ----

#' @rdname PseudoDualSimulationsSummary-class
#' @note Typically, end users will not use the `.DefaultPseudoDualSimulationsSummary()` function.
#' @export
.DefaultPseudoDualSimulationsSummary <- function() {
  stop(paste0("Class PseudoDualSimulationsSummary cannot be instantiated directly.  Please use one of its subclasses instead."))
}

## ---------------------------------------------------------------------------------------------

##' Class for the simulations output from DA based designs
##'
##' This class captures the trial simulations from DA based
##' designs. In comparison to the parent class \code{\linkS4class{Simulations}},
##' it contains additional slots to capture the time to DLT fits, additional
##' parameters and the trial duration.
##'
##' @slot trialduration the vector of trial duration values for all simulations.
##'
##' @export
##' @keywords classes
.DASimulations <-
  setClass(
    Class = "DASimulations",
    representation(trialduration = "numeric"),
    prototype(trialduration = rep(0, 2)),
    contains = "Simulations",
    validity = v_da_simulations
  )
validObject(.DASimulations())


##' Initialization function for `DASimulations`
##'
##' @param trialduration see \code{\linkS4class{DASimulations}}
##' @param \dots additional parameters from \code{\link{Simulations}}
##' @return the \code{\linkS4class{DASimulations}} object
##'
##' @export
##' @keywords methods
DASimulations <- function(trialduration,
                          ...) {
  start <- Simulations(...)
  .DASimulations(start,
    trialduration = trialduration
  )
}


## default constructor ----

#' @rdname DASimulations-class
#' @note Typically, end users will not use the `.DASimulations()` function.  This
#' function has a noticeable execution time.
#' @export
.DefaultDASimulations <- function() {
  design <- .DefaultDADesign()
  myTruth <- probFunction(design@model, alpha0 = 2, alpha1 = 3)
  exp_cond.cdf <- function(x, onset = 15) {
    a <- stats::pexp(28, 1 / onset, lower.tail = FALSE)
    1 - (stats::pexp(x, 1 / onset, lower.tail = FALSE) - a) / (1 - a)
  }

  simulate(
    design,
    args = NULL,
    truthTox = myTruth,
    truthSurv = exp_cond.cdf,
    trueTmax = 80,
    nsim = 2,
    seed = 819,
    mcmcOptions = .DefaultMcmcOptions(),
    firstSeparate = TRUE,
    deescalate = FALSE,
    parallel = FALSE
  )
}
# nolint end

# tidy

## tidy-Simulations ----

#' @rdname tidy
#' @aliases tidy-Simulations
#' @example examples/Simulations-method-tidy.R
#' @export
setMethod(
  f = "tidy",
  signature = signature(x = "Simulations"),
  definition = function(x, ...) {
    slot_names <- slotNames(x)
    rv <- list()
    for (nm in slot_names) {
      if (!is.function(slot(x, nm))) {
        if (nm %in% c("stop_reasons", "additional_stats")) {
        } else {
          rv[[nm]] <- h_tidy_slot(x, nm)
        }
      }
    }
    # Column bind of all list elements have the same number of rows
    if (length(rv) > 1 & length(unique(sapply(rv, nrow))) == 1) {
      rv <- rv %>% dplyr::bind_cols()
    }
    rv %>% h_tidy_class(x)
  }
)
