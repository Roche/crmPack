##' @include helpers.R
##' @include Model-class.R
NULL

# nolint start

##' Convert prior quantiles (lower, median, upper) to logistic (log)
##' normal model
##'
##' This function uses generalized simulated annealing to optimize
##' a \code{\linkS4class{LogisticNormal}} model to be as close as possible
##' to the given prior quantiles.
##'
##' @param dosegrid the dose grid
##' @param refDose the reference dose
##' @param lower the lower quantiles
##' @param median the medians
##' @param upper the upper quantiles
##' @param level the credible level of the (lower, upper) intervals (default:
##' 0.95)
##' @param logNormal use the log-normal prior? (not default) otherwise, the
##' normal prior for the logistic regression coefficients is used
##' @param parstart starting values for the parameters. By default, these
##' are determined from the medians supplied.
##' @param parlower lower bounds on the parameters (intercept alpha and the
##' slope beta, the corresponding standard deviations and the correlation.)
##' @param parupper upper bounds on the parameters
##' @param seed seed for random number generation
##' @param verbose be verbose? (default)
##' @param control additional options for the optimisation routine, see
##' \code{\link[GenSA]{GenSA}} for more details
##' @return a list with the best approximating \code{model}
##' (\code{\linkS4class{LogisticNormal}} or
##' \code{\linkS4class{LogisticLogNormal}}), the resulting \code{quantiles}, the
##' \code{required} quantiles and the \code{distance} to the required quantiles,
##' as well as the final \code{parameters} (which could be used for running the
##' algorithm a second time)
##'
##' @importFrom GenSA GenSA
##' @importFrom mvtnorm rmvnorm
##' @export
##' @keywords programming
Quantiles2LogisticNormal <- function(dosegrid,
                                     refDose,
                                     lower,
                                     median,
                                     upper,
                                     level = 0.95,
                                     logNormal = FALSE,
                                     parstart = NULL,
                                     parlower = c(-10, -10, 0, 0, -0.95),
                                     parupper = c(10, 10, 10, 10, 0.95),
                                     seed = 12345,
                                     verbose = TRUE,
                                     control =
                                       list(
                                         threshold.stop = 0.01,
                                         maxit = 50000,
                                         temperature = 50000,
                                         max.time = 600
                                       )) {
  ## extracts and checks
  nDoses <- length(dosegrid)

  assert_flag(logNormal)
  assert_flag(verbose)
  assert_probability(level, bounds_closed = FALSE)
  stopifnot(
    !is.unsorted(dosegrid, strictly = TRUE),
    ## the medians must be monotonically increasing:
    !is.unsorted(median),
    identical(length(lower), nDoses),
    identical(length(median), nDoses),
    identical(length(upper), nDoses),
    all(lower < median),
    all(upper > median),
    identical(length(parlower), 5L),
    identical(length(parupper), 5L),
    all(parlower < parstart),
    all(parstart < parupper)
  )

  ## put verbose argument in the control list
  control$verbose <- verbose

  ## parametrize in terms of the means for the intercept alpha and the
  ## (log) slope beta,
  ## the corresponding standard deviations and the correlation.
  ## Define start values for optimisation:
  startValues <-
    if (is.null(parstart)) {
      ## find approximate means for alpha and slope beta
      ## from fitting logistic model to medians:
      startAlphaBeta <-
        coef(lm(I(logit(median)) ~ I(log(dosegrid / refDose))))

      ## overall starting values:
      c(
        meanAlpha =
          startAlphaBeta[1],
        meanBeta =
          if (logNormal) log(startAlphaBeta[2]) else startAlphaBeta[2],
        sdAlpha =
          1,
        sdBeta =
          1,
        correlation =
          0
      )
    } else {
      parstart
    }

  ## what is the target function which we want to minimize?
  target <- function(param) {
    ## form the mean vector and covariance matrix
    mean <- param[1:2]
    cov <- matrix(
      c(
        param[3]^2,
        prod(param[3:5]),
        prod(param[3:5]),
        param[4]^2
      ),
      nrow = 2L, ncol = 2L
    )

    ## simulate from the corresponding normal distribution
    set.seed(seed)
    normalSamples <- mvtnorm::rmvnorm(
      n = 1e4L,
      mean = mean,
      sigma = cov
    )

    ## extract separate coefficients
    alphaSamples <- normalSamples[, 1L]
    betaSamples <- if (logNormal) exp(normalSamples[, 2L]) else normalSamples[, 2L]

    ## and compute resulting quantiles
    quants <- matrix(
      nrow = length(dosegrid),
      ncol = 3L
    )
    colnames(quants) <- c("lower", "median", "upper")

    ## process each dose after another:
    for (i in seq_along(dosegrid))
    {
      ## create samples of the probability
      probSamples <-
        plogis(alphaSamples + betaSamples * log(dosegrid[i] / refDose))

      ## compute lower, median and upper quantile
      quants[i, ] <-
        quantile(probSamples,
          probs = c((1 - level) / 2, 0.5, (1 + level) / 2)
        )
    }

    ## now we can compute the target value
    ret <- max(abs(quants - c(lower, median, upper)))
    return(structure(ret,
      mean = mean,
      cov = cov,
      quantiles = quants
    ))
  }

  set.seed(seed)
  ## now optimise the target
  genSAres <- GenSA::GenSA(
    par = startValues,
    fn = target,
    lower = parlower,
    upper = parupper,
    control = control
  )
  distance <- genSAres$value
  pars <- genSAres$par
  targetRes <- target(pars)

  ## and construct the model
  ret <-
    if (logNormal) {
      LogisticLogNormal(
        mean = attr(targetRes, "mean"),
        cov = attr(targetRes, "cov"),
        ref_dose = refDose
      )
    } else {
      LogisticNormal(
        mean = attr(targetRes, "mean"),
        cov = attr(targetRes, "cov"),
        ref_dose = refDose
      )
    }

  ## return it together with the resulting distance and the quantiles
  return(list(
    model = ret,
    parameters = pars,
    quantiles = attr(targetRes, "quantiles"),
    required = cbind(lower, median, upper),
    distance = distance
  ))
}

# nolint end

#' Helper for Minimal Informative Unimodal Beta Distribution
#'
#' As defined in Neuenschwander et al (2008), this function computes the
#' parameters of the minimal informative unimodal beta distribution, given the
#' request that the p-quantile should be q, i.e. `X ~ Be(a, b)` with
#' `Pr(X <= q) = p`.
#'
#' @param p (`number`)\cr the probability.
#' @param q (`number`)\cr the quantile.
#' @return A list with the two resulting beta parameters `a` and `b`.
#'
#' @keywords internal
h_get_min_inf_beta <- function(p, q) {
  assert_probability(p, bounds_closed = FALSE)
  assert_probability(q, bounds_closed = FALSE)

  if (q > p) {
    list(
      a = log(p) / log(q),
      b = 1
    )
  } else {
    list(
      a = 1,
      b = log(1 - p) / log(1 - q)
    )
  }
}

# nolint start

##' Construct a minimally informative prior
##'
##' This function constructs a minimally informative prior, which is captured in
##' a \code{\linkS4class{LogisticNormal}} (or
##' \code{\linkS4class{LogisticLogNormal}}) object.
##'
##' Based on the proposal by Neuenschwander et al (2008, Statistics in
##' Medicine), a minimally informative prior distribution is constructed. The
##' required key input is the minimum (\eqn{d_{1}} in the notation of the
##' Appendix A.1 of that paper) and the maximum value (\eqn{d_{J}}) of the dose
##' grid supplied to this function. Then \code{threshmin} is the probability
##' threshold \eqn{q_{1}}, such that any probability of DLT larger than
##' \eqn{q_{1}} has only 5% probability. Therefore \eqn{q_{1}} is the 95%
##' quantile of the beta distribution and hence \eqn{p_{1} = 0.95}. Likewise,
##' \code{threshmax} is the probability threshold \eqn{q_{J}}, such that any
##' probability of DLT smaller than \eqn{q_{J}} has only 5% probability
##' (\eqn{p_{J} = 0.05}). The probabilities \eqn{1 - p_{1}} and \eqn{p_{J}} can be
##' controlled with the arguments \code{probmin} and \code{probmax}, respectively.
##' Subsequently, for all doses supplied in the
##' \code{dosegrid} argument, beta distributions are set up from the assumption
##' that the prior medians are linear in log-dose on the logit scale, and
##' \code{\link{Quantiles2LogisticNormal}} is used to transform the resulting
##' quantiles into an approximating \code{\linkS4class{LogisticNormal}} (or
##' \code{\linkS4class{LogisticLogNormal}}) model. Note that the reference dose
##' is not required for these computations.
##'
##' @param dosegrid the dose grid
##' @param refDose the reference dose
##' @param threshmin Any toxicity probability above this threshold would
##' be very unlikely (see \code{probmin}) at the minimum dose (default: 0.2)
##' @param threshmax Any toxicity probability below this threshold would
##' be very unlikely (see \code{probmax}) at the maximum dose (default: 0.3)
##' @param probmin the prior probability of exceeding \code{threshmin} at the
##' minimum dose (default: 0.05)
##' @param probmax the prior probability of being below \code{threshmax} at the
##' maximum dose (default: 0.05)
##' @param \dots additional arguments for computations, see
##' \code{\link{Quantiles2LogisticNormal}}, e.g. \code{refDose} and
##' \code{logNormal=TRUE} to obtain a minimal informative log normal prior.
##' @return see \code{\link{Quantiles2LogisticNormal}}
##'
##' @example examples/MinimalInformative.R
##' @export
##' @keywords programming
MinimalInformative <- function(dosegrid,
                               refDose,
                               threshmin = 0.2,
                               threshmax = 0.3,
                               probmin = 0.05,
                               probmax = 0.05,
                               ...) {
  ## extracts and checks
  nDoses <- length(dosegrid)

  assert_probability(threshmin, bounds_closed = FALSE)
  assert_probability(threshmax, bounds_closed = FALSE)
  assert_probability(probmin, bounds_closed = FALSE)
  assert_probability(probmax, bounds_closed = FALSE)
  stopifnot(
    !is.unsorted(dosegrid, strictly = TRUE)
  )
  xmin <- dosegrid[1]
  xmax <- dosegrid[nDoses]

  ## derive the beta distributions at the lowest and highest dose
  betaAtMin <- h_get_min_inf_beta(
    q = threshmin,
    p = 1 - probmin
  )
  betaAtMax <- h_get_min_inf_beta(
    q = threshmax,
    p = probmax
  )

  ## get the medians of those beta distributions
  medianMin <- with(
    betaAtMin,
    qbeta(p = 0.5, a, b)
  )
  medianMax <- with(
    betaAtMax,
    qbeta(p = 0.5, a, b)
  )

  ## now determine the medians of all beta distributions
  beta <- (logit(medianMax) - logit(medianMin)) / (log(xmax) - log(xmin))
  alpha <- logit(medianMax) - beta * log(xmax / refDose)
  medianDosegrid <- plogis(alpha + beta * log(dosegrid / refDose))

  ## finally for all doses calculate 95% credible interval bounds
  ## (lower and upper)
  lower <- upper <- dosegrid
  for (i in seq_along(dosegrid))
  {
    ## get min inf beta distribution
    thisMinBeta <- h_get_min_inf_beta(
      p = 0.5,
      q = medianDosegrid[i]
    )

    ## derive required quantiles
    lower[i] <- with(
      thisMinBeta,
      qbeta(p = 0.025, a, b)
    )
    upper[i] <- with(
      thisMinBeta,
      qbeta(p = 0.975, a, b)
    )
  }

  ## now go to Quantiles2LogisticNormal
  Quantiles2LogisticNormal(
    dosegrid = dosegrid,
    refDose = refDose,
    lower = lower,
    median = medianDosegrid,
    upper = upper,
    level = 0.95,
    ...
  )
}

# nolint end
