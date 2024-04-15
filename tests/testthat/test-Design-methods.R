options(testthat.progress.max_fails = 0)

#' @include("helper-design.R")

# simulate ----

## Design ----

test_that("simulate produces consistent results with placebo data", {
  design <- h_get_design_data(TRUE)
  myTruth <- probFunction(design@model, alpha0 = 7, alpha1 = 8)
  options <- h_get_mcmc_options()

  result <- simulate(
    design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    mcmcOptions = options,
    parallel = FALSE
  )

  expect_snapshot(result)
})

test_that("simulate produces consistent results with sentinel patients", {
  design <- h_get_design_data()
  myTruth <- probFunction(design@model, alpha0 = 7, alpha1 = 8)
  options <- h_get_mcmc_options()

  result <- simulate(
    design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    mcmcOptions = options,
    parallel = FALSE,
    firstSeparate = TRUE
  )

  expect_snapshot(result)
})

test_that("simulate for the class design returns correct objects", {
  design <- h_get_design_data()
  myTruth <- probFunction(design@model, alpha0 = 7, alpha1 = 8)
  options <- h_get_mcmc_options()

  mySims <- simulate(
    design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    mcmcOptions = options,
    parallel = FALSE,
    derive = list(
      max_mtd = max,
      mean_mtd = mean,
      median_mtd = median
    )
  )

  expect_class(mySims, "Simulations") # check for correct class of returned object

  expect_equal(any(sapply(mySims@fit[[1]], is.numeric)), TRUE) # check if all elements in mySims@fit are numeric

  expect_equal(length(mySims@stop_report), 5) # check for length

  expect_logical(mySims@stop_report) # check for stop_report to be logical vector

  expect_list(mySims@data)

  expect_class(mySims@data[[1]], "Data") # check for data object has correct class

  expect_list(mySims@additional_stats)

  expect_list(mySims@additional_stats[[1]])

  expect_length(mySims@additional_stats[[1]], 3)

  expect_equal(mySims@doses, 1)
})

test_that("simulate for the class design with placebo returns correct objects", {
  design <- h_get_design_data(TRUE)
  myTruth <- probFunction(design@model, alpha0 = 7, alpha1 = 8)
  options <- h_get_mcmc_options()

  mySims <- simulate(
    design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    mcmcOptions = options,
    parallel = FALSE,
    derive = list(
      max_mtd = max,
      mean_mtd = mean,
      median_mtd = median
    )
  )

  expect_class(mySims, "Simulations") # check for correct class of returned object

  expect_equal(any(sapply(mySims@fit[[1]], is.numeric)), TRUE) # check if all elements in mySims@fit are numeric

  expect_equal(length(mySims@stop_report), 5) # check for length

  expect_logical(mySims@stop_report) # check for stop_report to be logical vector

  expect_list(mySims@data)

  expect_class(mySims@data[[1]], "Data") # check for data object has correct class

  expect_list(mySims@additional_stats)

  expect_list(mySims@additional_stats[[1]])

  expect_length(mySims@additional_stats[[1]], 3)

  expect_equal(mySims@doses, 1)
})

test_that("simulate for the class design with placebo and sentinel patients returns correct objects", {
  design <- h_get_design_data(TRUE)
  myTruth <- probFunction(design@model, alpha0 = 7, alpha1 = 8)
  options <- h_get_mcmc_options()

  mySims <- simulate(
    design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    mcmcOptions = options,
    parallel = FALSE,
    firstSeparate = TRUE,
    derive = list(
      max_mtd = max,
      mean_mtd = mean,
      median_mtd = median
    )
  )

  expect_class(mySims, "Simulations") # check for correct class of returned object

  expect_equal(any(sapply(mySims@fit[[1]], is.numeric)), TRUE) # check if all elements in mySims@fit are numeric

  expect_equal(length(mySims@stop_report), 5) # check for length

  expect_logical(mySims@stop_report) # check for stop_report to be logical vector

  expect_list(mySims@data)

  expect_class(mySims@data[[1]], "Data") # check for data object has correct class

  expect_list(mySims@additional_stats)

  expect_list(mySims@additional_stats[[1]])

  expect_length(mySims@additional_stats[[1]], 3)

  expect_equal(mySims@doses, 1)
})

## RuleDesign ----

test_that("simulate-RuleDesign produces consistent results", {
  design <- ThreePlusThreeDesign(doseGrid = c(1, 3, 5, 10, 15, 20, 25, 40, 50, 80, 100))
  myTruth <- function(x) seq(0.05, 0.55, length.out = length(design@data@doseGrid))
  options <- h_get_mcmc_options()

  result <- simulate(
    design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    mcmcOptions = options,
    parallel = FALSE,
    firstSeparate = TRUE
  )

  expect_snapshot(result)
})

## DualDesign ----

test_that("simulate-DualDesign produces consistent results", {
  design <- h_get_design_dualdata()

  # define scenarios for the TRUE toxicity and efficacy profiles
  betaMod <- function(dose, e0, eMax, delta1, delta2, scal) {
    maxDens <- (delta1^delta1) * (delta2^delta2) / ((delta1 + delta2)^(delta1 + delta2))
    dose <- dose / scal
    e0 + eMax / maxDens * (dose^delta1) * (1 - dose)^delta2
  }

  trueBiomarker <- function(dose) {
    betaMod(dose, e0 = 0.2, eMax = 0.6, delta1 = 5, delta2 = 5 * 0.5 / 0.5, scal = 100)
  }

  trueTox <- function(dose) {
    pnorm((dose - 60) / 10)
  }

  result <- simulate(
    design,
    trueTox = trueTox,
    trueBiomarker = trueBiomarker,
    sigma2W = 0.01,
    rho = 0,
    nsim = 1,
    parallel = FALSE,
    seed = 3,
    startingDose = 6,
    mcmcOptions = McmcOptions(
      burnin = 100,
      step = 1,
      samples = 300,
      rng_kind = "Mersenne-Twister",
      rng_seed = 1234
    )
  )

  expect_snapshot(result)
})

test_that("simulate-DualDesign produces consistent results with sentinel patients", {
  design <- h_get_design_dualdata()

  # define scenarios for the TRUE toxicity and efficacy profiles
  betaMod <- function(dose, e0, eMax, delta1, delta2, scal) {
    maxDens <- (delta1^delta1) * (delta2^delta2) / ((delta1 + delta2)^(delta1 + delta2))
    dose <- dose / scal
    e0 + eMax / maxDens * (dose^delta1) * (1 - dose)^delta2
  }

  trueBiomarker <- function(dose) {
    betaMod(dose, e0 = 0.2, eMax = 0.6, delta1 = 5, delta2 = 5 * 0.5 / 0.5, scal = 100)
  }

  trueTox <- function(dose) {
    pnorm((dose - 60) / 10)
  }

  result <- simulate(
    design,
    trueTox = trueTox,
    trueBiomarker = trueBiomarker,
    sigma2W = 0.01,
    rho = 0,
    nsim = 1,
    parallel = FALSE,
    seed = 3,
    startingDose = 6,
    firstSeparate = TRUE,
    mcmcOptions = McmcOptions(
      burnin = 100,
      step = 1,
      samples = 300,
      rng_kind = "Mersenne-Twister",
      rng_seed = 1234
    )
  )

  expect_equal(result@rho_est, 0.07991541, tolerance = 1e-7) # printed result

  expect_equal(result@rho_est, 0.079915412) # actual result

  expect_equal(result@sigma2w_est, 0.03177778, tolerance = 1e-7) # printed result

  expect_equal(result@sigma2w_est, 0.031777778) # actual result

  expect_equal(any(sapply(result@fit_biomarker[[1]], is.numeric)), TRUE) # all elements of fit are numeric

  expect_equal(dim(result@fit_biomarker[[1]])[1], 11)

  expect_equal(dim(result@fit_biomarker[[1]])[2], 3)

  expect_equal(length(result@stop_report), 3) # check for length

  expect_logical(result@stop_report) # check for stop_report to be logical vector

  expect_list(result@data)

  expect_class(result@data[[1]], "Data") # check for data object has correct class

  expect_list(result@additional_stats)

  expect_list(result@additional_stats[[1]])

  expect_length(result@additional_stats[[1]], 0)

  expect_equal(result@doses, 1)
})

test_that("simulate-DualDesign produces consistent results", {
  design <- h_get_design_dualdata(TRUE)

  # define scenarios for the TRUE toxicity and efficacy profiles
  betaMod <- function(dose, e0, eMax, delta1, delta2, scal) {
    maxDens <- (delta1^delta1) * (delta2^delta2) / ((delta1 + delta2)^(delta1 + delta2))
    dose <- dose / scal
    e0 + eMax / maxDens * (dose^delta1) * (1 - dose)^delta2
  }

  trueBiomarker <- function(dose) {
    betaMod(dose, e0 = 0.2, eMax = 0.6, delta1 = 5, delta2 = 5 * 0.5 / 0.5, scal = 100)
  }

  trueTox <- function(dose) {
    pnorm((dose - 60) / 10)
  }

  result <- simulate(
    design,
    trueTox = trueTox,
    trueBiomarker = trueBiomarker,
    sigma2W = 0.01,
    rho = 0,
    nsim = 1,
    parallel = FALSE,
    seed = 3,
    startingDose = 6,
    mcmcOptions = McmcOptions(
      burnin = 100,
      step = 1,
      samples = 300,
      rng_kind = "Mersenne-Twister",
      rng_seed = 1234
    )
  )

  expect_snapshot(result)
})


test_that("simulate-TDSamplesDesign produces consistent results", {
  data <- Data(doseGrid = seq(25, 300, 25))

  model <- LogisticIndepBeta(
    binDLE = c(1.05, 1.8),
    DLEweights = c(3, 3),
    DLEdose = c(25, 300),
    data = data
  )
  tdNextBest <- NextBestTDsamples(
    prob_target_drt = 0.35,
    prob_target_eot = 0.3,
    derive = function(samples) {
      as.numeric(quantile(samples, probs = 0.3))
    }
  )
  mySize <- CohortSizeConst(size = 3)
  myIncrements <- IncrementsRelative(
    intervals = c(min(data@doseGrid), max(data@doseGrid)),
    increments = c(2, 2)
  )
  myStopping <- StoppingMinPatients(nPatients = 36)

  design <- TDsamplesDesign(
    model = model,
    nextBest = tdNextBest,
    stopping = myStopping,
    increments = myIncrements,
    cohort_size = mySize,
    data = data, startingDose = 25
  )
  myTruth <- probFunction(model, phi1 = -53.66584, phi2 = 10.50499)
  options <- McmcOptions(burnin = 100, step = 2, samples = 200)
  result <- simulate(
    object = design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    mcmcOptions = options,
    parallel = FALSE
  )

  expect_snapshot(result)
})

test_that("simulate-TDDesign produces consistent results", {
  suppressWarnings({
    design <- h_get_design_tddesign()
  })

  myTruth <- probFunction(design@model, phi1 = -53.66584, phi2 = 10.50499)

  result <- simulate(
    object = design,
    args = NULL,
    truth = myTruth,
    nsim = 1,
    seed = 819,
    parallel = FALSE
  )
  expect_snapshot(result)
})

test_that("simulate-DualResponsesDesign produces consistent results", {
  design <- h_get_design_dualresponses()
  myTruthDLE <- probFunction(design@model, phi1 = -53.66584, phi2 = 10.50499)
  myTruthEff <- efficacyFunction(design@eff_model, theta1 = -4.818429, theta2 = 3.653058)

  myTruthGain <- function(dose) {
    return((myTruthEff(dose)) / (1 + (myTruthDLE(dose) / (1 - myTruthDLE(dose)))))
  }
  options <- McmcOptions(burnin = 100, step = 2, samples = 200)
  result <- simulate(
    object = design,
    args = NULL,
    trueDLE = myTruthDLE,
    trueEff = myTruthEff,
    trueNu = 1 / 0.025,
    nsim = 1,
    seed = 819,
    parallel = FALSE
  )

  expect_snapshot(result)
})

test_that("simulate-DualResponsesSamplesDesign produces consistent results", {
  data <- DataDual(doseGrid = seq(25, 300, 25), placebo = FALSE)
  DLEmodel <- LogisticIndepBeta(
    binDLE = c(1.05, 1.8),
    DLEweights = c(3, 3),
    DLEdose = c(25, 300),
    data = data
  )

  Effmodel <- Effloglog(
    eff = c(1.223, 2.513), eff_dose = c(25, 300),
    nu = c(a = 1, b = 0.025), data = data
  )

  mynextbest <- NextBestMaxGainSamples(
    prob_target_drt = 0.35,
    prob_target_eot = 0.3,
    derive = function(samples) {
      as.numeric(quantile(samples, prob = 0.3))
    },
    mg_derive = function(mg_samples) {
      as.numeric(quantile(mg_samples, prob = 0.5))
    }
  )

  myIncrements <- IncrementsRelative(
    intervals = c(25, 300),
    increments = c(2, 2)
  )
  mySize <- CohortSizeConst(size = 3)
  myStopping <- StoppingMinPatients(nPatients = 10)

  # Specified the design
  design <- DualResponsesSamplesDesign(
    nextBest = mynextbest,
    cohort_size = mySize,
    startingDose = 25,
    model = DLEmodel,
    eff_model = Effmodel,
    data = data,
    stopping = myStopping,
    increments = myIncrements
  )

  myTruthDLE <- probFunction(design@model, phi1 = -53.66584, phi2 = 10.50499)
  myTruthEff <- efficacyFunction(design@eff_model, theta1 = -4.818429, theta2 = 3.653058)

  myTruthGain <- function(dose) {
    return((myTruthEff(dose)) / (1 + (myTruthDLE(dose) / (1 - myTruthDLE(dose)))))
  }

  options <- McmcOptions(burnin = 10, step = 1, samples = 50)
  result <- simulate(design,
    args = NULL,
    trueDLE = myTruthDLE,
    trueEff = myTruthEff,
    trueNu = 1 / 0.025,
    nsim = 1,
    mcmcOptions = options,
    seed = 819,
    parallel = FALSE
  )

  expect_snapshot(result)
})

test_that("Test if simulate generate the expected output.", {
  data <- h_get_data(placebo = FALSE)
  model <- h_get_logistic_normal()
  increments <- h_increments_relative()
  next_best <- h_next_best_ncrm()
  size <- CohortSizeConst(size = 3)

  # Extreme truth function, which has constant probability 1 in dose grid range.
  truth <- probFunction(model, alpha0 = 175, alpha1 = 5)
  stop_rule <- StoppingMinPatients(nPatients = 5)
  design <- Design(
    model = model,
    stopping = stop_rule,
    increments = increments,
    nextBest = next_best,
    cohort_size = size,
    data = data,
    startingDose = 25
  )

  my_options <- McmcOptions(burnin = 100, step = 2, samples = 5, rng_kind = "Mersenne-Twister", rng_seed = 3)

  sim <- simulate(
    design,
    nsim = 1,
    truth = truth,
    seed = 819,
    mcmcOptions = my_options
  )

  # expect_snapshot_value doesn't work here, either in toto or slot-by-slot,
  # regardless of style
  expect_snapshot(sim)
})


## NextBestInfTheory ----

test_that("NextBestInfTheory produces consistent results for empty data", {
  emptydata <- Data(doseGrid = seq(from = 40, to = 200, by = 10))

  # Set up the model; sigma0 = 1.0278, sigma1 = 1.65, rho = 0.5.
  model <- LogisticLogNormal(
    mean = c(-4.47, 0.0033),
    cov = matrix(c(1.056373, 0.847935, 0.847935, 2.722500), nrow = 2)
  )

  stop_rule <- StoppingMinPatients(nPatients = 30)
  increments <- IncrementsRelative(interval = 0, increments = 1)
  new_my_next_best <- NextBestInfTheory(target = 0.25, asymmetry = 0.1)
  cohort <- CohortSizeConst(size = 3)
  my_truth <- probFunction(model, alpha0 = 175, alpha1 = 5)

  design <- Design(
    model = model,
    stopping = stop_rule,
    increments = increments,
    nextBest = new_my_next_best,
    cohort_size = cohort,
    data = emptydata,
    startingDose = 40
  )

  sim <- simulate(
    design,
    nsim = 5,
    truth = my_truth,
    mcmcOptions = h_get_mcmc_options()
  )

  result <- summary(sim, truth = my_truth, target = new_my_next_best@target)

  expect_equal(
    result@fit_at_dose_most_selected,
    c(0.985602, 0.985602, 0.985602, 0.985602, 0.985602),
    tolerance = 1e-07
  )
  expect_equal(result@prop_dlts, rep(1L, 5))
  expect_equal(result@mean_tox_risk, rep(1L, 5))
  expect_equal(result@dose_selected, rep(40, 5))
  expect_equal(result@tox_at_doses_selected, rep(1L, 5))
  # expect_snapshot_value doesn't work here regardless of style
  expect_snapshot(result@mean_fit)
})

test_that("NextBestInfTheory produces consistent results with a dataset", {
  my_data <- h_get_data(placebo = FALSE)

  # Set up the model; sigma0 = 1.0278, sigma1 = 1.65, rho = 0.5.
  model <- LogisticLogNormal(
    mean = c(-4.47, 0.0033),
    cov = matrix(c(1.056373, 0.847935, 0.847935, 2.722500), nrow = 2)
  )

  stop_rule <- StoppingMinPatients(nPatients = 5)
  increments <- IncrementsRelative(interval = 0, increments = 1)
  new_my_next_best <- NextBestInfTheory(target = 0.25, asymmetry = 0.1)
  cohort <- CohortSizeConst(size = 3)
  my_truth <- probFunction(model, alpha0 = 175, alpha1 = 5)

  design <- Design(
    model = model,
    stopping = stop_rule,
    increments = increments,
    nextBest = new_my_next_best,
    cohort_size = cohort,
    data = my_data,
    startingDose = 25
  )

  sim <- simulate(
    design,
    nsim = 5,
    truth = my_truth,
    mcmcOptions = h_get_mcmc_options()
  )

  result <- summary(sim, truth = my_truth, target = new_my_next_best@target)
  expect_equal(
    result@fit_at_dose_most_selected,
    c(0.222, 0.222, 0.222, 0.222, 0.222),
    tolerance = 1e-02
  )
  expect_equal(result@prop_dlts, rep(0.267, 5), tolerance = 1e-02)
  expect_equal(result@mean_tox_risk, rep(1L, 5))
  expect_equal(result@dose_selected, rep(50, 5))
  expect_equal(result@tox_at_doses_selected, rep(1L, 5))
  # expect_snapshot_value doesn't work here, regardless of style
  expect_snapshot(result@mean_fit)
})

## stop_reasons integration test ----

test_that("stop_reasons can be NA with certain stopping rule settings", {
  data <- h_get_data(placebo = FALSE)
  model <- h_get_logistic_normal()
  increments <- h_increments_relative()
  next_best <- h_next_best_ncrm()
  size <- CohortSizeConst(size = 3)
  # Extreme truth function, which has constant probability 1 in dose grid range.
  truth <- probFunction(model, alpha0 = 175, alpha1 = 5)
  stopping <- StoppingMissingDose()
  design <- Design(
    model = model,
    stopping = stopping,
    increments = increments,
    nextBest = next_best,
    cohort_size = size,
    data = data,
    startingDose = 25
  )
  sim <- simulate(
    design,
    nsim = 5,
    truth = truth,
    seed = 819,
    mcmcOptions = h_get_mcmc_options()
  )
  result <- sim@stop_reasons
  # In this case the trial always stops because no dose is deemed safe enough
  # to continue the trial. This is the default behavior of the
  # stopTrial() method.
  expected <- list(
    "Next dose is NA , i.e., no active dose is safe enough according to the NextBest rule.",
    "Next dose is NA , i.e., no active dose is safe enough according to the NextBest rule.",
    "Next dose is NA , i.e., no active dose is safe enough according to the NextBest rule.",
    "Next dose is NA , i.e., no active dose is safe enough according to the NextBest rule.",
    "Next dose is NA , i.e., no active dose is safe enough according to the NextBest rule."
  )
  expect_identical(result, expected)
})

## DesignGrouped ----

test_that("simulate for DesignGrouped works as expected", {
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 9),
      increments = IncrementsDoseLevels(levels = 5),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 1
    ),
    same_dose_for_all = TRUE,
    first_cohort_mono_only = TRUE
  )
  my_truth <- function(x) plogis(-4 + 0.2 * log(x / 0.1))
  my_combo_truth <- function(x) plogis(-4 + 0.5 * log(x / 0.1))

  result <- expect_silent(simulate(
    object,
    nsim = 2,
    seed = 123,
    truth = my_truth,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options()
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  mono_trial <- result$mono@data[[2L]]
  combo_trial <- result$combo@data[[2L]]

  # First cohort is only mono at the starting dose (lowest in dose grid).
  expect_true(all(mono_trial@xLevel[1:3] == 1))

  # We have the same dose for subsequent cohorts.
  expect_true(all(mono_trial@xLevel[4:6] == combo_trial@xLevel[1:3]))
  expect_true(all(mono_trial@xLevel[7:9] == combo_trial@xLevel[4:6]))
})

test_that("simulate for DesignGrouped works as expected with different doses, parallel first cohort", {
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 20),
      increments = IncrementsDoseLevels(levels = 5),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 1
    ),
    same_dose_for_all = FALSE,
    first_cohort_mono_only = FALSE
  )
  my_truth <- function(x) plogis(-4 + 0.2 * log(x / 0.1))
  my_combo_truth <- function(x) plogis(-4 + 0.9 * log(x / 0.1))

  result <- expect_silent(simulate(
    object,
    nsim = 1,
    seed = 123,
    truth = my_truth,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options()
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  mono_trial <- result$mono@data[[1L]]
  combo_trial <- result$combo@data[[1L]]

  # First cohort is joint at starting dose.
  expect_true(all(mono_trial@xLevel[1:3] == combo_trial@xLevel[1:3]))

  # We have different doses in subsequent cohorts.
  expect_false(all(mono_trial@xLevel[4:20] == combo_trial@xLevel[4:20]))
})

test_that("simulate for DesignGrouped works when first patient is dosed separately, different combo design", {
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 10),
      increments = IncrementsDoseLevels(levels = 3),
      cohort_size = CohortSizeConst(2),
      data = Data(doseGrid = 1:100),
      startingDose = 10
    ),
    combo = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 20),
      increments = IncrementsDoseLevels(levels = 5),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 1
    ),
    same_dose_for_all = FALSE,
    first_cohort_mono_only = FALSE
  )
  my_truth <- function(x) plogis(-4 + 0.2 * log(x / 0.1))
  my_combo_truth <- function(x) plogis(-2 + 0.9 * log(x / 0.1))

  result <- expect_silent(simulate(
    object,
    nsim = 1,
    seed = 123,
    truth = my_truth,
    firstSeparate = TRUE,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options(),
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  mono_trial <- result$mono@data[[1L]]
  combo_trial <- result$combo@data[[1L]]

  # We expect at least one cohort with just one patient in the combo arm here
  # because of the high toxicity.
  expect_true(any(table(combo_trial@cohort) == 1))

  # Check that we had the different cohort sizes between the two arms.
  expect_true(max(table(mono_trial@cohort)) == 2)
  expect_true(max(table(combo_trial@cohort)) == 3)

  # Check that we had different starting doses in the two arms.
  expect_true(mono_trial@x[1] == 10)
  expect_true(combo_trial@x[1] == 1)
})

test_that("simulate for DesignGrouped works with different starting doses and first mono", {
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 9),
      increments = IncrementsDoseLevels(levels = 5),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 10
    ),
    combo = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 9),
      increments = IncrementsRelative(c(0, 100), c(2, 1)),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 5
    ),
    same_dose_for_all = FALSE,
    first_cohort_mono_only = TRUE
  )
  my_truth <- function(x) plogis(-4 + 0.2 * log(x / 0.1))
  my_combo_truth <- function(x) plogis(-4 + 0.5 * log(x / 0.1))

  result <- expect_silent(simulate(
    object,
    nsim = 2,
    seed = 123,
    truth = my_truth,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options()
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  mono_trial <- result$mono@data[[2L]]
  combo_trial <- result$combo@data[[2L]]

  # First cohort is only mono at the starting dose.
  expect_true(all(mono_trial@xLevel[1:3] == 10))

  # In first combo cohort we have the expected starting dose.
  expect_true(all(combo_trial@xLevel[1:3] == 5))
})

test_that("simulate for DesignGrouped allows to stop mono when combo stops", {
  mono_arm <- Design(
    model = .LogisticNormal(),
    nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
    # With a custom label that we can check below.
    stopping = StoppingMinPatients(nPatients = 20, report_label = "my label"),
    increments = IncrementsDoseLevels(levels = 5),
    cohort_size = CohortSizeConst(3),
    data = Data(doseGrid = 1:100),
    startingDose = 10
  )
  combo_arm <- .Design(
    mono_arm,
    # Such that we stop after the first cohort.
    stopping = StoppingMinPatients(nPatients = 1)
  )
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = mono_arm,
    combo = combo_arm,
    same_dose_for_all = FALSE,
    first_cohort_mono_only = FALSE,
    stop_mono_with_combo = TRUE
  )
  my_truth <- function(x) plogis(-4 + 0.2 * log(x / 0.1))
  my_combo_truth <- function(x) plogis(-4 + 0.5 * log(x / 0.1))

  result <- expect_silent(simulate(
    object,
    nsim = 2,
    seed = 123,
    truth = my_truth,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options()
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  # We see the expected stop reasons.
  lapply(
    result$mono@stop_reasons,
    function(x) expect_subset("Based on external result stop", unlist(x))
  )
  expect_identical(
    result$combo@stop_reasons,
    rep(list("Number of patients is 3 and thus reached the prespecified minimum number 1"), 2)
  )

  # But mono still had the initial 3 patients in both simulations.
  expect_identical(
    lapply(result$mono@data, slot, "nObs"),
    rep(list(3L), 2)
  )

  # And we see the stop report includes the previous stopping rule too.
  expect_identical(
    colnames(result$mono@stop_report),
    c(NA, "my label", "Stop Mono with Combo")
  )
})

test_that("simulate for DesignGrouped reports correctly when mono is not stopped because of combo", {
  mono_arm <- Design(
    model = .LogisticNormal(),
    nextBest = NextBestNCRM(target = c(0.2, 0.4), overdose = c(0.4, 1), max_overdose_prob = 0.7),
    # With a custom label that we can check below.
    stopping = StoppingTargetProb(report_label = "my label"),
    increments = IncrementsDoseLevels(levels = 5),
    cohort_size = CohortSizeConst(3),
    data = Data(doseGrid = 1:100),
    startingDose = 10
  )
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = mono_arm,
    combo = mono_arm,
    same_dose_for_all = FALSE,
    first_cohort_mono_only = FALSE,
    stop_mono_with_combo = TRUE
  )
  my_truth <- function(x) plogis(-4 + 0.2 * log(x / 0.1))
  my_combo_truth <- function(x) plogis(-4 + 0.5 * log(x / 0.1))

  set.seed(123)
  result <- expect_silent(simulate(
    object,
    nsim = 2,
    seed = 123,
    truth = my_truth,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options()
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  # We see the stop report includes the previous stopping rule and the mono because combo thing too.
  expect_identical(
    colnames(result$mono@stop_report),
    c(NA, "my label", "Stop Mono with Combo")
  )
  # But not for the combo.
  expect_identical(
    colnames(result$combo@stop_report),
    "my label"
  )
})

test_that("simulate for DesignGrouped works with parallel start when first cohort mono only", {
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 9),
      increments = IncrementsDoseLevels(levels = 5),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 1
    ),
    combo = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 9),
      increments = IncrementsRelative(c(0, 100), c(2, 1)),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 1
    ),
    same_dose_for_all = FALSE,
    first_cohort_mono_only = TRUE,
    same_dose_for_start = TRUE
  )
  my_truth <- function(x) plogis(-4 + 0.2 * log(x / 0.1))
  my_combo_truth <- function(x) plogis(-4 + 0.5 * log(x / 0.1))

  result <- expect_silent(simulate(
    object,
    nsim = 2,
    seed = 123,
    truth = my_truth,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options()
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  mono_trial <- result$mono@data[[1L]]
  combo_trial <- result$combo@data[[1L]]

  # First cohort is only mono at the starting dose.
  expect_true(all(mono_trial@x[1:3] == 1))

  # Second cohort in mono is again the starting dose because of parallel start.
  expect_true(all(mono_trial@x[4:6] == 1))

  # In first combo cohort we have the expected starting dose.
  expect_true(all(combo_trial@x[1:3] == 1))
})

test_that("simulate for DesignGrouped works with parallel start when first cohort mono and combo", {
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 9),
      increments = IncrementsDoseLevels(levels = 5),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 1
    ),
    combo = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 9),
      increments = IncrementsRelative(c(0, 100), c(2, 1)),
      cohort_size = CohortSizeConst(3),
      data = Data(doseGrid = 1:100),
      startingDose = 3
    ),
    same_dose_for_all = FALSE,
    first_cohort_mono_only = FALSE,
    same_dose_for_start = TRUE
  )
  my_truth <- function(x) plogis(-4 + 0.2 * log(x / 0.1))
  my_combo_truth <- function(x) plogis(-4 + 0.5 * log(x / 0.1))

  result <- expect_silent(simulate(
    object,
    nsim = 2,
    seed = 123,
    truth = my_truth,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options()
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  mono_trial <- result$mono@data[[1L]]
  combo_trial <- result$combo@data[[1L]]

  # First cohort is only mono at the starting dose.
  expect_true(all(mono_trial@x[1:3] == 1))

  # In first combo cohort we have the lower, mono starting dose too, because
  # of parallel start.
  expect_true(all(combo_trial@x[1:3] == 1))
})

test_that("simulate for DesignGrouped uses DLT probabilities and cohort sizes correctly", {
  object <- DesignGrouped(
    model = LogisticLogNormalGrouped(mean = rep(-1, 4), cov = diag(5, 4), ref_dose = 1),
    mono = Design(
      model = .LogisticNormal(),
      nextBest = NextBestNCRM(target = c(0.3, 0.6), overdose = c(0.6, 1), max_overdose_prob = 0.7),
      stopping = StoppingMinPatients(nPatients = 9),
      increments = IncrementsDoseLevels(levels = 5),
      # Make sure that cohort size changes across doses.
      cohort_size = CohortSizeRange(intervals = c(0, 5, 10), cohort_size = c(1, 2, 3)),
      data = Data(doseGrid = 1:100),
      startingDose = 1
    ),
    same_dose_for_all = TRUE,
    first_cohort_mono_only = TRUE
  )
  # Make sure that the true DLT probabilities change heavily across doses.
  my_truth <- stats::stepfun(x = c(0, 5, 10), y = c(0, 0, 0.5, 0.99))
  my_combo_truth <- my_truth

  result <- expect_silent(simulate(
    object,
    nsim = 2,
    seed = 123,
    truth = my_truth,
    combo_truth = my_combo_truth,
    mcmcOptions = h_get_mcmc_options()
  ))

  expect_list(result)
  expect_names(names(result), identical.to = c("mono", "combo"))
  expect_valid(result$mono, "Simulations")
  expect_valid(result$combo, "Simulations")

  mono_trial <- result$mono@data[[2L]]
  combo_trial <- result$combo@data[[2L]]

  # We have seen doses of 5 or larger.
  expect_true(any(mono_trial@x >= 5))
  expect_true(any(combo_trial@x >= 5))

  # And therefore we must have seen DLTs.
  expect_true(any(mono_trial@y > 0))
  expect_true(any(combo_trial@y > 0))

  # And therefore we must also have had cohort sizes larger than 1.
  expect_true(any(table(mono_trial@cohort) > 1))
  expect_true(any(table(combo_trial@cohort) > 1))
})

# examine ----

## DADesign ----

test_that("examine for DADesign works as expected", {
  t_max <- 60

  emptydata <- DataDA(
    doseGrid = c(
      0.1, 0.5, 1, 1.5, 3, 6,
      seq(from = 10, to = 80, by = 2)
    ),
    Tmax = t_max
  )

  n_pieces <- 10
  lambda_prior <- function(k) {
    n_pieces / (t_max * (n_pieces - k + 0.5))
  }
  l_vector <- as.numeric(t(apply(
    as.matrix(c(1:n_pieces), 1, n_pieces),
    2,
    lambda_prior
  )))

  model <- DALogisticLogNormal(
    mean = c(-0.85, 1),
    cov = matrix(c(1, -0.5, -0.5, 1), nrow = 2),
    ref_dose = 56,
    npiece = n_pieces,
    l = l_vector
  )

  my_increments <- IncrementsRelative(
    intervals = c(0, 20),
    increments = c(1, 0.33)
  )
  my_next_best <- NextBestNCRM(
    target = c(0.2, 0.35),
    overdose = c(0.35, 1),
    max_overdose_prob = 0.25
  )
  my_size_1 <- CohortSizeRange(
    intervals = c(0, 30),
    cohort_size = c(1, 3)
  )
  my_size_2 <- CohortSizeDLT(
    intervals = c(0, 1),
    cohort_size = c(1, 3)
  )
  my_size <- maxSize(my_size_1, my_size_2)
  my_stopping_1 <- StoppingTargetProb(
    target = c(0.2, 0.35),
    prob = 0.5
  )
  my_stopping_2 <- StoppingMinPatients(nPatients = 50)
  my_stopping <- (my_stopping_1 | my_stopping_2)
  my_safety <- SafetyWindowConst(c(6, 2), 7, 7)

  design <- DADesign(
    model = model,
    increments = my_increments,
    nextBest = my_next_best,
    stopping = my_stopping,
    cohort_size = my_size,
    data = emptydata,
    safetyWindow = my_safety,
    startingDose = 3
  )

  options <- McmcOptions(
    burnin = 10,
    step = 1,
    samples = 100,
    rng_kind = "Mersenne-Twister",
    rng_seed = 12
  )
  expect_warning(
    result <- examine(design, mcmcOptions = options, maxNoIncrement = 2L),
    "Stopping because 2 times no increment"
  )
  expect_data_frame(result)
  expect_named(result, c("DLTsearly_1", "dose", "DLTs", "nextDose", "stop", "increment"))
})

## Design ----

test_that("examine produces consistent results", {
  design <- h_get_design_data()
  options <- h_get_mcmc_options()

  result <- examine(design, mcmcOptions = options)

  expect_snapshot(result)
})

test_that("examine produces consistent results with placebo data", {
  design <- h_get_design_data(TRUE)
  options <- h_get_mcmc_options()

  result <- examine(design, mcmcOptions = options)

  expect_snapshot(result)
})

## RuleDesign ----

test_that("simulate-RuleDesign produces consistent results", {
  design <- ThreePlusThreeDesign(doseGrid = c(1, 3, 5, 10, 15, 20, 25, 40, 50, 80, 100))
  options <- h_get_mcmc_options()

  result <- examine(design, mcmcOptions = options)

  expect_snapshot(result)
})

# tidy ----

## DualDesign ----
test_that("tidy-DualDesign works correctly", {
  obj <- .DefaultDualDesign()
  result <- tidy(obj)
  # style = "deparse" fails with Could not find function numeric
  expect_snapshot_value(result, style = "serialize")
})
