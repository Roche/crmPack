# This code mocks functions that have long execution times so that unit tests
# complete more quickly.  Initial tests suggest that the mocks need to be defined
# in the file in which the tests are executed.  `source`ing the mocks does not
# work.
#
# The persistent objects that are loaded are created by
# /testthat/fixtures/make_persistent_objects_for_mocked_constructors.R.

testthat::local_mocked_bindings(
  .DefaultDASimulations = function(...) {
    readRDS(testthat::test_path("fixtures", "default_da_simulations.Rds"))
  }
)

testthat::local_mocked_bindings(
  .DefaultSimulations = function(...) {
    readRDS(testthat::test_path("fixtures", "default_simulations.Rds"))
  }
)

testthat::local_mocked_bindings(
  .DefaultDualSimulationsSummary = function(...) {
    readRDS(testthat::test_path("fixtures", "default_dual_simulations_summary.Rds"))
  }
)
# End of mocks

# GeneralSimulations-class ----
test_that("GeneralSimulations generator function works as expected", {
  result <- expect_silent(.GeneralSimulations())
  expect_valid(result, "GeneralSimulations")
})

test_that("GeneralSimulations object can be created with the user constructor", {
  data <- list(
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
  )

  doses <- c(1, 2)

  seed <- as.integer(123)

  result <- expect_silent(
    GeneralSimulations(
      data,
      doses,
      seed
    )
  )

  expect_valid(result, "GeneralSimulations")
  expect_identical(result@data, data)
  expect_identical(result@doses, doses)
  expect_identical(result@seed, seed)
})

test_that("GeneralSimulations user constructor arguments names are as expected", {
  expect_function(
    GeneralSimulations,
    args = c("data", "doses", "seed"),
    ordered = TRUE
  )
})

# Simulations-class ----
test_that("Simulations generator function works as expected", {
  result <- expect_silent(.Simulations())
  expect_valid(result, "Simulations")
})

test_that("Simulations object can be created with the user constructor", {
  fit <- list(
    c(0.1, 0.2),
    c(0.3, 0.4)
  )
  stop_reasons <- list("A", "B")

  stop_report <- matrix(c(TRUE, FALSE), nrow = 2)

  additional_stats <- list(a = 1, b = 1)

  data <- list(
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
  )

  doses <- c(1, 2)

  seed <- as.integer(123)

  result <- expect_silent(
    Simulations(
      fit = fit,
      stop_reasons = stop_reasons,
      stop_report = stop_report,
      additional_stats = additional_stats,
      data,
      doses,
      seed
    )
  )

  expect_valid(result, "Simulations")
  expect_identical(result@fit, fit)
  expect_identical(result@stop_reasons, stop_reasons)
})

test_that("Simulations user constructor arguments names are as expected", {
  expect_function(
    Simulations,
    args = c("fit", "stop_reasons", "stop_report", "additional_stats", "..."),
    ordered = TRUE
  )
})

# DualSimulations-class ----
test_that("DualSimulations generator function works as expected", {
  result <- expect_silent(.DualSimulations())
  expect_valid(result, "DualSimulations")
})

test_that("DualSimulations object can be created with the user constructor", {
  rho_est <- c(0.25, 0.35)
  sigma2w_est <- c(0.15, 0.25)
  fit_biomarker <- list(c(0.3, 0.4), c(0.4, 0.5))

  data_list <- list(
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
  )

  doses <- c(1, 2)
  seed <- as.integer(123)

  fit <- list(
    c(0.1, 0.2),
    c(0.3, 0.4)
  )

  stop_report <- matrix(c(TRUE, FALSE), nrow = 2)

  stop_reasons <- list("A", "B")

  additional_stats <- list(a = 1, b = 1)

  result <- expect_silent(
    DualSimulations(
      rho_est = rho_est,
      sigma2w_est = sigma2w_est,
      fit_biomarker = fit_biomarker,
      fit = fit,
      stop_report = stop_report,
      stop_reasons = stop_reasons,
      additional_stats = additional_stats,
      data = data_list,
      doses = doses,
      seed = seed
    )
  )

  expect_valid(result, "DualSimulations")
  expect_identical(result@rho_est, rho_est)
  expect_identical(result@sigma2w_est, sigma2w_est)
  expect_identical(result@fit_biomarker, fit_biomarker)
})

test_that("DualSimulations user constructor arguments names are as expected", {
  expect_function(
    DualSimulations,
    args = c("rho_est", "sigma2w_est", "fit_biomarker", "..."),
    ordered = TRUE
  )
})

# GeneralSimulationsSummary ----
test_that("GeneralSimulationsSummary generates object correctly", {
  target_value <- 1
  target_dose_interval_value <- 2
  nsim_value <- 3L
  mean_tox_risk_value <- 4
  dose_selected_value <- 5

  result <- expect_silent(
    .GeneralSimulationsSummary(
      target = target_value,
      target_dose_interval = target_dose_interval_value,
      nsim = nsim_value,
      prop_dlts = list(),
      mean_tox_risk = mean_tox_risk_value,
      dose_selected = dose_selected_value,
      tox_at_doses_selected = 6,
      prop_at_target = 7,
      dose_most_selected = 8,
      obs_tox_rate_at_dose_most_selected = 9,
      n_obs = list(),
      n_above_target = 10L,
      dose_grid = 11,
      placebo = TRUE
    )
  )

  expect_valid(result, "GeneralSimulationsSummary")

  expect_identical(result@target, target_value)
  expect_identical(result@target_dose_interval, target_dose_interval_value)
  expect_identical(result@nsim, nsim_value)
  expect_identical(result@mean_tox_risk, mean_tox_risk_value)
  expect_identical(result@dose_selected, dose_selected_value)
})

test_that("GeneralSimulationsSummary cannot be instantiated directly", {
  expect_error(.DefaultGeneralSimulationsSummary(),
    "Class GeneralSimulationsSummary cannot be instantiated directly",
    fixed = FALSE
  )
})

# DualSimulationsSummary ----
test_that("DualSimulationsSummary object can be created with the user constructor", {
  biomarker_fit_at_dose_most_selected <- 0.3
  mean_biomarker_fit <- list(c(0.25, 0.5, 0.75)) # This should be a list

  result <- expect_silent(
    .DualSimulationsSummary(
      biomarker_fit_at_dose_most_selected = biomarker_fit_at_dose_most_selected,
      mean_biomarker_fit = mean_biomarker_fit
    )
  )

  expect_valid(result, "DualSimulationsSummary")
  expect_identical(result@biomarker_fit_at_dose_most_selected, biomarker_fit_at_dose_most_selected)
  expect_identical(result@mean_biomarker_fit, mean_biomarker_fit)
})

test_that("DualSimulationsSummary generator function works as expected", {
  result <- expect_silent(.DefaultDualSimulationsSummary())
  expect_valid(result, "DualSimulations")
})

# PseudoSimulations-class ----
test_that("PseudoSimulations generator function works as expected", {
  result <- expect_silent(.PseudoSimulations())
  expect_valid(result, "PseudoSimulations")
})

test_that("PseudoSimulations object can be created with the user constructor", {
  fit <- list(c(0.1, 0.2), c(0.3, 0.4))

  final_td_target_during_trial_estimates <- c(0.1, 0.2)
  final_td_target_end_of_trial_estimates <- c(0.1, 0.2)

  final_td_target_during_trial_at_dose_grid <- c(0.1, 0.2)
  final_td_target_end_of_trial_at_dose_grid <- c(0.1, 0.2)

  final_tdeot_cis <- list(c(0.1, 0.2), c(0.1, 0.2))
  final_tdeot_ratios <- c(0.1, 0.2)

  final_cis <- list(c(0.1, 0.2), c(0.1, 0.2))
  final_ratios <- c(0.1, 0.2)

  stop_report <- matrix(TRUE, nrow = 2)
  stop_reasons <- list("A", "B")

  doses <- c(100, 200)

  data <- list(
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
  )

  result <- expect_silent(
    PseudoSimulations(
      fit = fit,
      data = data,
      doses = doses,
      final_td_target_during_trial_estimates = final_td_target_during_trial_estimates,
      final_td_target_end_of_trial_estimates = final_td_target_end_of_trial_estimates,
      final_td_target_during_trial_at_dose_grid = final_td_target_during_trial_at_dose_grid,
      final_td_target_end_of_trial_at_dose_grid = final_td_target_end_of_trial_at_dose_grid,
      final_tdeot_cis = final_tdeot_cis,
      final_tdeot_ratios = final_tdeot_ratios,
      final_cis = final_cis,
      final_ratios = final_ratios,
      stop_report = stop_report,
      stop_reasons = stop_reasons,
      seed = 123
    )
  )

  expect_valid(result, "PseudoSimulations")
  expect_identical(result@fit, fit)
  expect_identical(result@stop_reasons, stop_reasons)
})

test_that("PseudoSimulations user constructor argument names are as expected", {
  expect_function(
    PseudoSimulations,
    args = c(
      "fit",
      "final_td_target_during_trial_estimates",
      "final_td_target_end_of_trial_estimates",
      "final_td_target_during_trial_at_dose_grid",
      "final_td_target_end_of_trial_at_dose_grid",
      "final_tdeot_cis",
      "final_tdeot_ratios",
      "final_cis",
      "final_ratios",
      "stop_report",
      "stop_reasons",
      "..."
    ),
    ordered = TRUE
  )
})

test_that(".DefaultPseudoSimulations cannot be instantiated directly", {
  expect_error(.DefaultPseudoSimulations(),
    "Class PseudoSimulations cannot be instantiated directly. Please use one of its subclasses instead.",
    fixed = FALSE
  )
})

# PseudoDualSimulations-class ----
test_that("PseudoDualSimulations generator does not throw error and validates", {
  result <- expect_silent(.PseudoDualSimulations())
  expect_valid(result, "PseudoDualSimulations")
})

test_that("PseudoDualSimulations object can be created with the user constructor", {
  fit <- list(c(0.1, 0.2), c(0.3, 0.4))

  fit_eff <- list(
    c(0.1, 0.2),
    c(0.3, 0.4)
  )

  final_gstar_estimates <- c(0.05, 0.06)
  final_gstar_at_dose_grid <- c(0.07, 0.08)
  final_gstar_cis <- list(
    c(0.1, 0.2),
    c(0.2, 0.3)
  )
  final_gstar_ratios <- c(0.2, 0.2)
  final_optimal_dose <- c(1, 2)
  final_optimal_dose_at_dose_grid <- c(3, 4)
  sigma2_est <- c(0.001, 0.002)

  final_td_target_during_trial_estimates <- c(0.1, 0.2)
  final_td_target_end_of_trial_estimates <- c(0.1, 0.2)

  final_td_target_during_trial_at_dose_grid <- c(0.1, 0.2)
  final_td_target_end_of_trial_at_dose_grid <- c(0.1, 0.2)

  final_tdeot_cis <- list(c(0.1, 0.2), c(0.1, 0.2))
  final_tdeot_ratios <- c(0.1, 0.2)

  final_cis <- list(c(0.1, 0.2), c(0.1, 0.2))
  final_ratios <- c(0.1, 0.2)

  stop_report <- matrix(TRUE, nrow = 2)
  stop_reasons <- list("A", "B")

  data <- list(
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
  )

  doses <- c(1, 2)

  seed <- as.integer(123)

  result <- expect_silent(
    PseudoDualSimulations(
      fit = fit,
      data = data,
      doses = doses,
      fit_eff = fit_eff,
      final_gstar_estimates = final_gstar_estimates,
      final_gstar_at_dose_grid = final_gstar_at_dose_grid,
      final_gstar_cis = final_gstar_cis,
      final_gstar_ratios = final_gstar_ratios,
      final_optimal_dose = final_optimal_dose,
      final_optimal_dose_at_dose_grid = final_optimal_dose_at_dose_grid,
      final_td_target_during_trial_estimates = final_td_target_during_trial_estimates,
      final_td_target_end_of_trial_estimates = final_td_target_end_of_trial_estimates,
      final_td_target_during_trial_at_dose_grid = final_td_target_during_trial_at_dose_grid,
      final_td_target_end_of_trial_at_dose_grid = final_td_target_end_of_trial_at_dose_grid,
      final_tdeot_cis = final_tdeot_cis,
      final_tdeot_ratios = final_tdeot_ratios,
      final_cis = final_cis,
      final_ratios = final_ratios,
      stop_report = stop_report,
      stop_reasons = stop_reasons,
      sigma2_est = sigma2_est,
      seed = seed
    )
  )

  expect_valid(result, "PseudoDualSimulations")
  expect_identical(result@fit_eff, fit_eff)
  expect_identical(result@final_gstar_estimates, final_gstar_estimates)
  expect_identical(result@final_gstar_at_dose_grid, final_gstar_at_dose_grid)
  expect_identical(result@final_gstar_cis, final_gstar_cis)
  expect_identical(result@final_gstar_ratios, final_gstar_ratios)
  expect_identical(result@final_optimal_dose, final_optimal_dose)
  expect_identical(result@final_optimal_dose_at_dose_grid, final_optimal_dose_at_dose_grid)
  expect_identical(result@sigma2_est, sigma2_est)
})

test_that("PseudoDualSimulations user constructor argument names are as expected", {
  expect_function(
    PseudoDualSimulations,
    args = c(
      "fit_eff", "final_gstar_estimates", "final_gstar_at_dose_grid", "final_gstar_cis",
      "final_gstar_ratios", "final_optimal_dose", "final_optimal_dose_at_dose_grid",
      "sigma2_est", "..."
    ),
    ordered = TRUE
  )
})

# PseudoDualFlexiSimulations-class ----
test_that("PseudoDualFlexiSimulations can be generated without error and return a valid object", {
  result <- expect_silent(.PseudoDualFlexiSimulations())
  expect_valid(result, "PseudoDualFlexiSimulations")
})

test_that("PseudoDualFlexiSimulations can be instantiated using the constructor", {
  fit_eff <- list(c(0.1, 0.2), c(0.3, 0.4))
  final_gstar_estimates <- c(0.1, 0.2)
  final_gstar_at_dose_grid <- c(0.3, 0.4)
  final_gstar_cis <- list(c(0.1, 0.2), c(0.3, 0.4))
  final_gstar_ratios <- c(0.1, 0.2)
  final_optimal_dose <- c(0.5, 0.6)
  final_optimal_dose_at_dose_grid <- c(0.7, 0.8)
  sigma2_est <- c(0.01, 0.02)
  sigma2_beta_west <- c(0.03, 0.04)

  fit <- list(c(0.1, 0.2), c(0.3, 0.4))
  final_td_target_during_trial_estimates <- c(0.5, 0.6)
  final_td_target_end_of_trial_estimates <- c(0.7, 0.8)
  final_td_target_during_trial_at_dose_grid <- c(0.9, 1.0)
  final_td_target_end_of_trial_at_dose_grid <- c(1.1, 1.2)
  final_tdeot_cis <- list(c(0.1, 0.2), c(0.3, 0.4))
  final_tdeot_ratios <- c(0.5, 0.6)
  final_cis <- list(c(0.7, 0.8), c(0.9, 1.0))
  final_ratios <- c(1.1, 1.2)
  stop_report <- matrix(TRUE, nrow = 2)
  stop_reasons <- list("A", "B")

  data <- list(
    Data(
      x = 1:3,
      y = c(0, 1, 0), # Adjusted values to meet the constraint
      doseGrid = 1:3,
      ID = 1L:3L,
      cohort = 1L:3L
    ),
    Data(
      x = 4:6,
      y = c(1, 0, 1), # Adjusted values to meet the constraint
      doseGrid = 4:6,
      ID = 1L:3L,
      cohort = 1L:3L
    )
  )

  doses <- c(1, 2)
  seed <- as.integer(123)

  sim_obj <- PseudoDualFlexiSimulations(
    fit_eff = fit_eff,
    final_gstar_estimates = final_gstar_estimates,
    final_gstar_at_dose_grid = final_gstar_at_dose_grid,
    final_gstar_cis = final_gstar_cis,
    final_gstar_ratios = final_gstar_ratios,
    final_optimal_dose = final_optimal_dose,
    final_optimal_dose_at_dose_grid = final_optimal_dose_at_dose_grid,
    sigma2_est = sigma2_est,
    sigma2_beta_west = sigma2_beta_west,
    fit = fit,
    data = data,
    doses = doses,
    final_td_target_during_trial_estimates = final_td_target_during_trial_estimates,
    final_td_target_end_of_trial_estimates = final_td_target_end_of_trial_estimates,
    final_td_target_during_trial_at_dose_grid = final_td_target_during_trial_at_dose_grid,
    final_td_target_end_of_trial_at_dose_grid = final_td_target_end_of_trial_at_dose_grid,
    final_tdeot_cis = final_tdeot_cis,
    final_tdeot_ratios = final_tdeot_ratios,
    final_cis = final_cis,
    final_ratios = final_ratios,
    stop_report = stop_report,
    stop_reasons = stop_reasons,
    seed = seed
  )

  expect_valid(sim_obj, "PseudoDualFlexiSimulations")
  expect_identical(sim_obj@sigma2_beta_west, sigma2_beta_west)
})

test_that("PseudoDualFlexiSimulations user constructor argument names", {
  expect_function(
    PseudoDualFlexiSimulations,
    args = c(
      "sigma2_beta_west", "..."
    ),
    ordered = TRUE
  )
})
