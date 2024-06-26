# GeneralSimulations ----

#' Internal Helper Functions for Validation of [`GeneralSimulations`] Objects
#'
#' @description `r lifecycle::badge("stable")`
#'
#' These functions are only used internally to validate the format of an input
#' [`GeneralSimulations`] or inherited classes and therefore not exported.
#'
#' @name v_general_simulations
#' @param object (`GeneralSimulations`)\cr object to validate.
#' @return A `character` vector with the validation failure messages,
#'   or `TRUE` in case validation passes.
NULL

#' @describeIn v_general_simulations validates that the [`GeneralSimulations`] object
#'   contains valid `data` object and valid `dose` simulations.

v_general_simulations <- function(object) {
  v <- Validate()

  nSims <- length(object@data)

  v$check(
    all(sapply(object@data, is, "Data")),
    "all data elements must be Data objects"
  )
  v$check(
    identical(length(object@doses), nSims),
    "doses must have same length as the data list"
  )

  v$result()
}

#' @describeIn v_general_simulations validates that the [`Simulations`] object
#'   contains valid object `fit`, `stop_reasons`, `stop_report`, and
#'   `additional_stats` compared to the general class [`GeneralSimulations`].
#'
v_simulations <- function(object) {
  v <- Validate()

  nSims <- length(object@data)

  v$check(
    identical(length(object@fit), nSims),
    "fit must have same length as data"
  )
  v$check(
    identical(length(object@stop_reasons), nSims),
    "stop_reasons must have same length as data"
  )

  v$check(
    checkmate::test_matrix(object@stop_report,
      mode = "logical",
      nrows = nSims,
      min.cols = 1,
      any.missing = FALSE
    ),
    "stop_report must be a matrix of mode logical in which the number of rows
    equals the number of simulations and which must not contain any missing values"
  )

  v$result()
}

#' @describeIn v_general_simulations validates that the [`DualSimulations`] object and
#' capture the dose-biomarker `fits`, and the `sigma2W` and `rho` estimates.
#'
v_dual_simulations <- function(object) {
  v <- Validate()

  nSims <- length(object@data)

  v$check(
    identical(length(object@fit_biomarker), nSims),
    "fit_biomarker list has to have same length as data"
  )
  v$check(
    identical(length(object@rho_est), nSims),
    "rho_est vector has to have same length as data"
  )
  v$check(
    identical(length(object@sigma2w_est), nSims),
    "sigma2w_est has to have same length as data"
  )

  v$result()
}

# PseudoSimulations ----

#' Internal Helper Functions for Validation of [`PseudoSimulations`] Objects
#'
#' @description `r lifecycle::badge("stable")`
#'
#' These functions are only used internally to validate the format of an input
#' [`PseudoSimulations`] or inherited classes and therefore not exported.
#'
#' @name v_pseudo_simulations
#' @param object (`PseudoSimulations`)\cr object to validate.
#' @return A `character` vector with the validation failure messages,
#'   or `TRUE` in case validation passes.
NULL

#' @describeIn v_pseudo_simulations validates that the [`PseudoSimulations`] object
#'   contains valid `fit`, `FinalTDtargetEndOfTrialEstimates` ,
#'   `FinalTDtargetDuringTrialAtDoseGrid`,`FinalTDtargetEndOfTrialAtDoseGrid` ,
#'   `FinalTDEOTCIs`, `FinalTDEOTRatios`, `FinalCIs`, `FinalRatios`,
#'     object and valid `stopReasons` simulations.

v_pseudo_simulations <- function(object) {
  v <- Validate()

  nSims <- length(object@data)
  v$check(
    identical(length(object@stop_reasons), nSims),
    "stopReasons must have same length as data"
  )

  v$result()
}

#' @describeIn v_pseudo_simulations validates that the [`PseudoDualSimulations`] object
#'   contains valid `fit_eff`, `final_gstar_estimates` , `final_gstar_at_dose_grid`,
#'    `final_gstar_cis` , `final_gstar_ratios`, `final_optimal_dose`, `final_optimal_dose_at_dose_grid`
#'     object and valid `sigma2_est` simulations.

v_pseudo_dual_simulations <- function(object) {
  v <- Validate()
  nSims <- length(object@data)
  v$check(
    identical(length(object@sigma2_est), nSims),
    "sigma2_est has to have same length as data"
  )
  v$result()
}

#' @describeIn v_pseudo_simulations validates that the [`PseudoDualFlexiSimulations`]
#' object contains valid `sigma2_beta_west` vector of the final posterior mean
#' sigma2betaW estimates.`FinalGstarEstimates` , `FinalGstarAtDoseGrid`,
#'
v_pseudo_dual_flex_simulations <- function(object) {
  v <- Validate()
  nSims <- length(object@data)
  v$check(
    identical(length(object@sigma2_beta_west), nSims),
    "sigma2_beta_west has to have same length as data"
  )
  v$result()
}

#' @describeIn v_general_simulations validates that the [`DASimulations`] object
#'   contains valid `trialduration` the vector of trial duration values for all
#'   simulations.

v_da_simulations <- function(object) {
  v <- Validate()

  nSims <- length(object@data)

  v$check(
    identical(length(object@trialduration), nSims),
    "trialduration vector has to have same length as data"
  )

  v$result()
}
