#' @include helpers.R
#' @include ModelParams-validity.R
#' @include CrmPackClass-class.R
NULL

# ModelParamsNormal ----

## class ----

#' `ModelParamsNormal`
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' [`ModelParamsNormal`] is the class for a bivariate normal model parameters,
#' i.e. the mean vector, covariance matrix and precision matrix.
#' The precision matrix is an inverse of the covariance matrix in the
#' `JAGS` and it is computed internally by the object constructor function.
#'
#' @slot mean (`numeric`)\cr the mean vector.
#' @slot cov (`matrix`)\cr the covariance matrix.
#' @slot prec (`matrix`)\cr the precision matrix, which is an inverse matrix of the `cov`.
#'
#' @seealso [`ModelLogNormal`], [`LogisticNormalMixture`].
#'
#' @aliases ModelParamsNormal
#' @export
#'
.ModelParamsNormal <- setClass(
  Class = "ModelParamsNormal",
  slots = c(
    mean = "numeric",
    cov = "matrix",
    prec = "matrix"
  ),
  contains = "CrmPackClass",
  validity = v_model_params_normal
)

## constructor ----

#' @rdname ModelParamsNormal-class
#'
#' @param mean (`numeric`)\cr the prior mean vector.
#' @param cov (`matrix`)\cr the prior covariance matrix. The precision matrix
#'   `prec` is internally calculated as an inverse of `cov`.
#'
#' @export
#' @examples
#' ModelParamsNormal(mean = c(1, 6), cov = diag(2))
ModelParamsNormal <- function(mean, cov) {
  assert_matrix(
    cov,
    mode = "numeric",
    any.missing = FALSE,
    nrows = length(mean),
    ncols = length(mean)
  )
  # To ensure that `cov` is invertible:
  assert_true(h_is_positive_definite(cov, length(mean)))

  .ModelParamsNormal(
    mean = mean,
    cov = cov,
    prec = solve(cov)
  )
}

## default constructor ----

#' @rdname ModelParamsNormal-class
#' @note Typically, end users will not use the `.ModelPAramsNormal()` function.
#' @export
.DefaultModelParamsNormal <- function() {
  ModelParamsNormal(
    mean = c(1, 0),
    cov = matrix(c(1, 0, 0, 1), nrow = 2)
  )
}
