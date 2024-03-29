% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model-class.R
\docType{class}
\name{TITELogisticLogNormal-class}
\alias{TITELogisticLogNormal-class}
\alias{.TITELogisticLogNormal}
\alias{TITELogisticLogNormal}
\alias{.DefaultTITELogisticLogNormal}
\title{\code{TITELogisticLogNormal}}
\usage{
TITELogisticLogNormal(weight_method = "linear", ...)

.DefaultTITELogisticLogNormal()
}
\arguments{
\item{weight_method}{(\code{string})\cr the weight function method: either linear
or adaptive. This was used in Liu, Yin and Yuan's paper.}

\item{...}{
  Arguments passed on to \code{\link[=LogisticLogNormal]{LogisticLogNormal}}
  \describe{
    \item{\code{mean}}{(\code{numeric})\cr the prior mean vector.}
    \item{\code{cov}}{(\code{matrix})\cr the prior covariance matrix. The precision matrix
\code{prec} is internally calculated as an inverse of \code{cov}.}
    \item{\code{ref_dose}}{(\code{number})\cr the reference dose \eqn{x*} (strictly positive
number).}
  }}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

\code{\link{TITELogisticLogNormal}} is the class for TITE-CRM based on a logistic
regression model using a bivariate normal prior on the intercept and log
slope parameters.

This class inherits from the \code{\link{LogisticLogNormal}}.
}
\section{Slots}{

\describe{
\item{\code{weight_method}}{(\code{string})\cr the weight function method: either linear
or adaptive. This was used in Liu, Yin and Yuan's paper.}
}}

\note{
Typically, end users will not use the \code{.DefaultTITELogisticLogNormal()} function.
}
\examples{
my_model <- TITELogisticLogNormal(
  mean = c(0, 1),
  cov = diag(2),
  ref_dose = 1,
  weight_method = "linear"
)

my_model1 <- TITELogisticLogNormal(
  mean = c(0, 1),
  cov = diag(2),
  ref_dose = 1,
  weight_method = "adaptive"
)
}
\seealso{
\code{\link{DALogisticLogNormal}}.
}
