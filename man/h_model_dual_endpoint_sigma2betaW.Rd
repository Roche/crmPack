% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_model.R
\name{h_model_dual_endpoint_sigma2betaW}
\alias{h_model_dual_endpoint_sigma2betaW}
\title{Update certain components of \code{\link{DualEndpoint}} model with regard to prior variance
factor of the random walk.}
\usage{
h_model_dual_endpoint_sigma2betaW(use_fixed, sigma2betaW, de)
}
\arguments{
\item{use_fixed}{(\code{flag})\cr indicates whether a fixed value for
\code{sigma2betaW} should be used or not. If \code{sigma2betaW} is not supposed
to be a fixed value, a prior distribution from the Inverse-Gamma distribution
will be used. See the details below, under \code{sigma2betaW} argument.}

\item{sigma2betaW}{(\code{numeric})\cr the prior variance factor of the random walk
prior for the biomarker model. Either a fixed value or Inverse-Gamma distribution
parameters, i.e. vector with two elements named \code{a} and \code{b}.}

\item{de}{(\code{DualEnpoint})\cr dual endpoint model whose slots will be updated.}
}
\value{
A \code{\link{DualEndpoint}} model with updated \code{priormodel}, \code{modelspecs},
\code{init}, \code{sample} slots.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

A simple helper function that takes \code{\link{DualEndpoint}} object and updates
\code{priormodel}, \code{modelspecs}, \code{init}, \code{sample} slots according to the random walk
variance.
}
\seealso{
\code{\link{DualEndpointRW}}.
}
