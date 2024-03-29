% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_jags.R
\name{h_jags_get_data}
\alias{h_jags_get_data}
\title{Getting Data for \code{JAGS}}
\usage{
h_jags_get_data(model, data, from_prior)
}
\arguments{
\item{model}{(\code{GeneralModel})\cr an input model.}

\item{data}{(\code{GeneralData})\cr an input data.}

\item{from_prior}{(\code{flag})\cr sample from the prior only? In this case
data will not be appended to the output, i.e. only the variables required
by the \code{model@priormodel} model will be returned in data.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

A simple helper function that prepares an object for \code{data} argument of
\code{\link[rjags:jags.model]{rjags::jags.model()}}, which is invoked by \code{\link[=mcmc]{mcmc()}} method.
}
\examples{
# Create some data from the class `Data`.
my_data <- Data(
  x = c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
  y = c(0, 0, 0, 0, 0, 0, 1, 0),
  doseGrid = c(0.1, 0.5, 1.5, 3, 6, seq(from = 10, to = 80, by = 2))
)

# Initialize the CRM model.
my_model <- LogisticLogNormal(
  mean = c(-0.85, 1),
  cov = matrix(c(1, -0.5, -0.5, 1), nrow = 2),
  ref_dose = 56
)

jags_data <- h_jags_get_data(my_model, my_data, from_prior = FALSE)
jags_data
}
