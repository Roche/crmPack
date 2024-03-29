% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data-methods.R
\name{getEff}
\alias{getEff}
\alias{getEff,DataDual-method}
\alias{getEff-DataDual}
\title{Extracting Efficacy Responses for Subjects Categorized by the DLT}
\usage{
getEff(object, ...)

\S4method{getEff}{DataDual}(object, no_dlt = FALSE)
}
\arguments{
\item{object}{(\code{DataDual})\cr object from which the responses and dose levels
are extracted.}

\item{...}{further arguments passed to class-specific methods.}

\item{no_dlt}{(\code{flag})\cr should only no DLT responses be returned? Otherwise,
all responses are returned.}
}
\value{
\code{list} with efficacy responses categorized by the DLT value.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

A method that extracts efficacy responses for subjects and categorizes it
with respect to DLT, i.e. DLT or no DLT. The efficacy responses
are reported together with their corresponding dose levels.
}
\examples{
# Example data.
data <- DataDual(
  x = c(25, 50, 25, 50, 75, 300, 250, 150),
  y = c(0, 0, 0, 0, 0, 1, 1, 0),
  w = c(0.31, 0.42, 0.59, 0.45, 0.6, 0.7, 0.6, 0.52),
  doseGrid = seq(25, 300, 25)
)

# Get the efficacy response and their corresponding dose levels
# categorized by the DLT.
getEff(data)
}
