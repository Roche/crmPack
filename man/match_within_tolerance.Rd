% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{match_within_tolerance}
\alias{match_within_tolerance}
\title{Helper function for value matching with tolerance}
\usage{
match_within_tolerance(x, table)
}
\arguments{
\item{x}{the values to be matched}

\item{table}{the values to be matched against}
}
\value{
A vector of the same length as \code{x} or
empty vector if \code{table} is empty.
}
\description{
This is a modified version of \code{match} that supports tolerance.
}
\keyword{programming}
