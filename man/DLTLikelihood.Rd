% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Samples-methods.R
\name{DLTLikelihood}
\alias{DLTLikelihood}
\title{Likelihood of DLTs in each interval}
\usage{
DLTLikelihood(lambda, Tmax)
}
\arguments{
\item{lambda}{the vector of piecewise hazards}

\item{Tmax}{the end of the time interval for DLTs}
}
\value{
vector with the probabilities for DLTs within the intervals.
}
\description{
This is a helper function for the \code{fitPEM} methods below.
}
\keyword{internal}
