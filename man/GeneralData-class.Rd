% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data-class.R
\docType{class}
\name{GeneralData-class}
\alias{GeneralData-class}
\alias{.GeneralData}
\alias{GeneralData}
\alias{.DefaultDataGeneral}
\title{\code{GeneralData}}
\usage{
.DefaultDataGeneral()
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

\code{\link{GeneralData}} is a class for general data input.
}
\section{Slots}{

\describe{
\item{\code{ID}}{(\code{integer})\cr unique patient IDs.}

\item{\code{cohort}}{(\code{integer})\cr the cohort (non-negative sorted) indices.}

\item{\code{nObs}}{(\code{integer})\cr number of observations, a single value.}
}}

\note{
Typically, end users will not use the \code{.DefaultDataGeneral()} function.
}
