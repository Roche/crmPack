% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rules-class.R
\docType{class}
\name{Stopping-class}
\alias{Stopping-class}
\alias{Stopping}
\title{\code{Stopping}}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

\code{\link{Stopping}} is a class for stopping rules.
}
\section{Slots}{

\describe{
\item{\code{report_label}}{(\code{string})\cr a label for the stopping report. The meaning
of this parameter is twofold. If it is equal to \code{NA_character_} (default),
the \code{report_label} will not be used in the report at all. Otherwise, if it
is specified as an empty character (i.e. \code{character(0)}) in a user constructor,
then a default, class-specific label will be created for this slot.
Finally, for the remaining cases, a user can provide a custom label.}
}}

\seealso{
\code{\link{StoppingList}}, \code{\link{StoppingCohortsNearDose}}, \code{\link{StoppingPatientsNearDose}},
\code{\link{StoppingMinCohorts}}, \code{\link{StoppingMinPatients}}, \code{\link{StoppingTargetProb}},
\code{\link{StoppingMTDdistribution}}, \code{\link{StoppingTargetBiomarker}}, \code{\link{StoppingHighestDose}}
\code{\link{StoppingMTDCV}}, \code{\link{StoppingLowestDoseHSRBeta}}, \code{\link{StoppingSpecificDose}}.
}
