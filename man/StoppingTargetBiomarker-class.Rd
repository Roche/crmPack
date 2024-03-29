% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rules-class.R
\docType{class}
\name{StoppingTargetBiomarker-class}
\alias{StoppingTargetBiomarker-class}
\alias{.StoppingTargetBiomarker}
\alias{StoppingTargetBiomarker}
\alias{.DefaultStoppingTargetBiomarker}
\title{\code{StoppingTargetBiomarker}}
\usage{
StoppingTargetBiomarker(
  target = c(0.9, 1),
  prob = 0.3,
  is_relative = TRUE,
  report_label = NA_character_
)

.DefaultStoppingTargetBiomarker()
}
\arguments{
\item{target}{(\code{numeric})\cr see slot definition.}

\item{prob}{(\code{proportion})\cr see slot definition.}

\item{is_relative}{(\code{flag})\cr see slot definition.}

\item{report_label}{(\code{string} or \code{NA})\cr see slot definition.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

\code{\link{StoppingTargetBiomarker}} is a class for stopping based on probability of
target biomarker.
}
\section{Slots}{

\describe{
\item{\code{target}}{(\code{numeric})\cr the biomarker target range that needs to be
reached. For example, \code{target = c(0.8, 1.0)} with \code{is_relative = TRUE}
means that we target a dose with at least 80\% of maximum biomarker level.}

\item{\code{is_relative}}{(\code{flag})\cr is target relative? If it so (default), then
the \code{target} is interpreted relative to the maximum, so it must be a
probability range. Otherwise, the \code{target} is interpreted as absolute
biomarker range.}

\item{\code{prob}}{(\code{proportion})\cr required target probability (except 0 or 1) for
reaching sufficient precision.}
}}

\note{
Typically, end users will not use the \code{.DefaultStoppingTargetBiomarker()} function.
}
\examples{

# Stopping the study if there is at least 0.5 probability that the biomarker
# (efficacy) is within the biomarker target range of [0.9, 1.0] (relative to the
# maximum for the biomarker).

my_stopping <- StoppingTargetBiomarker(target = c(0.9, 1), prob = 0.5)
}
