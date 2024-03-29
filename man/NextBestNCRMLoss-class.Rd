% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rules-class.R
\docType{class}
\name{NextBestNCRMLoss-class}
\alias{NextBestNCRMLoss-class}
\alias{.NextBestNCRMLoss}
\alias{NextBestNCRMLoss}
\alias{.DefaultNextBestNCRMLoss}
\title{\code{NextBestNCRMLoss}}
\usage{
NextBestNCRMLoss(
  target,
  overdose,
  unacceptable = c(1, 1),
  max_overdose_prob,
  losses
)

.DefaultNextBestNCRMLoss()
}
\arguments{
\item{target}{(\code{numeric})\cr see slot definition.}

\item{overdose}{(\code{numeric})\cr see slot definition.}

\item{unacceptable}{(\code{numeric})\cr see slot definition.}

\item{max_overdose_prob}{(\code{proportion})\cr see slot definition in \code{\link{NextBestNCRM}}.}

\item{losses}{(\code{numeric})\cr see slot definition.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

\code{\link{NextBestNCRMLoss}} is the class based on NCRM rule and loss function.
This class is similar to \code{\link{NextBestNCRM}} class, but differences are the
addition of loss function and re-defined toxicity intervals, see each
toxicity interval documentation and the note for details. As in NCRM rule, first admissible doses are found,
which are those with probability to fall in overdose category being below
\code{max_overdose_prob}. Next, within the admissible doses, the loss function is
calculated, i.e. \code{losses} \%*\% \code{target}. Finally, the corresponding
dose with lowest loss function (Bayes risk) is recommended for the next dose.
}
\section{Slots}{

\describe{
\item{\code{target}}{(\code{numeric})\cr the target toxicity interval (limits included).
It has to be a probability range excluding 0 and 1.}

\item{\code{overdose}}{(\code{numeric})\cr the overdose toxicity interval (lower limit
excluded, upper limit included) or the excessive toxicity interval (lower
limit excluded, upper limit included) if unacceptable is not provided.
It has to be a probability range. It is used to filter probability samples.}

\item{\code{unacceptable}}{(\code{numeric})\cr an unacceptable toxicity
interval (lower limit excluded, upper limit included). This must be
specified if the \code{overdose} does not include 1. Otherwise, it is \code{c(1, 1)}
(default), which is essentially a scalar equals 1. It has to be a
probability range.}

\item{\code{losses}}{(\code{numeric})\cr a vector specifying the loss function. If the
\code{unacceptable} is provided, the vector length must be 4, otherwise 3.}
}}

\note{
The loss function should be a vector of either 3 or 4 values.
This is because the loss function values must be specified for each
interval, that is under-dosing, target toxicity, and overdosing toxicity or
under-dosing, target toxicity, overdosing (excessive) toxicity, and
unacceptable toxicity intervals.

Typically, end users will not use the \code{.DefaultNextBestnCRMLoss()} function.
}
\examples{
# In the example below, the target toxicity interval [0.2, 0.35] while the
# overdose interval is (0.35, 1]. We would like to constrain the probability
# of overdosing below 25\%. The loss function is c(1, 0, 1, 2).
my_next_best <- NextBestNCRMLoss(
  target = c(0.2, 0.35),
  overdose = c(0.35, 0.6),
  unacceptable = c(0.6, 1),
  max_overdose_prob = 0.25,
  losses = c(1, 0, 1, 2)
)
}
