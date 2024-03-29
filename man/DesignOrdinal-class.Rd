% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Design-class.R
\docType{class}
\name{DesignOrdinal-class}
\alias{DesignOrdinal-class}
\alias{.DesignOrdinal}
\alias{DesignOrdinal}
\alias{.DefaultDesignOrdinal}
\title{\code{DesignOrdinal}}
\usage{
DesignOrdinal(
  model,
  stopping,
  increments,
  pl_cohort_size = CohortSizeOrdinal(1L, CohortSizeConst(0L)),
  ...
)

.DefaultDesignOrdinal()
}
\arguments{
\item{model}{(\code{LogisticLogNormalOrdinal})\cr see slot definition.}

\item{stopping}{(\code{StoppingOrdinal})\cr see slot definition.}

\item{increments}{(\code{IncrementsOrdinal})\cr see slot definition.}

\item{pl_cohort_size}{(\code{CohortSizeOrdinal})\cr see slot definition.}

\item{...}{
  Arguments passed on to \code{\link[=RuleDesignOrdinal]{RuleDesignOrdinal}}
  \describe{
    \item{\code{next_best}}{(\code{NextBestOrdinal})\cr see slot definition.}
    \item{\code{cohort_size}}{(\code{CohortSizeOrdinal})\cr see slot definition.}
    \item{\code{data}}{(\code{DataOrdinal})\cr see slot definition.}
    \item{\code{starting_dose}}{(\code{number})\cr see slot definition.}
  }}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

\code{\link{DesignOrdinal}} is the class for rule-based ordinal designs. The difference
between this class and its parent \code{\link{RuleDesignOrdinal}} class is that the
\code{\link{DesignOrdinal}} class contains additional \code{model}, \code{stopping},
\code{increments} and \code{pl_cohort_size} slots.
}
\section{Slots}{

\describe{
\item{\code{model}}{(\code{LogisticLogNormalOrdinal})\cr the model to be used.}

\item{\code{stopping}}{(\code{StoppingOrdinal})\cr stopping rule(s) for the trial.}

\item{\code{increments}}{(\code{IncrementsOrdinal})\cr how to control increments between dose levels.}

\item{\code{pl_cohort_size}}{(\code{CohortSizeOrdinal})\cr rules for the cohort sizes for placebo,
if any planned (defaults to constant 0 placebo patients).}
}}

\note{
Typically, end users will not use the \code{.DefaultDesignOrdinal()} function.
}
\examples{
my_size1 <- CohortSizeRange(
  intervals = c(0, 30),
  cohort_size = c(1, 3)
)
my_size2 <- CohortSizeDLT(
  intervals = c(0, 1),
  cohort_size = c(1, 3)
)
my_size <- CohortSizeOrdinal(1L, maxSize(my_size1, my_size2))

my_stopping1 <- StoppingMinCohorts(nCohorts = 3)
my_stopping2 <- StoppingTargetProb(
  target = c(0.2, 0.35),
  prob = 0.5
)
my_stopping3 <- StoppingMinPatients(nPatients = 20)
my_stopping <- StoppingOrdinal(1L, (my_stopping1 & my_stopping2) | my_stopping3)

# Initialize the design.
design <- DesignOrdinal(
  model = LogisticLogNormalOrdinal(
    mean = c(-3, -4, 1),
    cov = diag(c(3, 4, 1)),
    ref_dose = 50
  ),
  next_best = NextBestOrdinal(
    1L,
    NextBestNCRM(
      target = c(0.2, 0.35),
      overdose = c(0.35, 1),
      max_overdose_prob = 0.25
    )
  ),
  stopping = my_stopping,
  increments = IncrementsOrdinal(
    1L,
    IncrementsRelative(
      intervals = c(0, 20),
      increments = c(1, 0.33)
    )
  ),
  cohort_size = my_size,
  data = DataOrdinal(
    doseGrid = c(1, 3, 5, 10, 15, 20, 25, 40, 50, 80, 100),
    yCategories = c("No tox" = 0L, "Sub-tox AE" = 1L, "DLT" = 2L)
  ),
  starting_dose = 3
)
}
