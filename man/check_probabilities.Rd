% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkmate.R
\name{check_probabilities}
\alias{check_probabilities}
\alias{assert_probabilities}
\alias{test_probabilities}
\alias{expect_probabilities}
\title{Check if an argument is a probability vector}
\usage{
check_probabilities(
  x,
  bounds_closed = TRUE,
  len = NULL,
  unique = FALSE,
  sorted = FALSE
)

assert_probabilities(
  x,
  bounds_closed = TRUE,
  len = NULL,
  unique = FALSE,
  sorted = FALSE,
  .var.name = checkmate::vname(x),
  add = NULL
)

test_probabilities(
  x,
  bounds_closed = TRUE,
  len = NULL,
  unique = FALSE,
  sorted = FALSE
)

expect_probabilities(
  x,
  bounds_closed = TRUE,
  len = NULL,
  unique = FALSE,
  sorted = FALSE,
  info = NULL,
  label = vname(x)
)
}
\arguments{
\item{x}{(\code{numeric})\cr vector or matrix with numerical values to check.}

\item{bounds_closed}{(\code{logical})\cr should bounds be closed? This can be a
scalar or vector of length two. If it is a scalar, then its value applies
equally to lower bound \eqn{0} and upper bound \eqn{1}. If this is a vector
with two flags, the first flag corresponds to the lower bound \eqn{0}
only, and the second to the upper bound \eqn{1} only.}

\item{len}{[\code{integer(1)}]\cr
Exact expected length of \code{x}.}

\item{unique}{[\code{logical(1)}]\cr
Must all values be unique? Default is \code{FALSE}.}

\item{sorted}{[\code{logical(1)}]\cr
Elements must be sorted in ascending order. Missing values are ignored.}

\item{.var.name}{[\code{character(1)}]\cr
Name of the checked object to print in assertions. Defaults to
the heuristic implemented in \code{\link[checkmate]{vname}}.}

\item{add}{[\code{AssertCollection}]\cr
Collection to store assertion messages. See \code{\link[checkmate]{AssertCollection}}.}

\item{info}{[\code{character(1)}]\cr
Extra information to be included in the message for the testthat reporter.
See \code{\link[testthat]{expect_that}}.}

\item{label}{[\code{character(1)}]\cr
Name of the checked object to print in messages. Defaults to
the heuristic implemented in \code{\link[checkmate]{vname}}.}
}
\value{
\code{TRUE} if successful, otherwise a string with the error message.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Check if every element in a given numerical vector or matrix represents a
probability, that is a number within (0, 1) interval, that can optionally be
closed at any side.
}
\note{
If there are any missing or non-finite values in \code{x}, this function
returns \code{FALSE}, regardless of the values of other elements in \code{x}.
}
\examples{
x <- c(0, 0.2, 0.1, 0.3, 1)
check_probabilities(x)
check_probabilities(x, bounds_closed = FALSE)
check_probabilities(x, bounds_closed = c(FALSE, TRUE))
}
\seealso{
\code{\link{assertions}} for more details.
}
