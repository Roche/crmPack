% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data-methods.R
\name{update,DataDual-method}
\alias{update,DataDual-method}
\alias{update-DataDual}
\title{Updating \code{DataDual} Objects}
\usage{
\S4method{update}{DataDual}(object, w, ..., check = TRUE)
}
\arguments{
\item{object}{(\code{DataDual})\cr object you want to update.}

\item{w}{(\code{numeric})\cr the continuous vector of biomarker values
for all the patients in this update.}

\item{...}{further arguments passed to \code{Data} update method \code{\link{update-Data}}.}

\item{check}{(\code{flag})\cr whether the validation of the updated object
should be conducted. See help for \code{\link{update-Data}} for more details
on the use case of this parameter.}
}
\value{
The new, updated \code{\link{DataDual}} object.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

A method that updates existing \code{\link{DataDual}} object with new data.
}
\examples{
# Create some data of class 'DataDual'.
my_data <- DataDual(
  x = c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
  y = c(0, 0, 0, 0, 0, 0, 1, 0),
  w = rnorm(8),
  doseGrid = c(0.1, 0.5, 1.5, 3, 6, seq(from = 10, to = 80, by = 2))
)

# Update the data with a new cohort.
my_data1 <- update(my_data, x = 20, y = c(0, 1, 1), w = c(0.4, 1.2, 2.2))
my_data1
}
