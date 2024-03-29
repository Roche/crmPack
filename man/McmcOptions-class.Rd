% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/McmcOptions-class.R
\docType{class}
\name{McmcOptions-class}
\alias{McmcOptions-class}
\alias{.McmcOptions}
\alias{McmcOptions}
\alias{.DefaultMcmcOptions}
\title{\code{McmcOptions}}
\usage{
McmcOptions(
  burnin = 10000L,
  step = 2L,
  samples = 10000L,
  rng_kind = NA_character_,
  rng_seed = NA_integer_
)

.DefaultMcmcOptions()
}
\arguments{
\item{burnin}{(\code{count})\cr number of burn-in iterations which are not saved.}

\item{step}{(\code{count})\cr only every step-th iteration is saved after
the burn-in.}

\item{samples}{(\code{count})\cr number of resulting samples.}

\item{rng_kind}{(\code{string})\cr the name of the RNG type. Possible types are:
\code{Wichmann-Hill}, \code{Marsaglia-Multicarry}, \code{Super-Duper}, \code{Mersenne-Twister}.
If it is \code{NA} (default), then the RNG kind will be chosen by \verb{[rjags}].}

\item{rng_seed}{(\code{number})\cr RNG seed corresponding to chosen \code{rng_kind}.
It must be an integer value or \code{NA} (default), which means that the seed
will be chosen by \verb{[rjags}].}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

\code{\link{McmcOptions}} is a class for the three canonical MCMC options as well as
Random Number Generator settings.
}
\section{Slots}{

\describe{
\item{\code{iterations}}{(\code{count})\cr number of MCMC iterations.}

\item{\code{burnin}}{(\code{count})\cr number of burn-in iterations which are not saved.}

\item{\code{step}}{(\code{count})\cr only every \code{step}-th iteration is saved after
the \code{burnin}. In other words, a sample from iteration
\verb{i = 1,...,iterations}, is saved if and only if
\verb{(i - burnin) mod step = 0}.\cr
For example, for \code{iterations = 6}, \code{burnin = 0} and \code{step = 2}, only
samples from iterations \verb{2,4,6} will be saved.}

\item{\code{rng_kind}}{(\code{string})\cr a Random Number Generator (RNG) type used by
\code{\link{rjags}}. It must be one out of the following four values:
\code{base::Wichmann-Hill}, \code{base::Marsaglia-Multicarry},
\code{base::Super-Duper}, \code{base::Mersenne-Twister}, or \code{NA_character_}.
If it is \code{NA_character_} (default), then the RNG kind will be chosen by
\code{\link{rjags}}.}

\item{\code{rng_seed}}{(\code{number})\cr a Random Number Generator (RNG) seed
used by \code{\link{rjags}} for a chosen \code{rng_kind}. It must be an integer scalar or
\code{NA_integer_}, which means that the seed will be chosen by \code{\link{rjags}}.}
}}

\note{
Typically, end users will not use the \code{.DefaultMcmcOptions()} function.
}
\examples{
# Set up MCMC option in order to have a burn-in of 10000 iterations and
# then take every other iteration up to a collection of 10000 samples.
McmcOptions(burnin = 10000, step = 2, samples = 10000)
}
