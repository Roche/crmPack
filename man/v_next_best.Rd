% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rules-validity.R
\name{v_next_best}
\alias{v_next_best}
\alias{v_next_best_mtd}
\alias{v_next_best_ncrm}
\alias{v_next_best_ncrm_loss}
\alias{v_next_best_dual_endpoint}
\alias{v_next_best_min_dist}
\alias{v_next_best_inf_theory}
\alias{v_next_best_td}
\alias{v_next_best_td_samples}
\alias{v_next_best_max_gain_samples}
\alias{v_next_best_prob_mtd_lte}
\alias{v_next_best_prob_mtd_min_dist}
\alias{v_next_best_ordinal}
\title{Internal Helper Functions for Validation of \code{\link{NextBest}} Objects}
\usage{
v_next_best_mtd(object)

v_next_best_ncrm(object)

v_next_best_ncrm_loss(object)

v_next_best_dual_endpoint(object)

v_next_best_min_dist(object)

v_next_best_inf_theory(object)

v_next_best_td(object)

v_next_best_td_samples(object)

v_next_best_max_gain_samples(object)

v_next_best_prob_mtd_lte(object)

v_next_best_prob_mtd_min_dist(object)

v_next_best_ordinal(object)
}
\arguments{
\item{object}{(\code{NextBest})\cr object to validate.}
}
\value{
A \code{character} vector with the validation failure messages,
or \code{TRUE} in case validation passes.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

These functions are only used internally to validate the format of an input
\code{\link{NextBest}} or inherited classes and therefore not exported.
}
\section{Functions}{
\itemize{
\item \code{v_next_best_mtd()}: validates that the \code{\link{NextBestMTD}} object
contains valid \code{target} probability and \code{derive} function.

\item \code{v_next_best_ncrm()}: validates that the \code{\link{NextBestNCRM}} object
contains valid \code{target} probability, \code{overdose} and \code{max_overdose_prob}
probability ranges.

\item \code{v_next_best_ncrm_loss()}: validates that the \code{\link{NextBestNCRMLoss}} object
contains valid objects.

\item \code{v_next_best_dual_endpoint()}: validates that the \code{\link{NextBestDualEndpoint}} object
contains valid probability objects.

\item \code{v_next_best_min_dist()}: validates that the \code{\link{NextBestMinDist}} object
contains valid \code{target} object.

\item \code{v_next_best_inf_theory()}: validates that the \code{\link{NextBestInfTheory}} object
contains valid \code{target} and \code{asymmetry} objects.

\item \code{v_next_best_td()}: validates that the \code{\link{NextBestTD}} object
contains valid \code{prob_target_drt} and \code{prob_target_eot} probabilities.

\item \code{v_next_best_td_samples()}: validates that the \code{\link{NextBestTDsamples}} object
contains valid \code{derive} function.

\item \code{v_next_best_max_gain_samples()}: validates that the \code{\link{NextBestMaxGainSamples}} object
contains valid \code{derive} and \code{mg_derive} functions.

\item \code{v_next_best_prob_mtd_lte()}: validates that the \code{\link{NextBestProbMTDLTE}} object
contains valid \code{target} probability and \code{method} string value.

\item \code{v_next_best_prob_mtd_min_dist()}: validates that the \code{\link{NextBestProbMTDMinDist}} object
contains valid \code{target} probability and \code{method} string value.

\item \code{v_next_best_ordinal()}: validates that the \code{\link{NextBestOrdinal}} object
contains valid \code{grade} and standard \code{NextBest} rule.

}}
