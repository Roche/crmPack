% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rules-class.R
\docType{class}
\name{StoppingTargetBiomarker-class}
\alias{StoppingTargetBiomarker-class}
\alias{.StoppingTargetBiomarker}
\title{Stop based on probability of target biomarker}
\description{
Stop based on probability of target biomarker
}
\section{Slots}{

\describe{
\item{\code{target}}{the biomarker target range, that
needs to be reached. For example, (0.8, 1.0) and \code{scale="relative"} 
means we target a dose with at least 80\% of maximum biomarker level.}

\item{\code{scale}}{either \code{relative} (default, then the \code{target} is interpreted 
relative to the maximum, so must be a probability range) or \code{absolute}
(then the \code{target} is interpreted as absolute biomarker range)}

\item{\code{prob}}{required target probability for reaching sufficient precision}
}}

\examples{
# As example, here is the rule for: 
#   stopping the study if there is at least 0.5 probability that the biomarker 
#   (efficacy) is within the biomarker target range of [0.9, 1.0] (relative to the
#   maximum for the biomarker).

myStopping <- StoppingTargetBiomarker(target = c(0.9, 1),
                                      prob = 0.5)

}
\keyword{classes}
