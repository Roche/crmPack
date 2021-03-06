% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model-class.R
\docType{class}
\name{LogisticLogNormalMixture-class}
\alias{LogisticLogNormalMixture-class}
\alias{.LogisticLogNormalMixture}
\title{Standard logistic model with online mixture of two bivariate log normal priors}
\description{
This model can be used when data is arising online from the informative
component of the prior, at the same time with the data of the trial of
main interest. Formally, this is achieved by assuming that the probability
of a DLT at dose \eqn{x} is given by
}
\details{
\deqn{p(x) = \pi p_{1}(x) + (1 - \pi) p_{2}(x)}

where \eqn{\pi} is the probability for the model \eqn{p(x)} being the same
as the model \eqn{p_{1}(x)} - this is 
the informative component of the prior. From this model data arises in 
parallel: at doses \code{xshare}, DLT information \code{yshare} is observed, 
in total \code{nObsshare} data points, see \code{\linkS4class{DataMixture}}.
On the other hand, \eqn{1 - \pi}
is the probability of a separate model \eqn{p_{2}(x)}. Both components 
have the same log normal prior distribution, which can be specified by the
user, and which is inherited from the \code{\linkS4class{LogisticLogNormal}}
class.
}
\section{Slots}{

\describe{
\item{\code{shareWeight}}{the prior weight for sharing the same model \eqn{p_{1}(x)}}
}}

\examples{
## decide on the dose grid:
doseGrid <- 1:80

## and MCMC options:
options <- McmcOptions()

## the classic model would be:
model <- LogisticLogNormal(mean = c(-0.85, 1),
                           cov = matrix(c(1, -0.5, -0.5, 1), nrow = 2),
                           refDose = 50)

nodata <- Data(doseGrid=doseGrid)

priorSamples <- mcmc(nodata, model, options)
plot(priorSamples, model, nodata)

## set up the mixture model and data share object:
modelShare <- LogisticLogNormalMixture(shareWeight=0.1,
                                       mean = c(-0.85, 1),
                                       cov = matrix(c(1, -0.5, -0.5, 1), nrow = 2),
                                       refDose = 50)

nodataShare <- DataMixture(doseGrid=doseGrid,
                           xshare=
                             c(rep(10, 4),
                               rep(20, 4),
                               rep(40, 4)),
                           yshare=
                             c(rep(0L, 4),
                               rep(0L, 4),
                               rep(0L, 4)))

## now compare with the resulting prior model:
priorSamplesShare <- mcmc(nodataShare, modelShare, options)
plot(priorSamplesShare, modelShare, nodataShare)
}
\seealso{
the \code{\linkS4class{DataMixture}} class for use with this model
}
\keyword{classes}
