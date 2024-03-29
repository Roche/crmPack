% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rules-methods.R
\name{or-StoppingAny-Stopping}
\alias{or-StoppingAny-Stopping}
\alias{|,Stopping,StoppingAny-method}
\title{The method combining an atomic and a stopping list}
\usage{
\S4method{|}{Stopping,StoppingAny}(e1, e2)
}
\arguments{
\item{e1}{\code{\linkS4class{Stopping}} object}

\item{e2}{\code{\linkS4class{StoppingAny}} object}
}
\value{
The modified \code{\linkS4class{StoppingAny}} object
}
\description{
The method combining an atomic and a stopping list
}
\examples{

## Example of combining an atomic stopping rule with a list of stopping rules
## with an OR ('|') operator

myStopping1 <- StoppingMinCohorts(nCohorts=3)
myStopping2 <- StoppingTargetProb(target=c(0.2, 0.35),
                                  prob=0.5)

myStopping3 <- StoppingMinPatients(nPatients=20)

myStopping <-  myStopping3 | (myStopping1 & myStopping2 ) 




}
\keyword{methods}
