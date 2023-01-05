# nolint start

##If only DLE responses are considered in the simulations
##Specified your simulations when no DLE samples are used
## data set with dose levels from 25 to 300 with increments 25
data <- Data(doseGrid=seq(25,300,25))

##The design only incorporate DLE responses and DLE samples are involved
##Specified the model of 'ModelTox' class eg 'LogisticIndepBeta' class model
model<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),DLEdose=c(25,300),data=data)
##Then the escalation rule
tdNextBest <- NextBestTD(prob_target_drt=0.35,
                         prob_target_eot=0.3)

##Then the starting data, an empty data set
emptydata<-Data(doseGrid=seq(25,300,25))
## The cohort size, size of 3 subjects
mySize <-CohortSizeConst(size=3)
##Deifne the increments for the dose-escalation process
##The maximum increase of 200% for doses up to the maximum of the dose specified in the doseGrid
##The maximum increase of 200% for dose above the maximum of the dose specified in the doseGrid
##This is to specified a maximum of 3-fold restriction in dose-esclation
myIncrements<-IncrementsRelative(intervals=c(min(data@doseGrid),max(data@doseGrid)),
                                 increments=c(2,2))
##Specified the stopping rule e.g stop when the maximum sample size of 36 patients has been reached
myStopping <- StoppingMinPatients(nPatients=36)


##Specified the design(for details please refer to the 'TDDesign' example)
design <- TDDesign(model=model,
                   nextBest=tdNextBest,
                   stopping=myStopping,
                   increments=myIncrements,
                   cohortSize=mySize,
                   data=data,startingDose=25)

##Specify the truth of the DLE responses
myTruth <- probFunction(model, phi1 = -53.66584, phi2 = 10.50499)

##(Please refer to desgin-method 'simulate TDDesign' examples for details)
##For illustration purpose only 1 simulation is produced (nsim=1).
mySim <- simulate(design,
                  args=NULL,
                  truth=myTruth,
                  nsim=1,
                  seed=819,
                  parallel=FALSE)
##Then produce a summary of your simulations
summary(mySim,
        truth=myTruth)

##If DLE samples are involved
##specify the next best
tdNextBest <- NextBestTDsamples(
  prob_target_drt = 0.35,
  prob_target_eot = 0.3,
  derive = function(samples) {
    as.numeric(quantile(samples, probs = 0.3))
  }
)
##The design
design <- TDsamplesDesign(model=model,
                          nextBest=tdNextBest,
                          stopping=myStopping,
                          increments=myIncrements,
                          cohortSize=mySize,
                          data=data,startingDose=25)
##options for MCMC
##For illustration purpose, we will use 50 burn-ins to generate 200 samples
options<-McmcOptions(burnin=50,step=2,samples=200)
##The simulations
## For illustration purpose we will only generate 2 trials (nsim=2)
mySim <- simulate(design,
                  args=NULL,
                  truth=myTruth,
                  nsim=2,
                  seed=819,
                  mcmcOptions=options,
                  parallel=FALSE)
##Then produce a summary of your simulations
summary(mySim,
        truth=myTruth)

# nolint end