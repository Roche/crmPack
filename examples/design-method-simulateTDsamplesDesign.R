# nolint start

## Simulate dose-escalation procedure based only on DLE responses with DLE samples involved

## The design comprises a model, the escalation rule, starting data,
## a cohort size and a starting dose
## Define your data set first using an empty data set
## with dose levels from 25 to 300 with increments 25
data <- Data(doseGrid = seq(25, 300, 25))

## The design only incorporate DLE responses and DLE samples are involved
## Specified the model of 'ModelTox' class eg 'LogisticIndepBeta' class model
model <- LogisticIndepBeta(
  binDLE = c(1.05, 1.8),
  DLEweights = c(3, 3),
  DLEdose = c(25, 300),
  data = data
)
## Then the escalation rule
tdNextBest <- NextBestTDsamples(
  prob_target_drt = 0.35,
  prob_target_eot = 0.3,
  derive = function(samples) {
    as.numeric(quantile(samples, probs = 0.3))
  }
)

## The cohort size, size of 3 subjects
mySize <- CohortSizeConst(size = 3)
## Deifne the increments for the dose-escalation process
## The maximum increase of 200% for doses up to the maximum of the dose specified in the doseGrid
## The maximum increase of 200% for dose above the maximum of the dose specified in the doseGrid
## This is to specified a maximum of 3-fold restriction in dose-esclation
myIncrements <- IncrementsRelative(
  intervals = c(min(data@doseGrid), max(data@doseGrid)),
  increments = c(2, 2)
)
## Specified the stopping rule e.g stop when the maximum sample size of 36 patients has been reached
myStopping <- StoppingMinPatients(nPatients = 36)

## Specified the design(for details please refer to the 'TDsamplesDesign' example)
design <- TDsamplesDesign(
  model = model,
  nextBest = tdNextBest,
  stopping = myStopping,
  increments = myIncrements,
  cohort_size = mySize,
  data = data, startingDose = 25
)

## Specify the truth of the DLE responses
myTruth <- probFunction(model, phi1 = -53.66584, phi2 = 10.50499)

## then plot the truth to see how the truth dose-DLE curve look like
curve(myTruth(x), from = 0, to = 300, ylim = c(0, 1))

## Then specified the simulations and generate the trial
## options for MCMC
options <- McmcOptions(burnin = 100, step = 2, samples = 200)
## The simulations
## For illustration purpose only 1 simulation is produced (nsim=1).
mySim <- simulate(
  object = design,
  args = NULL,
  truth = myTruth,
  nsim = 1,
  seed = 819,
  mcmcOptions = options,
  parallel = FALSE
)

# nolint end
