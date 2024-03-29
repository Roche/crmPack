# nolint start

## Simulate dose-escalation procedure based on DLE and efficacy responses where no DLE
## and efficacy samples are used
## we need a data object with doses >= 1:
data <- DataDual(doseGrid = seq(25, 300, 25), placebo = FALSE)
## First for the DLE model
## The DLE model must be of 'ModelTox' (e.g 'LogisticIndepBeta') class
DLEmodel <- LogisticIndepBeta(
  binDLE = c(1.05, 1.8),
  DLEweights = c(3, 3),
  DLEdose = c(25, 300),
  data = data
)

## The efficacy model of 'ModelEff' (e.g 'Effloglog') class
Effmodel <- Effloglog(
  eff = c(1.223, 2.513), eff_dose = c(25, 300),
  nu = c(a = 1, b = 0.025), data = data
)

## The escalation rule using the 'NextBestMaxGain' class
mynextbest <- NextBestMaxGain(
  prob_target_drt = 0.35,
  prob_target_eot = 0.3
)


## The increments (see Increments class examples)
## 200% allowable increase for dose below 300 and 200% increase for dose above 300
myIncrements <- IncrementsRelative(
  intervals = c(25, 300),
  increments = c(2, 2)
)
## cohort size of 3
mySize <- CohortSizeConst(size = 3)
## Stop only when 36 subjects are treated
myStopping <- StoppingMinPatients(nPatients = 36)
## Now specified the design with all the above information and starting with a dose of 25

## Specified the design(for details please refer to the 'DualResponsesDesign' example)
design <- DualResponsesDesign(
  nextBest = mynextbest,
  model = DLEmodel,
  eff_model = Effmodel,
  stopping = myStopping,
  increments = myIncrements,
  cohort_size = mySize,
  data = data, startingDose = 25
)
## Specify the true DLE and efficacy curves
myTruthDLE <- probFunction(DLEmodel, phi1 = -53.66584, phi2 = 10.50499)
myTruthEff <- efficacyFunction(Effmodel, theta1 = -4.818429, theta2 = 3.653058)

## The true gain curve can also be seen
myTruthGain <- function(dose) {
  return((myTruthEff(dose)) / (1 + (myTruthDLE(dose) / (1 - myTruthDLE(dose)))))
}


## Then specified the simulations and generate the trial
## For illustration purpose only 1 simulation is produced (nsim=1).
options <- McmcOptions(burnin = 100, step = 2, samples = 200)
mySim <- simulate(
  object = design,
  args = NULL,
  trueDLE = myTruthDLE,
  trueEff = myTruthEff,
  trueNu = 1 / 0.025,
  nsim = 1,
  seed = 819,
  parallel = FALSE
)

# nolint end
