# nolint start

# Define the dose-grid and PEM parameters
emptydata <- DataDA(doseGrid = c(
  0.1, 0.5, 1, 1.5, 3, 6,
  seq(from = 10, to = 80, by = 2)
), Tmax = 60)
# Initialize the mDA-CRM model
npiece_ <- 10
Tmax_ <- 60

lambda_prior <- function(k) {
  npiece_ / (Tmax_ * (npiece_ - k + 0.5))
}

model <- DALogisticLogNormal(
  mean = c(-0.85, 1),
  cov = matrix(c(1, -0.5, -0.5, 1), nrow = 2),
  ref_dose = 56,
  npiece = npiece_,
  l = as.numeric(t(apply(as.matrix(c(1:npiece_), 1, npiece_), 2, lambda_prior))),
  c_par = 2
)
# Choose the rule for dose increments
myIncrements <- IncrementsRelative(
  intervals = c(0, 20),
  increments = c(1, 0.33)
)

myNextBest <- NextBestNCRM(
  target = c(0.2, 0.35),
  overdose = c(0.35, 1),
  max_overdose_prob = 0.25
)

# Choose the rule for the cohort-size
mySize1 <- CohortSizeRange(
  intervals = c(0, 30),
  cohort_size = c(1, 3)
)
mySize2 <- CohortSizeDLT(
  intervals = c(0, 1),
  cohort_size = c(1, 3)
)
mySize <- maxSize(mySize1, mySize2)

# Choose the rule for stopping
myStopping1 <- StoppingTargetProb(
  target = c(0.2, 0.35),
  prob = 0.5
)
myStopping2 <- StoppingMinPatients(nPatients = 50)

myStopping <- (myStopping1 | myStopping2)

# Choose the safety window
mysafetywindow <- SafetyWindowConst(c(6, 2), 7, 7)

# Initialize the design
design <- DADesign(
  model = model,
  increments = myIncrements,
  nextBest = myNextBest,
  stopping = myStopping,
  cohort_size = mySize,
  data = emptydata,
  safetyWindow = mysafetywindow,
  startingDose = 3
)

set.seed(4235)
# MCMC parameters are set to small values only to show this example. They should be
# increased for a real case.
# This procedure will take a while.
options <- McmcOptions(
  burnin = 10,
  step = 1,
  samples = 100,
  rng_kind = "Mersenne-Twister",
  rng_seed = 12
)
testthat::expect_warning(
  result <- examine(design, mcmcOptions = options, maxNoIncrement = 2),
  "Stopping because 2 times no increment"
)

# nolint end
