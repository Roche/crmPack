# nolint start

# Create the data
data <- DataDual(
  x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10,
      20, 20, 20, 40, 40, 40, 50, 50, 50),
  y=c(0, 0, 0, 0, 0, 0, 1, 0,
      0, 1, 1, 0, 0, 1, 0, 1, 1),
  w=c(0.31, 0.42, 0.59, 0.45, 0.6, 0.7, 0.55, 0.6,
      0.52, 0.54, 0.56, 0.43, 0.41, 0.39, 0.34, 0.38, 0.21),
  doseGrid=c(0.1, 0.5, 1.5, 3, 6,
             seq(from=10, to=80, by=2)))

# Initialize the Dual-Endpoint model (in this case RW1)
model <- DualEndpointRW(mean = c(0, 1),
                        cov = matrix(c(1, 0, 0, 1), nrow=2),
                        sigma2betaW = 0.01,
                        sigma2W = c(a=0.1, b=0.1),
                        rho = c(a=1, b=1),
                        rw1 = TRUE)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=500)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# In this case target a dose achieving at least 0.9 of maximum biomarker level (efficacy)
# and with a probability below 0.25 that prob(DLT)>0.35 (safety)
myNextBest <- NextBestDualEndpoint(target=c(0.9, 1),
                                   overdose=c(0.35, 1),
                                   max_overdose_prob=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples,
                               model=model,
                               data=data)

# Define the stopping rule such that the study would be stopped if if there is at
# least 0.5 posterior probability that the biomarker (efficacy) is within the
# biomarker target range of [0.9, 1.0] (relative to the maximum for the biomarker).
myStopping <- StoppingTargetBiomarker(target = c(0.9, 1),
                                      prob = 0.5)

# Evaluate if to stop the trial
stopTrial(stopping=myStopping,
          dose=doseRecommendation$value,
          samples=samples,
          model=model,
          data=data)

# nolint end
