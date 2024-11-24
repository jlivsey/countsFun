# load libraries
library(countsFun)
library(tictoc)
library(optimx)
library(ltsa)
library(itsmr)
library(numDeriv)
library(MASS)
library(parallel)
library(doParallel)

# Specify model and methods
n              = 50
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 3
ARParm         = 0.5
MAParm         = NULL
ARMAModel      = c(length(ARParm), length(MAParm))
ParticleNumber = 5
epsilon        = 0.5
EstMethod      = "PFR"
TrueParam      = c(MargParm, ARParm, MAParm)
Task           = 'Simulation'
SampleSize     = 50
nsim           = 2
no_cores       = 2
OptMethod      = "bobyqa"
OptMethod      = "L-BFGS-B"
OutputType     = "list"
ParamScheme    = NULL
maxdiff        = 10 ^ (-6)

DependentVar =
  c(5, 4, 2, 4, 1, 3, 3, 4, 4, 7, 4, 5, 7, 5, 1, 2, 2, 3, 3, 4,
    4, 5, 4, 2, 2, 1, 1, 1, 1, 2, 0, 0, 3, 4, 6, 5, 4, 3, 2, 1, 4,
    3, 5, 3, 1, 1, 3, 4, 6, 3)


# save the data in a data frame
df = data.frame(DependentVar)

# specify the regression model
formula = DependentVar~0

# call the wrapper
a = lgc(
  formula,
  df,
  EstMethod,
  CountDist,
  ARMAModel,
  ParticleNumber,
  epsilon,
  initialParam = TrueParam,
  TrueParam,
  Task = "Simulation",
  SampleSize,
  nsim ,
  no_cores,
  OptMethod,
  OutputType,
  ParamScheme,
  maxdiff
)
