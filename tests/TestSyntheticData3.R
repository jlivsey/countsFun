#-------------------------------------------------------------------------#
# Purpose: Fit Synthetic data using the lgc wrapper
#
#
# Date: Jan 2023
#-------------------------------------------------------------------------#

# load necessary libraries
library(optimx)
library(ltsa)
require(countsFun)
library(itsmr)

# Specify model and methods
n              = 100
# Regressor      = cbind(rep(1,n),rbinom(n,1,0.25))
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 3
ARParm         = c(0.8, -0.25)
MAParm         = NULL
#ARParm         = NULL
#MAParm         = c(0.5,0.25)
ARMAModel      = c(length(ARParm),length(MAParm))
ParticleNumber = 10
epsilon        = 0.5
EstMethod      = "PFR"
initialParam   = NULL
TrueParam      = NULL
Task           = 'Optimization'
SampleSize     = NULL
nsim           = NULL
no_cores       = NULL
OptMethod      = "bobyqa"
OutputType     = "data.frame"
ParamScheme    = NULL
Algorithm      = "IA"

# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)


t0 = tic()
# Fit the model using the lgc wrapper
mylgc1 = lgc(DependentVar   = DependentVar,
            Regressor      = Regressor,
            EstMethod      = EstMethod,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            ParticleNumber = ParticleNumber,
            epsilon        = epsilon,
            initialParam   = initialParam,
            TrueParam      = TrueParam,
            Task           = Task,
            SampleSize     = SampleSize,
            nsim           = nsim,
            no_cores       = no_cores,
            OptMethod      = OptMethod,
            OutputType     = OutputType,
            ParamScheme    = ParamScheme
)


