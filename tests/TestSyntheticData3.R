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

# Specify model and methods
n              = 100
Regressor      = cbind(rep(1,n),rbinom(n,1,0.25))
CountDist      = "Poisson"
beta0          = 2
beta1          = 0.4
MargParm       = c(beta0, beta1)
ARParm         = c(0.75)
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))
ParticleNumber = 10
epsilon        = 0.5

# simulate data
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

# Fit the model using the lgc wrapper
mylgc = lgc(DependentVar   = DependentVar,
            Regressor      = Regressor,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            ParticleNumber = ParticleNumber,
            epsilon        = epsilon
)























