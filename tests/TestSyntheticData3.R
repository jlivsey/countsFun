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
# Regressor      = NULL
Regressor      = cbind(rep(1,n),rbinom(n,1,0.25))
CountDist      = "Poisson"
# MargParm       = 10
beta0          = 2
beta1          = 0.4
MargParm       = c(beta0, beta1)
ARParm         = c(0.75)
MAParm         = NULL
EstMethod      = "PFR"
ARMAorder      = c(length(ARParm),length(MAParm))
ParticleNumber = 10
epsilon        = 0.5
# initialParam   = c(MargParm*0.5, 0*ARParm, 0.5*MAParm)
initialParam   = NULL
TrueParam      = c(MargParm, ARParm, MAParm)
theta          = initialParam
OptMethod      = "bobyqa"
Task           = 'Optimization'
OutputType     = "data.frame"
ParamScheme    = 1

# simulate data
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm,MAParm, Regressor)

# Fit the model using the lgc wrapper
mylgc = lgc(DependentVar   = DependentVar,
            Regressor      = Regressor,
            EstMethod      = EstMethod,
            CountDist      = CountDist,
            ARMAorder      = ARMAorder,
            ParticleNumber = ParticleNumber,
            epsilon        = epsilon,
            initialParam   = initialParam,
            TrueParam      = TrueParam,
            Task           = Task,
            OptMethod      = OptMethod,
            OutputType     = OutputType,
            ParamScheme    = ParamScheme
)



