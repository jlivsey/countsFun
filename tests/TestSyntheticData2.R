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
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 10
ARParm         = c(0.75)
MAParm         = NULL
EstMethod      = "PFR"
ARMAorder      = c(length(ARParm),length(MAParm))
ParticleNumber = 10
epsilon        = 0.5
initialParam   = c(MargParm*0.5, 0*ARParm, 0.5*MAParm)
initialParam   = NULL
theta          = initialParam
OptMethod      = NULL
Optimization   = 0
OutputType     = "matrix"

# simulate data
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm,MAParm)

# Fit the model using the lgc wrapper
mylgc = lgc(DependentVar   = DependentVar,
            Regressor      = NULL,
            EstMethod      = EstMethod,
            CountDist      = CountDist,
            ARMAorder      = ARMAorder,
            ParticleNumber = ParticleNumber,
            epsilon        = epsilon,
            initialParam   = initialParam,
            Optimization   = Optimization,
            OptMethod      = OptMethod,
            OutputType     = OutputType
)




