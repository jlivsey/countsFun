
#============================================================================================#
# PURPOSE: Compare the synthesis of GenPois with regressor. I am looking at two things:
#          1. compare sim_lgc_old vs sim_lgc (a more compact version of sim_lgc)
#          2. compare the results I using the distribution functions of VGAM
#             It seems VGAM functions are less likely to produce NA
#
# Author: Stefanos Kechagias
# Date:   August 2024
#============================================================================================#


# the following ios an example of how this would affect synthesis
SampleSize     = 100
n              = SampleSize
b0             = 0.5
b1             = 2
beta           = c(b0,b1)
alpha          = 0.5
MargParm       = c(b0,b1,alpha)
set.seed(3)
e              = rbinom(n,1,0.1)
Regressor      = cbind(rep(1,n),e)
ARParm         = c(0.5, 0.2)
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))
CountDist      = "Generalized Poisson 2"
TrueParam      = c(MargParm,ARParm, MAParm)
Task           = "Synthesis"


# Pars the model
mod = ModelScheme(CountDist      = CountDist,
                  ARMAModel      = ARMAModel,
                  TrueParam      = TrueParam,
                  Regressor      = Regressor,
                  Task           = Task,
                  SampleSize     = SampleSize)

# Simulate with the consice new version
set.seed(1)
x2 = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm, Regressor)
x2

# simulate with the older version
set.seed(1)
x1   = sim_lgc_old(SampleSize, CountDist, MargParm, ARParm, MAParm, Regressor)
x1

# Simulate with our own implementation of distribution functions
set.seed(1)
CountDist      = "Generalized Poisson"
x3   = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm, Regressor)
x3
