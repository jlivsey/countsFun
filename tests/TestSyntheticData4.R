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
library(tictoc)

# Specify model and methods
n              = 200
# Regressor      = cbind(rep(1,n),rbinom(n,1,0.25))
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 3
#ARParm         = c(0.5,0.25, -0.4)
#MAParm         = NULL
ARParm         = NULL
MAParm         = c(0.5,0.25, -0.4)
ARMAModel      = c(length(ARParm),length(MAParm))
ParticleNumber = 10
epsilon        = 0.5
EstMethod      = "PFR"
initialParam   = NULL
initialParam   = c(MargParm,ARParm,MAParm)
TrueParam      = NULL
Task           = 'Optimization'
SampleSize     = NULL
nsim           = NULL
no_cores       = NULL
OptMethod      = "bobyqa"
OutputType     = "data.frame"
ParamScheme    = NULL


# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)


# parse all the parameters and the data into a list called mod
mod = ModelScheme(DependentVar, Regressor, EstMethod, ARMAModel, CountDist,ParticleNumber, epsilon,
                  initialParam, TrueParam, Task,SampleSize, OptMethod, OutputType, ParamScheme)

if (is.null(initialParam)){
  theta  = InitEst = InitialEstimates(mod)
  mod$initialParam = InitEst
}else{
  theta  = InitEst = mod$initialParam
}

t0 = tic()
for (i in 1:300){
  a1 = ParticleFilter_Res_AR(theta,mod)
}
t0 = tic() - t0

t1 = tic()
for (i in 1:300){
a2 = ParticleFilter_Res_MA(theta,mod)
}
t1 = tic() - t1

t2 = tic()
for (i in 1:300){
a3 = ParticleFilter_Res_ARMA(theta,mod)
}
t2 = tic() - t2
t0
t1
t2





















