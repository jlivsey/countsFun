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
ARParm         = c(0.5,0.2)
# MAParm         = 0.5
#ARParm         = NULL
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))
ParticleNumber = 1
epsilon        = 0.5
EstMethod      = "PFR"
TrueParam      = c(MargParm,ARParm,MAParm)
initialParam   = TrueParam
Task           = 'Optimization'
SampleSize     = NULL
nsim           = NULL
no_cores       = NULL
OptMethod      = "bobyqa"
OutputType     = "data.frame"
ParamScheme    = NULL
maxdiff        = 10^(-8)
# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

mod = ModelScheme(DependentVar, Regressor, EstMethod, ARMAModel, CountDist,ParticleNumber, epsilon,
                  initialParam, TrueParam, Task,SampleSize, OptMethod, OutputType, ParamScheme, maxdiff)


if (is.null(initialParam)){
  theta  = InitEst = InitialEstimates(mod)
  mod$initialParam = InitEst
}else{
  theta  = InitEst = mod$initialParam
}



t1 = tic()
for (i in 1:1){
  a1 = ParticleFilter_Res_ARMA(theta,mod)
}
t1 = tic() - t1

t2 = tic()
for (i in 1:1){
  a2 = ParticleFilter_Res_AR(theta,mod)
}
t2 = tic() - t2

t1
t2
a1
a2
















