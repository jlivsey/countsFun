# load necessary libraries
library(optimx)
library(ltsa)
library(countsFun)
library(itsmr)
library(tictoc)

# Specify model and methods
n              = 1000
# Regressor      = cbind(rep(1,n),rbinom(n,1,0.25))
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 3
ARParm         = c(0.8, -0.25)
#MAParm         = c(0.2, 0.5,0.2, 0.1)
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
OptMethod      = "L-BFGS-B"
OutputType     = "data.frame"
ParamScheme    = NULL
maxdiff        = 10^(-6)
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


t3 = tic()
for (i in 1:1){
  a3 = ParticleFilter_Res_ARMA_old(theta,mod)
}
t3 = tic() - t3

t4 = tic()
for (i in 1:1){
  a4 = ParticleFilter_Res_AR(theta,mod)
}
t4 = tic() - t4

t5 = tic()
for (i in 1:1){
  a5 = ParticleFilter_Res_AR_old(theta,mod)
}
t5 = tic() - t5

t1
a1
a3
t3
a4
t4
a5
t5

