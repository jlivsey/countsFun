
library(countsFun)
library(optimx)
library(ltsa)
library(tictoc)

n              = 1000
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 2
ARParm         = c(0.75)
EstMethod      = "PFR"
ARMAorder      = c(1,0)
ParticleNumber = 10
epsilon        = 0.5
initialParam   = c(3, 0.1)
OptMethod      = "bobyqa"

theta = initialParam
data           = as.numeric(read.csv("C:/Users/statha/OneDrive - SAS/Desktop/countsFun/data/SimPois2_AR1_075_N1k.csv")[,1])
mod            = ModelScheme(data, Regressor, EstMethod, ARMAorder, CountDist,ParticleNumber, epsilon, initialParam)


# tic()
# test(initialParam, mod)
# toc()
#
# tic()
# ParticleFilter_Res_AR(initialParam, mod)
# toc()


optim.output1 <- optimx(
  par     = theta,
  fn      = ParticleFilter_Res_AR,
  lower   = mod$LB,
  upper   = mod$UB,
  hessian = TRUE,
  method  = mod$OptMethod,
  mod     = mod,
  control = list(trace = 6))

optim.output2 <- optimx(
  par     = theta,
  fn      = test,
  lower   = mod$LB,
  upper   = mod$UB,
  hessian = TRUE,
  method  = mod$OptMethod,
  mod     = mod)
optim.output1
optim.output2

