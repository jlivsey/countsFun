# load libraries
# library(countsFun)
# library(tictoc)
# library(optimx)
# library(ltsa)
# library(itsmr)
# library(numDeriv)
# library(MASS)

test_that("LGC wrapper works for ARMA(1,1)", {

# Specify model and methods
n              = 100
# Regressor      = cbind(rep(1,n),rbinom(n,1,0.25))
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 3
ARParm         = 0.5
MAParm         = 0.2
#ARParm         = NULL
#MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))
ParticleNumber = 5
epsilon        = 0.5
EstMethod      = "PFR"
TrueParam      = c(MargParm,ARParm,MAParm)
Task           = 'Optimization'
SampleSize     = NULL
nsim           = NULL
no_cores       = NULL
OptMethod      = "bobyqa"
OptMethod      = "L-BFGS-B"

OutputType     = "list"
ParamScheme    = NULL
maxdiff        = 10^(-6)
# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

initialParam   = c(10,0.1,0.2)
# populate the model scheme
mod = ModelScheme(DependentVar, Regressor, EstMethod, ARMAModel, CountDist,ParticleNumber, epsilon,
                  initialParam, TrueParam, Task,SampleSize, OptMethod, OutputType, ParamScheme, maxdiff)

# call the wrapper
a = lgc(DependentVar, Regressor, EstMethod, CountDist, ARMAModel, ParticleNumber,
         epsilon, initialParam, TrueParam,  Task, SampleSize, nsim, no_cores, OptMethod,
         OutputType, ParamScheme, maxdiff)
a$ParamEstimates

# for the set.seed(2) I will get the following likelihood: 302.57734
expect_equal(a$ParamEstimates[1], 2.8172329)
expect_equal(a$ParamEstimates[2], 0.43318923)
expect_equal(a$ParamEstimates[3], 0.17216145)

})
