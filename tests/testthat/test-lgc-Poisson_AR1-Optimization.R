# load libraries
# library(countsFun)
# library(tictoc)
# library(optimx)
# library(ltsa)
# library(itsmr)
# library(numDeriv)
# library(MASS)

test_that("LGC wrapper works for ARMA(1,1) with bad initial param.", {

# Specify model and methods
n              = 50
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

# set initial parameters
initialParam   = c(10,0.1,0.2)

# save the data in a data frame
df = data.frame(DependentVar)

# specify the regression model
formula = DependentVar~0


# call the wrapper
a = lgc(formula,df, EstMethod, CountDist, ARMAModel, ParticleNumber,
         epsilon, initialParam, TrueParam,  Task, SampleSize, nsim, no_cores, OptMethod,
         OutputType, ParamScheme, maxdiff)

# for the set.seed(2) I will get the following likelihood: 302.57734
expect_equal(a$ParamEstimates[1], 3.057192, tolerance = 10^(-4))
expect_equal(a$ParamEstimates[2], 0.4400711, tolerance = 10^(-4))
expect_equal(a$ParamEstimates[3], 0.07076272, tolerance = 10^(-4))

nobs = nobs(a)
expect_equal(nobs, n)

})
