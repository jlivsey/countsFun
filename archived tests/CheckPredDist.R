#================================================================================#
# PURPOSE: Check if PIT plot works
#
# Notes:  Focus on AR(2) Poisson
#
#
# DATE: July 2020
# Author: Stefanos Kechagias
#================================================================================#


rm(list=ls())
library(countsFun)
library(optimx)
library(FitAR)

# list with true ARMA parameters
n = 800
p=2
q=0
ARMAorder = c(p,q)
# list with true ARMA parameters
ARMAmodel = list(NULL,NULL)
ARMAmodel[[1]] = c(0.5, 0.2)

set.seed(120)
data = sim_poisson(n, ARMAmodel, 4)

#other parameters
OptMethod = "bobyqa"
CountDist = "Poisson"
epsilon   = 0.5
MaxCdf    = NULL
nHC       = NULL
Model     = NULL

ParticleNumber = 400
Regressor = NULL
LB = c(0.001, -0.99, -0.99)
UB = c(100, 0.99, 0.99)

theta = c(4,0.5, 0.2)


# fit the model
mod = FitMultiplePF_Res(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon, LB, UB, OptMethod)

# compute rped Dist
#mod =c (3.765643, 0.450734, 0.1917529 )
predDist = PDvaluesARp(mod[1:3], data, Regressor, ARMAorder, ParticleNumber,CountDist, epsilon)




PIT = PITvalues(10, predDist)

barplot(PIT)


