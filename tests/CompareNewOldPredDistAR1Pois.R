#================================================================================#
# PURPOSE: Compare new and old code for computing predictive distributions
#
# Notes:  Focus on AR(1) Poisson
#
#
#
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
n = 200
p=1
q=0
ARMAorder = c(p,q)
# list with true ARMA parameters
ARMAmodel = list(NULL,NULL)
ARMAmodel[[1]] = 0.5

set.seed(3)
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
LB = c(0.001, -0.99)
UB = c(100, 0.99)

theta = c(4,0.5)


mod = FitMultiplePF_Res(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon, LB, UB, OptMethod)


#set.seed(1)
predDist1 = PDvaluesAR1(mod[1:2], as.numeric(data), Regressor, ARMAorder, ParticleNumber,CountDist)
predDist2 = PDvaluesPoissonAR1Old(mod[1], mod[2], as.numeric(data),   ParticleNumber)
sum(predDist1 - predDist2)






#PIT = PITvalues(10, predDist)

#barplot(PIT)


