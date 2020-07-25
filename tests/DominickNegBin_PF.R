#======================================================================================================#
#Purpose:   Fit several NegBin PF models for Dominick data
#
# Author:   Stefanos Kechagias
# Team:     Vladas Pipiras, James Livsey, Stefanos Kechagias, Robert Lund, Yisu Jia
# Date:     July 2020
#=====================================================================================================#

# load libraries
library(countsFun)
library(tictoc)
library(optimx)
library(FitAR)


# load the data
mysales = read.csv("/Users/stef/Desktop/countsFun/data/MySelectedSeries.csv")

# attach the datafrmae
n = 104
Smallsales  = mysales[1:n,]
attach(Smallsales)

# regressor variable with intercept
Regressor = cbind(rep(1,length(Buy)),Buy)

#other parameters
OptMethod = "bobyqa"
CountDist = "Negative Binomial"
epsilon   = 0.5
Particles = 200
MaxCdf    = NULL
nHC       = NULL
Model     = NULL
#===================================   AR(1)  ========================================================================================#
ARMAorder = c(1,0)
theta = c(2, 1, 0.5, -0.3)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf)
UB = c(100, 100, Inf, Inf)
tic()
mod1 = FitMultiplePFResReg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
#====================================================================================================================================#


#===================================   AR(2)  ========================================================================================#
ARMAorder = c(2,0)
theta = c(2, 1, 0.5, -0.3, 0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf)
tic()
mod2 = FitMultiplePFResReg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
#====================================================================================================================================#



#===================================   AR(3)  ========================================================================================#
ARMAorder = c(3,0)
theta = c(2, 1, 0.5, -0.3, 0, 0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf)
tic()
mod3 = FitMultiplePFResReg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
#====================================================================================================================================#


#===================================   MA(1)  ========================================================================================#
ARMAorder = c(0,1)
theta = c(2, 1, 0.5, -0.3)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf)
UB = c(100, 100, Inf, Inf)
tic()
mod4 = FitMultiplePFResReg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
#====================================================================================================================================#


#===================================   MA(2)  ========================================================================================#
ARMAorder = c(0,2)
theta = c(2, 1, 0.5, -0.3,0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf)
tic()
mod5 = FitMultiplePFResReg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
#====================================================================================================================================#

#===================================   MA(3)  ========================================================================================#
ARMAorder = c(0,3)
theta = c(2, 1, 0.5, -0.3,0,0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf)
tic()
mod6 = FitMultiplePFResReg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
#====================================================================================================================================#


