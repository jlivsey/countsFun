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
MOVE = Smallsales$MOVE
Buy = Smallsales$Buy


# regressor variable with intercept
Regressor = cbind(rep(1,length(Buy)),Buy)

#other parameters
OptMethod = "bobyqa"
CountDist = "Generalized Poisson"
epsilon   = 0.5
Particles = 400
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
mod4 = FitMultiplePFResReg(theta, MOVE, Regressor, CountDist, 100, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
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


#===================================   WN    ========================================================================================#
ARMAorder = c(0,0)
theta = c(2, 1, 0.5)
# lower and upper bounds
LB = c(-100, -100, 0.001)
UB = c(100, 100, Inf)
data = MOVE
tic()
mod0 = FitMultiplePFResReg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
#====================================================================================================================================#

GenPoisResultsPF = list(mod0,mod1,mod2,mod3,mod4,mod5,mod6)

Models = c("WN","AR(1)","AR(2)","AR(3)","MA(1)","MA(2)","MA(3)")
AIC = rbind(mod0[,"AIC"], mod1[,"AIC"],mod2[,"AIC"],mod3[,"AIC"],mod4[,"AIC"],mod5[,"AIC"],mod6[,"AIC"])
BIC = rbind(mod0[,"BIC"], mod1[,"BIC"],mod2[,"BIC"],mod3[,"BIC"],mod4[,"BIC"],mod5[,"BIC"],mod6[,"BIC"])
AICc = rbind(mod0[,"AICc"], mod1[,"AICc"],mod2[,"AICc"],mod3[,"AICc"],mod4[,"AICc"],mod5[,"AICc"],mod6[,"AICc"])
All = data.frame(Models, AIC, BIC, AICc)

save.image(file = "C:/Users/Stef/Desktop/countsFun/tests/GenPois_PF.RData")


# Neg Binomial exact GLM
# summary(m1 <- glm.nb(MOVE ~ Buy))
