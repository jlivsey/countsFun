#======================================================================================================#
#Purpose:   Fit several NegBin GL models for Dominick data
#
# Author:   Stefanos Kechagias
# Team:     Vladas Pipiras, James Livsey, Stefanos Kechagias, Robert Lund, Yisu Jia
# Date:     July 2020
#=====================================================================================================#

# load libraries
library(countsFun)
library(tictoc)
library(lavaSearch2)
library(optimx)
library(FitAR)
symmetrize <- lavaSearch2:::symmetrize

# load the data
mysales = read.csv("/Users/stef/Desktop/countsFun/data/MySelectedSeries.csv")

# attach the datafrmae
n = 104
Smallsales  = mysales[1:n,]
attach(Smallsales)

# regressor variable with intercept
Regressor = cbind(rep(1,length(Buy)),Buy)

#other parameters
MaxCdf = 5000
nHC = 30
Model = 0
OptMethod = "nlminb"
CountDist = "Generalized Poisson"
#===================================   AR(1)  ========================================================================================#
ARMAorder = c(1,0)
theta = c(2, 1, 0.5, -0.3)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf)
UB = c(100, 100, Inf, Inf)
tic()
mod1 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()


#
#====================================================================================================================================#



#===================================   AR(2)  ========================================================================================#
ARMAorder = c(2,0)
theta = c(2, 1, 0.5, -0.3, 0.1)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf)
tic()
mod2 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()
#
#====================================================================================================================================#


#===================================   AR(3)  ========================================================================================#
ARMAorder = c(3,0)
theta = c(2, 1, 0.5, -0.3, 0.1, 0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf)
tic()
mod3 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()
#
#====================================================================================================================================#



#===================================   MA(1)  ========================================================================================#
ARMAorder = c(0,1)
theta = c(2, 1, 0.5, -0.3)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf)
UB = c(100, 100, Inf, Inf)
GaussLogLikGP_Reg(theta, MOVE, Regressor, ARMAorder, MaxCdf, nHC, CountDist)
tic()
mod4 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()
#
#====================================================================================================================================#




#===================================   MA(2)  ========================================================================================#
ARMAorder = c(0,2)
theta = c(2, 1, 0.5, -0.3, 0.1)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf)
tic()
mod5 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()
#
#====================================================================================================================================#


#===================================   MA(3)  ========================================================================================#
ARMAorder = c(0,3)
theta = c(2, 1, 0.5, -0.3, 0.1, 0.3)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf)
tic()
mod6 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()
#
#====================================================================================================================================#


#===================================   MA(4)  ========================================================================================#
ARMAorder = c(0,4)
theta = c(2, 0.5, 1, -0.3, 0.4, 0.3,0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf, Inf)
tic()
mod7 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()
#
#====================================================================================================================================#



#===================================   White Noise  ========================================================================================#
ARMAorder = c(0,0)
theta = c(2, 1, 0.5)
# lower and upper bounds
LB = c(-100, -100, 0.001)
UB = c(100, 100, Inf)
tic()
mod0 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()
#
#====================================================================================================================================#



#===================================   ARMA(1,1)  ========================================================================================#
ARMAorder = c(1,1)
theta = c(2, 1, 0.5, -0.3, 0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf)
tic()
mod8 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist)
toc()
#
#====================================================================================================================================#




Models = c("WN", "AR(1)","AR(2)","AR(3)","MA(1)","MA(2)","MA(3)", "MA(4)")
AIC = rbind(mod0[,"AIC"],mod1[,"AIC"],mod2[,"AIC"],mod3[,"AIC"],mod4[,"AIC"],mod5[,"AIC"],mod6[,"AIC"],mod7[,"AIC"])
BIC = rbind(mod0[,"BIC"],mod1[,"BIC"],mod2[,"BIC"],mod3[,"BIC"],mod4[,"BIC"],mod5[,"BIC"],mod6[,"BIC"],mod7[,"BIC"])
AICc = rbind(mod0[,"AICc"],mod1[,"AICc"],mod2[,"AICc"],mod3[,"AICc"],mod4[,"AICc"],mod5[,"AICc"],mod6[,"AICc"],mod7[,"AICc"])
All = data.frame(Models, AIC, BIC, AICc)

save.image(file = "C:/Users/Stef/Desktop/countsFun/tests/GenPois_GL.RData")


# optim.output4 <- optimx(par       = theta,
#                         fn        = GaussLogLikNB_Reg,
#                         data      = MOVE,
#                         Regressor = Regressor,
#                         CountDist = CountDist,
#                         ARMAorder = ARMAorder,
#                         MaxCdf    = 5000,
#                         nHC       = nHC,
#                         hessian   = TRUE,
#                         lower     = LB,
#                         upper     = UB,
#                          control = list(all.methods = TRUE,trace = 2)
# )










