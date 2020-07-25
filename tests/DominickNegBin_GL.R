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
MaxCdf    = 5000
nHC       = 30
Model     = 0
OptMethod = "bobyqa"
CountDist = "Negative Binomial"
epsilon   = NULL
Particles = NULL
#===================================   AR(1)  ========================================================================================#
ARMAorder = c(1,0)
theta = c(2, 1, 0.5, -0.3)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf)
UB = c(100, 100, Inf, Inf)
tic()
mod1 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
#2.435945 0.5343325 0.9538895 -0.5525136 0.08895981 0.1083391 0.1999505 0.1118167 313.995 827.1292 837.7067 827.5332      0    1    1
#====================================================================================================================================#

GaussLogLikNB_Reg(theta, MOVE, Regressor, ARMAorder, MaxCdf, nHC, CountDist)


#===================================   AR(2)  ========================================================================================#
ARMAorder = c(2,0)
theta = c(2, 1, 0.5, -0.3, 0.1)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf)
tic()
mod2 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
# 2.42525 0.5304039 1.051025 -0.8428225 -0.2595247 0.08729126 0.1058773 0.1956355 0.2155361 0.2066433 313.3422 827.8236 841.0456 828.4359      0    1    1
#====================================================================================================================================#


#===================================   AR(3)  ========================================================================================#
ARMAorder = c(3,0)
theta = c(2, 1, 0.5, -0.3, 0.1, 0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf)
tic()
mod3 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
# 2.43246 0.5691376 0.884361 -0.4471752 0.1453361 0.2078561 0.09579149 0.115404 0.2071525 0.1748645 0.1710018 0.130477 312.9453 829.0298 844.8961 829.8957      0    1    1
#====================================================================================================================================#



#===================================   MA(1)  ========================================================================================#
ARMAorder = c(0,1)
theta = c(2, 1, 0.5, -0.3)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf)
UB = c(100, 100, Inf, Inf)
tic()
mod4 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
# 2.436593 0.5196907 0.9183684 -0.4009871 0.08437518 0.1125265 0.1954814 0.1452232 317.3437 833.8266 844.4041 834.2306      0    1    1
#====================================================================================================================================#




#===================================   MA(2)  ========================================================================================#
ARMAorder = c(0,2)
theta = c(2, 1, 0.5, -0.3, 0.1)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf)
tic()
mod5 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
# 2.405191 0.5852359 0.954474 -0.5594147 0.4405829 0.09592033 0.1053513 0.1991286 7.007475e-06 5.518935e-06 311.8423 824.8237 838.0457 825.436      0    0    0
#====================================================================================================================================#


#===================================   MA(3)  ========================================================================================#
ARMAorder = c(0,3)
theta = c(2, 1, 0.5, -0.3, 0.1, 0.3)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf)
tic()
mod6 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
# 2.402127 0.6474095 0.9226765 -0.6471064 0.6855718 0.3326952 0.01898355 0.0969907 0.1416 0.01401753 0.01485076 0.007206799 310.7309 824.6011 840.4674 825.467      0    0    0
#====================================================================================================================================#


#===================================   MA(4)  ========================================================================================#
ARMAorder = c(0,4)
theta = c(2, 0.5, 1, -0.3, 0.4, 0.3,0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf, Inf)
tic()
mod7 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
# 2.413522 0.6342288 0.9022111 -0.5932237 0.5940733 0.2381122 0.05081465 0.1004769 0.1016517 0.18668 2.171341e-09 1.959557e-09 8.791535e-10 1.719617e-10 310.6195 826.3781 844.8889 827.5448      0    0    0
#====================================================================================================================================#



#===================================   White Noise  ========================================================================================#
ARMAorder = c(0,0)
theta = c(2, 1, 0.5)
# lower and upper bounds
LB = c(-100, -100, 0.001)
UB = c(100, 100, Inf)
tic()
mod0 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
# 2.397635 0.6174012 0.9747719 0.1138878 0.1233839 0.255244 323.4319 844.003 851.9362 844.243      0    1    1
#====================================================================================================================================#



#===================================   ARMA(1,1)  ========================================================================================#
ARMAorder = c(1,1)
theta = c(2, 1, 0.5, -0.3, 0)
# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf)
tic()
mod8 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)
toc()
# 2.421937 0.5335326 1.115948 -0.5003381 -0.9999782 0.08809474 0.1069163 0.1819092 0.09258593 3.425106 313.5258 828.1908 841.4127 828.803      0    1    1
#====================================================================================================================================#


NegBinResults = list(mod0,mod1,mod2,mod3,mod4,mod5,mod6,mod7)


Models = c("WN", "AR(1)","AR(2)","AR(3)","MA(1)","MA(2)","MA(3)", "MA(4)")
AIC = rbind(mod0[,"AIC"],mod1[,"AIC"],mod2[,"AIC"],mod3[,"AIC"],mod4[,"AIC"],mod5[,"AIC"],mod6[,"AIC"],mod7[,"AIC"])
BIC = rbind(mod0[,"BIC"],mod1[,"BIC"],mod2[,"BIC"],mod3[,"BIC"],mod4[,"BIC"],mod5[,"BIC"],mod6[,"BIC"],mod7[,"BIC"])
AICc = rbind(mod0[,"AICc"],mod1[,"AICc"],mod2[,"AICc"],mod3[,"AICc"],mod4[,"AICc"],mod5[,"AICc"],mod6[,"AICc"],mod7[,"AICc"])
All = data.frame(Models, AIC, BIC, AICc)
save.image(file = "C:/Users/Stef/Desktop/countsFun/data/NegBin_GL.RData")

# optim.output4 <- optimx(par       = theta,
#                         fn        = GaussLogLikNB_Reg,
#                         data      = MOVE,
#                         Regressor = Regressor,
#                         ARMAorder = ARMAorder,
#                         MaxCdf    = 5000,
#                         nHC       = nHC,
#                         hessian   = TRUE,
#                         lower     = LB,
#                         upper     = UB,
#                          control = list(all.methods = TRUE,trace = 2)
# )










