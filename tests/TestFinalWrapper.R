# load libraries
library(countsFun)
library(tictoc)
library(optimx)
library(FitAR)
library(itsmr)
library(lavaSearch2)
library(numDeriv)
library(sandwich)

symmetrize <- lavaSearch2:::symmetrize

# load the data
mysales = read.csv("data/MySelectedSeries.csv")

# attach the datafrmae
n = 104
Smallsales  = mysales[1:n,]
MOVE = Smallsales$MOVE
Buy = Smallsales$Buy


#other parameters
epsilon   = 0.5
MaxCdf    = 1000
nHC       = 30
ParticleNumber = 10
epsilon = 0.5
data = MOVE


# test sand
ARMAorder = c(1,0)
Regressor = NULL
CountDist = "Poisson"
EstMethod="GL"
initialParam = NULL
OptMethod = "bobyqa"
#Regressor = cbind(rep(1,length(Buy)),Buy)
theta     = c(2, 0.5)
initialParam = theta
mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon, initialParam)

sand(theta, data, Regressor, mod)





#countC(data, Regressor, CountDist, EstMethod, ARMAorder,nHC, MaxCdf,ParticleNumber, epsilon, initialParam, OptMethod)









#Regressor = cbind(rep(1,length(Buy)),Buy)
theta     = c(2, 0.5)
initialParam = theta
mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon, initialParam)

sand(theta, data, Regressor, mod)

#FitGaussianLogLik (theta, data, Regressor, mod, OptMethod)










#
#
#
#
#
#
#
#
# GaussianLogLik(theta, data, Regressor, mod)
# GaussLogLik = function(theta, data)
#
#
#
#
# optim.output
# # xt = data
# # nparms  = length(theta)
# # n = length(xt)
# # # allocate memory to save parameter estimates, hessian values, and loglik values
# # ParmEst = matrix(0,nrow=1,ncol=nparms)
# # se =  matrix(NA,nrow=1,ncol=nparms)
# # loglik = rep(0,1)
# # convcode = rep(0,1)
# # kkt1 = rep(0,1)
# # kkt2 = rep(0,1)
# #
# #
# # optim.output <- optimx(par            = theta,
# #                        fn             = GaussianLogLik,
# #                        data           = xt,
# #                        Regressor      = Regressor,
# #                        mod            = mod,
# #                        lower          = mod$LB,
# #                        upper          = mod$UB,
# #                        method         = OptMethod,
# #                        hessian        = TRUE)
# #
# #
#
# out = FitGaussianLogLik(theta, data, Regressor, mod, "L-BFGS-B")
#
# s = sand(theta, data, Regressor, mod)
#
#
#
#
#
# #====================================== PF =============================================#
#
#
# # test PF without regressor + AR
# ARMAorder = c(2,0)
# Regressor = NULL
# CountDist = "Generalized Poisson"
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
# theta     = c(2,0.5, 0.5, -0.3)
# ParticleFilter_Res(theta, data, Regressor, mod)
# LB = c(0.01, 0.01, -Inf, -Inf)
# UB = c(Inf, Inf, Inf, Inf)
# FitMultiplePF_Res(theta, data, Regressor, mod, OptMethod)
#
#
#
#
# # test PF with regressor + MA
# ARMAorder = c(0,1)
# Regressor = cbind(rep(1,length(Buy)),Buy)
# CountDist = "Negative Binomial"
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
# theta     = c(2, 1, 0.5, -0.3)
# ParticleFilter_Res(theta, data, Regressor, mod)
#
#
#
# #====================================== GL =============================================#
# # test GL with regressor + AR
# ARMAorder = c(2,0)
# Regressor = cbind(rep(1,length(Buy)),Buy)
# CountDist = "Generalized Poisson"
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC, ParticleNumber, epsilon )
# theta     = c(2,0.5, 0.5,0.5, -0.3)
# GaussianLogLik(theta, data, Regressor, mod)
# OptMethod = "bobyqa"
# LB = c(-100, -100, 0.001, -Inf, -Inf)
# UB = c(100, 100, Inf, Inf, Inf)
# GaussLogLikGP_Reg(theta, data, Regressor, ARMAorder, MaxCdf, nHC, CountDist)
# # FitGaussianLogLik(theta, data, Regressor, mod, OptMethod)
#
#
# # test GL without regressor + AR
# ARMAorder = c(2,0)
# Regressor = NULL
# CountDist = "Generalized Poisson"
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC, ParticleNumber, epsilon )
# theta     = c(2, 0.5,0.5, -0.3)
# GaussianLogLik(theta, data, Regressor, mod)
# OptMethod = "bobyqa"
# GaussLogLikGP(theta, data, ARMAorder, MaxCdf, nHC)
#
#
#
# # test GL with regressor + MA
# ARMAorder = c(0,1)
# Regressor = cbind(rep(1,length(Buy)),Buy)
# CountDist = "Negative Binomial"
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
# theta     = c(2, 1, 0.5, -0.3)
# GaussianLogLik(theta, data, Regressor, mod)
#
#
#
# # test Gen Pois with regressor and WN
# ARMAorder = c(0,0)
# Regressor = cbind(rep(1,length(Buy)),Buy)
# CountDist = "Poisson"
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
# theta     = c(2,1)
# GaussianLogLik(theta, data, Regressor, mod)
#
