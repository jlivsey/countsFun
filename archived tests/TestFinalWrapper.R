#====================================================================================#
# PURPOSE      Test of calling the main wrapper function of the timecounts package.
#
#
# Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras, Jiajie Kong
# Date         Fall 2022
# Version      4.2.1
#====================================================================================#

# load libraries
library(countsFun)
library(tictoc)
library(optimx)
library(ltsa)
library(itsmr)
library(lavaSearch2)
library(numDeriv)
library(sandwich)
library(MASS)

# FIX ME: Check Where is this used?
symmetrize     = lavaSearch2:::symmetrize

# load the data
mysales        = read.csv("data/MySelectedSeries.csv")

# specify parameters
n              = 104
epsilon        = 0.5
MaxCdf         = 1000
nHC            = 30
ParticleNumber = 10
data           = mysales$MOVE[1:n]
ARMAorder      = c(3,0)
Regressor      = cbind(rep(1,length(mysales$Buy[1:n] )),mysales$Buy[1:n] )
CountDist      = "Negative Binomial"
# the Poisson will yield an infinite value.
# CountDist      = "Poisson"
EstMethod      = "PFR"
OptMethod      = "bobyqa"
maxit          = 0
# initialParam   = c(2.597666, 1.15373, 1.197833, -0.3748205, 0.226, 0.227)
initialParam   = NULL

#initialParam   = c(2.13258844,  1.16177357, -0.39141331,  0.08239859,  0.05745040)
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon, initialParam, EstMethod, maxit)

# stop if there was an error in model specification
# if(mod$error) stop(mod$errorMsg)
#
# # fix me: I need a function that computes initial parameters
# if (is.null(initialParam)){
#   theta  = InitialEstimates(mod)
# }else{
#   theta  = mod$initialParam
# }
# theta = c(2.264, 1.01,-0.341, 0.223, 0.291)
# ParticleFilter_Res_AR(theta, mod)

# call the wrapper
a = countC(data, Regressor, CountDist, EstMethod, ARMAorder,
                   nHC, MaxCdf, ParticleNumber, epsilon, initialParam ,
                   OptMethod, maxit)







#
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon, initialParam, EstMethod)
#
#
# theta = c(2.264, 1.01, 1.21, -0.341, 0.223, 0.291)
# optim.output <- optimx(par            = theta,
#                        fn             = GaussianLogLik,
#                        data           = data,
#                        Regressor      = Regressor,
#                        mod            = mod,
#                        lower          = mod$LB,
#                        upper          = mod$UB,
#                        method         = OptMethod,
#                        hessian        = TRUE)
#
#
# ParmEst  = as.numeric(optim.output[1:mod$nparms])
# loglik   = optim.output$value
#
#
# # optim.output <- optimx(par            = theta,
# #                        fn             = GaussianLogLik,
# #                        data           = data,
# #                        Regressor      = Regressor,
# #                        mod            = mod,
# #                        lower          = mod$LB,
# #                        upper          = mod$UB,
# #                        control=list(all.methods=TRUE),
# #                        hessian        = TRUE)
#
# # optim.output <- optimx(par            = theta,
# #                        fn             = GaussianLogLik,
# #                        data           = data,
# #                        Regressor      = Regressor,
# #                        mod            = mod,
# #                        lower          = mod$LB,
# #                        upper          = mod$UB,
# #                        method         = "Rvmmin",
# #                        hessian        = TRUE)
#
#
# # save estimates, loglik value and diagonal hessian
# ParmEst  = as.numeric(optim.output[5,1:mod$nparms])
# loglik   = optim.output$value
# convcode = optim.output$convcode
# kkt1     = optim.output$kkt1
# kkt2     = optim.output$kkt2
#
#
# # compute sandwich standard errors
# se = sand(ParmEst, data, Regressor, mod)
#
#
#
# h <- gHgen(fn        = GaussianLogLik,
#            par       = theta,
#            data      = data,           # additional arg for GaussLogLik
#            Regressor = Regressor,      # additional arg for GaussLogLik
#            mod       = mod)            # additional arg for GaussLogLik
# SE.hess <- sqrt(diag(solve(h$Hn)))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# #sand(theta, data, Regressor, mod)
#
#
# countC(data, Regressor, CountDist, EstMethod, ARMAorder,nHC, MaxCdf,ParticleNumber, epsilon, initialParam, OptMethod, maxit)
#
#
#
#
#
#
#
#
#
# #Regressor = cbind(rep(1,length(Buy)),Buy)
# theta     = c(2, 0.5)
# initialParam = theta
# mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon, initialParam)
#
# sand(theta, data, Regressor, mod)
#
# #FitGaussianLogLik (theta, data, Regressor, mod, OptMethod)
#
#
#







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
