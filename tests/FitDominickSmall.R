#======================================================================================================#
#Purpose:   Load a series from the Dominick retail data that satisfies
#           1) Has negative lag 1 acf
#           2) Has small mean (<30)
#           3) Is cross correlated with a Sale dummy variable
#           4) After a regression with the dummy variable the residuals have significant lag1 acf.
#
# Notes:    The following steps briefly describe the procedure we followed
#
#           Step 1: Filter intermittent and retired series
#           Step 2: Focus on series with significant lag 1 negative acf
#           Step 3: Select series with mean < 30
#           Step 4: Run a regression on a Buy one get on free dummy variable and examine residuals
#           Step 5: Trial and error to select a series whose residuals have significant lag 1 neg acf
#
# Author:   Stefanos Kechagias
# Team:     Vladas Pipiras, James Livsey, Stefanos Kechagias, Robert Lund, Yisu Jia
# Date:     March 2020
#=====================================================================================================#

# load libraries
library(gcmr)
library(countsFun)
#library(ggplot2)
library(tictoc)
library(lavaSearch2)
library(optimx)
library(VGAM)

library(FitAR)
symmetrize <- lavaSearch2:::symmetrize


# load the data
#mysales = read.csv("C:/Users/Stef/Desktop/countsFun/data/MySelectedSeries.csv")
mysales = read.csv("/Users/stef/Desktop/countsFun/data/MySelectedSeries.csv")


# attach the datafrmae
n = 104
Smallsales  = mysales[1:n,]
attach(Smallsales)


# plot the series
# plot.ts(MOVE, ylab = "sales", xlab = "time")

# use GCMR to fit the data
#tic()
#mod <- gcmr(MOVE~Buy, marginal = negbin.marg, cormat = arma.cormat(3, 0), no.se = FALSE)
#toc()

# plot residuals
# plot(mod)


#get the estimates
#est = as.numeric(mod$estimate)


# # use gcmr as initiual ewstimates
# b0 = est[1]
# b1 = est[2]
# b = c(b0,b1)
# k  = est[3]
# phi = est[4:6]
# theta = c(b,k,phi)

# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf)

#other parameters
MaxCdf = 5000
nHC = 30
ARMAorder = c(3,0)
epsilon = 1
ParticleNumber  = 100

# regressor variable with intercept
Regressor = cbind(rep(1,length(Buy)),Buy)


#========================= FIT GEN POISSON ================================#

# lower and upper bounds
LB = c(-100, -100, 0.001, -Inf, -Inf, -Inf)
UB = c(100, 100, 0.999, Inf, Inf, Inf)
theta = c(2, 1, 0.5, -0.3,0.2,0.1)
ARMAorder = c(3,0)
CountDist = "Generalized Poisson"
ParticleNumber  = 100
epsilon = 0.5
OptMethod = "bobyqa"

#- fit via GL
tic()
mod4 = FitGaussianLikGP_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, 400, nHC, OptMethod)
toc()

#- fit via PF
tic()
mod5 = FitMultiplePFResReg(theta, MOVE, Regressor,CountDist, ParticleNumber, LB, UB, ARMAorder, epsilon)
toc()



optim.output4 <- optimx(par       = theta,
                      fn        = GaussLogLikGP_Reg,
                      data      = MOVE,
                      Regressor = Regressor,
                      ARMAorder = ARMAorder,
                      MaxCdf    = 300,
                      nHC       = nHC,
                      hessian   = TRUE,
                      lower     = LB,
                      upper     = UB,
                      control = list(all.methods = TRUE,trace = 2)
)


optim.output5 <- optimx(par            = theta,
                       fn             = ParticleFilterRes_Reg,
                       data           = MOVE,
                       Regressor      = Regressor,
                       ARMAorder      = ARMAorder,
                       ParticleNumber = ParticleNumber,
                       CountDist      = CountDist,
                       epsilon        = epsilon,
                       lower          = LB,
                       upper          = UB,
                       hessian        = TRUE,
                       control = list(all.methods = TRUE,trace = 2))





#========================= FIT Neg Bin  ================================#
# GL with r depending on t
# tic()
# mod1 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC,1)
# toc()
#

# GL with p depending on t
tic()
theta = c(2, 1, 0.5, -0.3,0.2,0.1)
mod2 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, LB, UB, ARMAorder, MaxCdf, nHC,0)
toc()

# PF with p depending on t
CountDist = "Negative Binomial"
theta = c(2, 1, 0.5, -0.3,0.2,0.1)
tic()
mod3 = FitMultiplePFResReg(theta, MOVE, Regressor,CountDist, ParticleNumber, LB, UB, ARMAorder, 0.5)
toc()
#==========================================================================#





optim.output2 <- optimx(par       = theta,
                       fn        = GaussLogLikNB_Reg,
                       data      = MOVE,
                       Regressor = Regressor,
                       ARMAorder = ARMAorder,
                       MaxCdf    = 300,
                       nHC       = nHC,
                       hessian   = TRUE,
                       lower     = LB,
                       upper     = UB,
                       control = list(all.methods = TRUE,trace = 2)
)





# theta = c(1,0, 0.5,0,0,0)
optim.output <- optimx(par            = theta,
                       fn             = ParticleFilterRes_Reg,
                       data           = MOVE,
                       Regressor      = Regressor,
                       ARMAorder      = ARMAorder,
                       ParticleNumber = ParticleNumber,
                       CountDist      = CountDist,
                       epsilon        = epsilon,
                       lower          = LB,
                       upper          = UB,
                       hessian        = TRUE,
                       control = list(all.methods = TRUE,trace = 2))













