# PURPOSE: Call the NegBinMA1_GL function that performs the NegBin-MA(1) GL simulations for our paper.
#          Here we call it multiple times for different simulation schemes

# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# set directory to  save results
setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/GL/RData")

# ---- Load libraries ----
library(parallel)
library(doParallel)
library(countsFun)
library(optimx)

# fixed parameters across all simulation schemes
CountDist       = "Mixed Poisson"
nsim            = 200
no_cores <- detectCores()-1
LB = c(0.001, 0.01, 0.01, -0.995)
UB = c(0.499, Inf, Inf, 0.995)
MaxCdf = 5000
nHC = 30
p = 1
q = 0
ARParm = 0.75

#-----------------------------------------------Positive MA parameter--------------------------------------------------#

lam1            = 2
lam2            = 5
prob            = 0.25
MargParm        = c(prob, lam1, lam2)
trueParam = c(MargParm, ARParm)


n =100
df7 = MixedPoisson_GL(trueParam, p, q, LB, UB, MaxCdf, nHC, n, nsim, no_cores)
save(df7, file = sprintf("MixedPoisson%s_%s_%sAR%s_GL_N%s_NS%s_NotTrue.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

n =200
df8 = MixedPoisson_GL(trueParam, p, q, LB, UB, MaxCdf, nHC, n, nsim, no_cores)
save(df8, file = sprintf("MixedPoisson%s_%s_%sAR%s_GL_N%s_NS%s_NotTrue.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

n =400
df9 = MixedPoisson_GL(trueParam, p, q, LB, UB, MaxCdf, nHC, n, nsim, no_cores)
save(df9, file = sprintf("MixedPoisson%s_%s_%sAR%s_GL_N%s_NS%s_NotTrue.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))


#-----------------------------------------------Negative MA parameter--------------------------------------------------#
lam1            = 2
lam2            = 10
prob            = 0.25
MargParm        = c(prob, lam1, lam2)
trueParam = c(MargParm, ARParm)

n =100
df10 = MixedPoisson_GL(trueParam, p, q, LB, UB, MaxCdf, nHC, n, nsim, no_cores)
save(df10, file = sprintf("MixedPoisson%s_%s_%sAR%s_GL_N%s_NS%s_NotTrue.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

n =200
df11 = MixedPoisson_GL(trueParam, p, q, LB, UB, MaxCdf, nHC, n, nsim, no_cores)
save(df11, file = sprintf("MixedPoisson%s_%s_%sAR%s_GL_N%s_NS%s_NotTrue.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

n =400
df12 = MixedPoisson_GL(trueParam, p, q, LB, UB, MaxCdf, nHC, n, nsim, no_cores)
save(df12, file = sprintf("MixedPoisson%s_%s_%sAR%s_GL_N%s_NS%s_NotTrue.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))














# ARMAorder = c(1,0)
# # # list with true ARMA parameters
# ARMAmodel = list(NULL,NULL)
# if(p>0){ARMAmodel[[1]] = trueParam[4:(3+p)]}
# if(q>0){ARMAmodel[[2]] = trueParam[(4+p):length(trueParam)]}
# set.seed(1)
# data = sim_mixedPoisson(n, ARMAmodel, MargParm[1], MargParm[2], MargParm[3] )
# theta = trueParam
# x0 = theta
# X = data
#GaussLogLikMP(theta, data, ARMAorder, MaxCdf, nHC)

#
#
#
# optim.output <- optimx(par       = x0,
#                       fn        = GaussLogLikMP,
#                       data      = X,
#                       ARMAorder = ARMAorder,
#                       MaxCdf    = MaxCdf,
#                       nHC       = nHC,
#                       method    = "L-BFGS-B",
#                       hessian   = TRUE,
#                       lower     = LB,
#                       upper     = UB
# )
#
# ParmEst   = c(optim.output$p1,optim.output$p2,optim.output$p3,optim.output$p4)
# H = gHgen(par   = ParmEst,
#           fn    = GaussLogLikMP,
#           data  = X,
#           ARMAorder = ARMAorder,
#           MaxCdf    = MaxCdf,
#           nHC       = nHC
# )$Hn
#
#
# sqrt(abs(diag(solve(H))))







