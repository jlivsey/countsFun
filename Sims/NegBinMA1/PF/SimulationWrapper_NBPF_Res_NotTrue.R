# PURPOSE: Call the NegBinAR1_GL function that performs the NegBin-MA(1) GL simulations for our paper.
#          Here we call it multiple times for different simulation schemes

# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# set directory to  save results
setwd("C:/Users/Stef/Desktop/countsFun/Sims/NegBinMA1/PF/RData")

# ---- Load libraries ----
library(parallel)
library(doParallel)
library(countsFun)

# fixed parameters across all simulation schemes
CountDist       = "Negative Binomial"
r               = 3
prob            = 0.8
MargParm        = c(r,prob)
nsim            = 200
no_cores <- detectCores()-1
LB = c(0.01, 0.01, -0.99)
UB = c(10, 0.99, 0.99)
p = 0
q = 1
ARMAorder = c(p,q)
Particles = 10
epsilon = 1
nHC = 30
useTrueInit = 1
UseDEOptim = 1

MAParm = -0.75
ThetaSign = ifelse(MAParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param
trueParam = c(MargParm,MAParm)

n = 400
df6 = NegBin_PF(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon,  no_cores, UseDEOptim, nHC, useTrueInit)
save(df6, file = sprintf("NegBin%s_%sMA%s_PF_N%s_NS%s_Part%s_Theta%s_e%s_Deopt.RData", MargParm[1], MargParm[2], 1, n, nsim,Particles, ThetaSign, epsilon))





MAParm = 0.75
ThetaSign = ifelse(MAParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param
trueParam = c(MargParm,MAParm)

n = 400
df3 = NegBin_PF(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon,  no_cores, UseDEOptim, nHC, useTrueInit)
save(df3, file = sprintf("NegBin%s_%sMA%s_PF_N%s_NS%s_Part%s_Theta%s_e%s_Deopt.RData", MargParm[1], MargParm[2], 1, n, nsim,Particles, ThetaSign, epsilon))
