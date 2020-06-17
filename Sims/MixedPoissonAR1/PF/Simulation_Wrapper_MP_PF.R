# PURPOSE: Call the NegBinAR1_GL function that performs the NegBin-MA(1) GL simulations for our paper.
#          Here we call it multiple times for different simulation schemes

# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# set directory to  save results
setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/PF/RData")
# ---- Load libraries ----
library(parallel)
library(doParallel)
library(countsFun)
library(optimx)

# fixed parameters across all simulation schemes
CountDist   = "Mixed Poisson"
nsim        = 200
no_cores    = detectCores()-1
LB          = c(0.001, 0.01, 0.01, -0.995)
UB          = c(0.499, Inf, Inf, 0.995)
MaxCdf      = 5000
nHC         = 30
p           = 1
q           = 0
ARParm      = 0.75
Particles   = 100
epsilon     = 0.5
lam1        = 2
prob        = 0.25
useTrueInit = 0
#-----------------------------------------------lam2 = 5--------------------------------------------------#

lam2      = 5
MargParm  = c(prob, lam1, lam2)
trueParam = c(MargParm, ARParm)

n = 100
df1 = MixedPoisson_PF(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0, nHC, useTrueInit)
save(df1, file = sprintf("MixedPoisson%s_%s_%sAR%s_PF_N%s_NS%s_Part%s_e%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim, Particles, epsilon))

n = 200
df2 = MixedPoisson_PF(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0, nHC, useTrueInit)
save(df2, file = sprintf("MixedPoisson%s_%s_%sAR%s_PF_N%s_NS%s_Part%s_e%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim, Particles, epsilon))

n = 400
df3 = MixedPoisson_PF(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0, nHC, useTrueInit)
save(df3, file = sprintf("MixedPoisson%s_%s_%sAR%s_PF_N%s_NS%s_Part%s_e%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim, Particles, epsilon))


#-----------------------------------------------lam2 = 10--------------------------------------------------#
lam2      = 10
MargParm  = c(prob, lam1, lam2)
trueParam = c(MargParm, ARParm)

n = 100
df4 = MixedPoisson_PF(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0, nHC, useTrueInit)
save(df4, file = sprintf("MixedPoisson%s_%s_%sAR%s_PF_N%s_NS%s_Part%s_e%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim, Particles, epsilon))

n = 200
df5 = MixedPoisson_PF(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0, nHC, useTrueInit)
save(df5, file = sprintf("MixedPoisson%s_%s_%sAR%s_PF_N%s_NS%s_Part%s_e%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim, Particles, epsilon))

n = 400
df6 = MixedPoisson_PF(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0, nHC, useTrueInit)
save(df6, file = sprintf("MixedPoisson%s_%s_%sAR%s_PF_N%s_NS%s_Part%s_e%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim, Particles, epsilon))


