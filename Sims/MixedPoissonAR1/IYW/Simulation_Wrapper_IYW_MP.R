# PURPOSE: Call the NegBinAR1_GL function that performs the NegBin-MA(1) GL simulations for our paper.
#          Here we call it multiple times for different simulation schemes

# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# set directory to  save results
setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/IYW/RData")
# ---- Load libraries ----
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
lam1        = 2
prob        = 0.25
#-----------------------------------------------lam2 = 5--------------------------------------------------#

lam2      = 5
MargParm  = c(prob, lam1, lam2)
trueParam = c(MargParm, ARParm)

n = 100
df13 = MixedPoisson_IYW(trueParam, p, q, LB, UB, n, nsim, nHC)
save(df13, file = sprintf("MixedPoisson%s_%s_%sAR%s_IYW_N%s_NS%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

n = 200
df14 = MixedPoisson_IYW(trueParam, p, q, LB, UB, n, nsim, nHC)
save(df14, file = sprintf("MixedPoisson%s_%s_%sAR%s_IYW_N%s_NS%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

n = 400
df15 = MixedPoisson_IYW(trueParam, p, q, LB, UB, n, nsim, nHC)
save(df15, file = sprintf("MixedPoisson%s_%s_%sAR%s_IYW_N%s_NS%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))


#-----------------------------------------------lam2 = 10--------------------------------------------------#

lam2      = 10
MargParm  = c(prob, lam1, lam2)
trueParam = c(MargParm, ARParm)

n = 100
df16 = MixedPoisson_IYW(trueParam, p, q, LB, UB, n, nsim, nHC)
save(df16, file = sprintf("MixedPoisson%s_%s_%sAR%s_IYW_N%s_NS%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

n = 200
df17 = MixedPoisson_IYW(trueParam, p, q, LB, UB, n, nsim, nHC)
save(df17, file = sprintf("MixedPoisson%s_%s_%sAR%s_IYW_N%s_NS%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

n = 400
df18 = MixedPoisson_IYW(trueParam, p, q, LB, UB, n, nsim, nHC)
save(df18, file = sprintf("MixedPoisson%s_%s_%sAR%s_IYW_N%s_NS%s.RData", MargParm[1], MargParm[2],MargParm[3], 1, n, nsim))

