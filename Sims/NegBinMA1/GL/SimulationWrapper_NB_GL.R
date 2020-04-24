# PURPOSE: Call the NegBinMA1_GL function that performs the NegBin-MA(1) GL simulations for our paper.
#          Here we call it multiple times for different simulation schemes

# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# set directory to  save results
setwd("C:/Users/Stef/Desktop/countsFun/Sims/NegBinMA1/GL/RData")

# ---- Load libraries ----
library(parallel)
library(doParallel)
library(countsFun)

# fixed parameters across all simulation schemes
CountDist       = "NegBin"
r               = 3
p               = 0.2
MargParm        = c(r,p)
nsim            = 5
no_cores <- detectCores() -1

#-----------------------------------------------Positive MA parameter--------------------------------------------------#
MAParm = 0.75
ThetaSign = ifelse(MAParm > 0, 'Pos', 'Neg')   # SIGN OF MA(1) param

n=100
df7 = NegBinMA1_GL(CountDist, MargParm, MAParm, n, nsim, no_cores)

