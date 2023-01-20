# PURPOSE: Call the PoisAR1_GL function that performs the Poisson-Ar(1) GL simulations for our paper.
#          Here we call it multiple times for different simulation schemes

# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# set directory to  save results
setwd("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/GL/RData")

# ---- Load libraries ----
library(parallel)
library(doParallel)
library(countsFun)

# fixed parameters across all simulation schemes
CountDist       = "Poisson"
MargParm        = 10
nsim            = 200
no_cores <- detectCores() -1

#-----------------------------------------------Positive AR parameter--------------------------------------------------#
ARParm = 0.75
PhiSign = ifelse(ARParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param

n=100
df7 = PoisAR1_GL(CountDist, MargParm, ARParm, n, nsim, no_cores)
save(df7, file = sprintf("Pois%sAR%s_GL_N%s_NS%s_Phi%s.RData", MargParm, 1, n, nsim, PhiSign))

n=200
df8 = PoisAR1_GL(CountDist, MargParm, ARParm, n, nsim, no_cores)
save(df8, file = sprintf("Pois%sAR%s_GL_N%s_NS%s_Phi%s.RData", MargParm, 1, n, nsim, PhiSign))

n=400
df9 = PoisAR1_GL(CountDist, MargParm, ARParm, n, nsim, no_cores)
save(df9, file = sprintf("Pois%sAR%s_GL_N%s_NS%s_Phi%s.RData", MargParm, 1, n, nsim, PhiSign))


#-----------------------------------------------Negative AR parameter--------------------------------------------------#
ARParm = -0.75
PhiSign = ifelse(ARParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param

n=100
df10 = PoisAR1_GL(CountDist, MargParm, ARParm, n, nsim, no_cores)
save(df10, file = sprintf("Pois%sAR%s_GL_N%s_NS%s_Phi%s.RData", MargParm, 1, n, nsim, PhiSign))

n=200
df11 = PoisAR1_GL(CountDist, MargParm, ARParm, n, nsim, no_cores)
save(df11, file = sprintf("Pois%sAR%s_GL_N%s_NS%s_Phi%s.RData", MargParm, 1, n, nsim, PhiSign))

n=400
df12 = PoisAR1_GL(CountDist, MargParm, ARParm, n, nsim, no_cores)
save(df12, file = sprintf("Pois%sAR%s_GL_N%s_NS%s_Phi%s.RData", MargParm, 1, n, nsim, PhiSign))




