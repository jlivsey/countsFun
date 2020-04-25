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
nsim            = 200
no_cores <- detectCores()/2

#-----------------------------------------------Positive MA parameter--------------------------------------------------#
MAParm = -0.75
ThetaSign = ifelse(MAParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param

n=100
df7 = NegBinMA1_GL(CountDist, MargParm, MAParm, n, nsim, no_cores)
save(df7, file = sprintf("NegBin%s_%sMA%s_GL_N%s_NS%s_Theta%s.RData", MargParm[1], MargParm[2], 1, n, nsim, ThetaSign))

n=200
df8 = NegBinMA1_GL(CountDist, MargParm, MAParm, n, nsim, no_cores)
save(df8, file = sprintf("NegBin%s_%sMA%s_GL_N%s_NS%s_Theta%s.RData", MargParm[1], MargParm[2], 1, n, nsim, ThetaSign))

n=400
df9 = NegBinMA1_GL(CountDist, MargParm, MAParm, n, nsim, no_cores)
save(df9, file = sprintf("NegBin%s_%sMA%s_GL_N%s_NS%s_Theta%s.RData", MargParm[1], MargParm[2], 1, n, nsim, ThetaSign))




#-----------------------------------------------Negative MA parameter--------------------------------------------------#
MAParm = 0.75
ThetaSign = ifelse(MAParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param

n=100
df10 = NegBinMA1_GL(CountDist, MargParm, MAParm, n, nsim, no_cores)
save(df10, file = sprintf("NegBin%s_%sMA%s_GL_N%s_NS%s_Theta%s.RData", MargParm[1], MargParm[2], 1, n, nsim, ThetaSign))

n=200
df11 = NegBinMA1_GL(CountDist, MargParm, MAParm, n, nsim, no_cores)
save(df11, file = sprintf("NegBin%s_%sMA%s_GL_N%s_NS%s_Theta%s.RData", MargParm[1], MargParm[2], 1, n, nsim, ThetaSign))

n=400
df12 = NegBinMA1_GL(CountDist, MargParm, MAParm, n, nsim, no_cores)
save(df12, file = sprintf("NegBin%s_%sMA%s_GL_N%s_NS%s_Theta%s.RData", MargParm[1], MargParm[2], 1, n, nsim, ThetaSign))


