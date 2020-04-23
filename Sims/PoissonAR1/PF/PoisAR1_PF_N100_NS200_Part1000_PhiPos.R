# PURPOSE: Produce Poisson AR(1) PArticle filter estimates for revised Figure 3 (see aarxiv)
#
#
# AUTHORS: Stefanos Kechagias, James Livsey
#
# DATE:    April 2020
#
# R version 3.6.3


# load PF likelihood function and other necessary functions from UNC cluster directory
setwd("C:/Users/Stef/Desktop/countsFun/R")
source('LikSISGenDist_ARp_Res.R')
source('LikSIS_ARpGenDist_functions.R')

# load necessary libraries.
library(itsmr)
library(FitAR)
library(foreach)
library(doParallel)
library(tictoc)

# Setup the simulation scheme
CountDist = "Poisson"                        # Distribution
MargParm = 2                                 # marginal parameter
ARParm = 0.75                                # AR parameters
PhiSign = ifelse(ARParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param
ARorder = length(ARParm)                     # AR parameters
n = 100                                      # sample size
nsim = 200                                   # number of realizations
nfit = 1                                     # number of times that we fit the same realization
ParticleSchemes = 1000                        # number of particles used in likelihood approximation
initial.param = c(MargParm, ARParm)          # Initial PArameters
no_cores <- detectCores() - 1                # Select the number of cores
###########################################################################

# generate all the realizations and save in a list
l <- list()
for(i in 1:nsim) {
  set.seed(i)
  l[[i]] = sim_pois_ar(n, ARParm, MargParm )
}

t0 = tic()
# initiate and register the cluster
cl <- makeCluster(no_cores)

#clusterSetRNGStream(cl, 1001) #make the bootstrapping exactly the same as above to equate computation time
registerDoParallel(cl)

# run foreach
all = foreach(index = 1:nsim,
              .combine = rbind,
              .packages = "FitAR")  %dopar%  {
                FitMultiplePF(initial.param, l[[index]], CountDist, nfit, ParticleSchemes)
              }

stopCluster(cl)
toc(t0)


# Prepare results for the plot.
df1 = data.frame(matrix(ncol = 8, nrow = 200))

#Create columns lam.est, phi.est, estim.method, n, phi, phi.se, lam, lam.se
names(df1) = c('lam.est', 'phi.est', 'estim.method', 'n', 'phi.true', 'phi.se', 'lam.true', 'lam.se')

df1[,1:2] = all[,1:2]
df1[,3] = 'particle'
df1[,4] = n
df1[,5] = ARParm
df1[,6] = all[,4]
df1[,7] = MargParm
df1[,8] = all[,3]


# save results
save(df1, file = sprintf("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/PF/PoisAR%s_PF_N%s_NS%s_Part%s_Phi%s.RData", ARorder, n, nsim,ParticleSchemes, PhiSign))
rm(list = ls())









