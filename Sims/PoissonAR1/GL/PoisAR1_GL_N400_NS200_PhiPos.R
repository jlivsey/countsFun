# PURPOSE: Produce Poisson AR(1) Gaussian likelihood estimates for revised Figure 3 (see aarxiv)
#
#
# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# ---- Load libraries ----
library(parallel)
library(doParallel)
library(countsFun)

# ---- setup parameters for Poisson(lam)-AR(1) series ----
CountDist = "Poisson"                       # Distribution
MargParm = 2                                # marginal parameter
ARParm  = 0.75
PhiSign = ifelse(ARParm > 0, 'Pos', 'Neg')  # SIGN OF ar(1) param
ARorder = length(ARParm)                    # AR parameters
n = 400                                     # sample size
nsim = 200                                  # number of realizations
initial.param = c(MargParm, ARParm)         # Initial PArameters
no_cores <- detectCores() -1                # Select the number of cores

# Generate all the data and save in a list
l <- list()
for(r in 1:nsim){
  set.seed(r)
  l[[r]] = sim_pois_ar(n, ARParm, MargParm )
}

# initiate and register the cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# fit the gaussian log likelihood using foreach
all = foreach(index = 1:nsim,
              .combine = rbind,
              .packages = c("MASS", "countsFun")) %dopar%
  FitGaussianLik(initial.param, l[[index]])

stopCluster(cl)

# Prepare results for the plot.
df11 = data.frame(matrix(ncol = 8, nrow = 200))

#Create columns lam.est, phi.est, estim.method, n, phi, phi.se, lam, lam.se
names(df11) = c('lam.est', 'phi.est', 'estim.method', 'n', 'phi.true', 'phi.se', 'lam.true', 'lam.se')

df11[,1:2] = all[,1:2]
df11[,3] = 'gaussianLik'
df11[,4] = n
df11[,5] = ARParm
df11[,6] = all[,4]
df11[,7] = MargParm
df11[,8] = all[,3]


# save results
save(df11, file = sprintf("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/GL/PoisAR%s_GL_N%s_NS%s_Phi%s.RData", ARorder, n, nsim, PhiSign))
rm(list = ls())



