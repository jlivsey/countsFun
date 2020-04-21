# PURPOSE: Produce Poisson AR(1) IWY estimates for revised Figure 3 (see aarxiv)
#
#
# AUTHORS: Stefanos Kechagias, James Livsey
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
ARParm  = -0.75
PhiSign = ifelse(ARParm > 0, 'Pos', 'Neg')  # SIGN OF ar(1) param
ARorder = length(ARParm)                    # AR parameters
n = 400                                     # sample size
nsim = 200                                  # number of realizations
initial.param = c(MargParm, ARParm)         # Initial PArameters
no_cores <- detectCores() -1                # Select the number of cores

# Print File name
fileName <- sprintf("PoisAR%s_IWY_N%s_NS%s_Phi%s", ARorder, n, nsim, PhiSign)
print(fileName)

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
  fit_IWY(x = l[[index]], p = ARorder)

stopCluster(cl)

# Prepare results for the plot.
df = data.frame(matrix(ncol = 8, nrow = 200))

#Create columns lam.est, phi.est, estim.method, n, phi, phi.se, lam, lam.se
names(df) = c('lam.est', 'phi.est', 'estim.method', 'n', 'phi.true', 'phi.se', 'lam.true', 'lam.se')

df[,1:2] = all[,1:2]
df[,3] = 'IYW'
df[,4] = n
df[,5] = ARParm
df[,6] = NA
df[,7] = MargParm
df[,8] = NA

# save results
df15 <- df
saveFileName <- paste0(fileName, ".RData")
save(df15, file = saveFileName)
rm(list = ls())



