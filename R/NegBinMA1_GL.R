NegBinMA1_GL = function(CountDist,MargParm,MAParm,
                      n, nsim, no_cores) {

# PURPOSE: Fit Gaussian likelihood to many realizations of synthetic Poisson data with
#          an underlying MA(1). We want to investigate if having a large lambda
#          and large sample size reduce the bias in parameter estimates.
#
# NOTES:
#
# AUTHORS: James Livsey, Stefanos Kechagias, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.5.3

# ---- Load libraries ----
library(parallel)
library(doParallel)
library(countsFun)
library(ggplot2)

# ---- setup parameters for Poisson(lam)-AR(1) series ----

initial.param = c(MargParm, MAParm)    # Initial Parameters


# Generate all the data and save in a list
l <- list()
for(r in 1:nsim){
  set.seed(r)
  l[[r]] = sim_negbin_ma(n, MAParm, MargParm[1], MargParm[2] )
}


# initiate and register the cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# fit the gaussian log lik using foreach
all = foreach(index = 1:nsim,
              .combine = rbind,
              .packages = c("countsFun")) %dopar%
  FitGaussianLikNB_MA(initial.param, l[[index]])

stopCluster(cl)

# Prepare results for the plot.
df = data.frame(matrix(ncol = 11, nrow = nsim))


names(df) = c('r.est',
              'p.est',
              'theta.est',
              'estim.method',
              'n',
              'theta.true',
              'theta.se',
              'r.true',
              'r.se',
              'p.true',
              'p.se' )

df[,1:3] = all[,1:3]
df[,4]   = 'gaussianLik'
df[,5]   = n
df[,6]   = MAParm
df[,7]   = all[,6]
df[,8]   = MargParm[1]
df[,9]   = all[,4]
df[,10]  = MargParm[2]
df[,11]  = all[,5]



return(df)
