PoisAR1_GL = function(CountDist,MargParm,ARParm, n, nsim, no_cores) {

  # PURPOSE: Wrapper that performs simulation and produces Poisson AR(1) Gaussian likelihood
  # estimates for revised Figure 3 (see aarxiv)
  #
  #
  # AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
  #
  # DATE:    April 2020
  #
  # R version 3.6.3


# ---- setup parameters for Poisson(lam)-AR(1) series ----
PhiSign = ifelse(ARParm > 0, 'Pos', 'Neg')  # SIGN OF ar(1) param
ARorder = length(ARParm)                    # AR parameters
initial.param = c(MargParm, ARParm)         # Initial PArameters

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
df = data.frame(matrix(ncol = 8, nrow = nsim))

#Create columns lam.est, phi.est, estim.method, n, phi, phi.se, lam, lam.se
names(df) = c('lam.est', 'phi.est', 'estim.method', 'n', 'phi.true', 'phi.se', 'lam.true', 'lam.se')

df[,1:2] = all[,1:2]
df[,3] = 'gaussianLik'
df[,4] = n
df[,5] = ARParm
df[,6] = all[,4]
df[,7] = MargParm
df[,8] = all[,3]

return(df)
}



