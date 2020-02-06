library(countsFun)

# ---- setup parameters for Poisson(lam)-AR(p) series ----
lam = 5
phi = .75
n = 100     # sample size
nsim = 50  # number of realizations

# ---- Variables needed that are functions of manual inputs ----
p = length(phi)  # AR order
nparms = p+1     # total number of parameters

# ---- allocate memory to save the following: ----
ParmEst = matrix(0, nrow=nsim, ncol=nparms)

# for-loop for each realization
for (r in 1:nsim){

  # Print 1st 10 iterations numbers and then every 20th
  if(r < 10) print(r)
  if(r %% 20 == 0) print(r)

  # generate Poisson-AR(p) data
  x <- sim_pois_ar(n, phi, lam )

  # select initial parameters as true ones
  initial.param <- c(lam, phi)

  # run optimization for our model
  optim.output <- FitGaussianLik(initialParam = initial.param, x = x)

  # Store parameter estimates
  ParmEst[r, ] <- optim.output[1:2]

  # reset the seed
  set.seed(r)
}




# ---- Save output to external file ----

# save estimates and std errors in stef MAC
#dir = sprintf('~/Dropbox/latentGaussCounts/simulations/poisson ar1/PFrevision/PoisAR%.0f_N%.0f_Nfit%.0f.csv',p,n,nfit )

# pc directory
#dir = sprintf("C:/Users/stefa/Dropbox/latentGaussCounts/simulations/poisson ar1/PFRevision/PoisAR%.0f_N%.0f_Nfit%.0f.csv",p,n/100,nfit)

# cluster directory
# dir = sprintf('/nas/longleaf/home/kechagia/LCG-PF-Rcode/PoisAR%.0f_NH%.0f_Nfit%.0f_TrueInit.csv',p,n/100,nfit )

# Jim MAC
# dir = sprintf('~/Dropbox/jim/latentGaussCounts/simulations/Revisions/GaussianLik/Poisson/PoisAR%.0f_N%.0f.csv',p,n)
#
# write.csv(All,dir)
