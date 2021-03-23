PoisAR1_PF = function(CountDist,MargParm,ARParm,
                       n, nsim, ParticleSchemes) {

  # PURPOSE: Wrapper that performs simulation and produces Poisson AR(1) Particle filter
  # estimates for revised Figure 3 (see aarxiv)
  #
  #
  # AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
  #
  # DATE:    April 2020
  #
  # R version 3.6.3


  # load PF likelihood function and other necessary functions from UNC cluster directory

  #source('C:/Users/Stef/Desktop/countsFun/R/LikSISGenDist_ARp_Res.R')
  #source('C:/Users/Stef/Desktop/countsFun/R/LikSIS_ARpGenDist_functions.R')

  # load necessary libraries.
  library(itsmr)
  library(FitAR)
  library(foreach)
  library(doParallel)
  library(tictoc)

  # Simulation scheme details
  PhiSign       = ifelse(ARParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param
  ARorder       = length(ARParm)                     # AR parameters
  nfit          = 1                                  # number of times that we fit the same realization
  initial.param = c(MargParm, ARParm)                # Initial PArameters
  no_cores      = detectCores() - 1                  # Select the number of cores
  ##-------------------------------------------------------------------------------------------------#

  # generate all the realizations and save in a list
  l <- list()
  for (i in 1:nsim) {
    set.seed(i)
    l[[i]] = sim_pois_ar(n, ARParm, MargParm)
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
  df = data.frame(matrix(ncol = 8, nrow = nsim))

  #Create columns lam.est, phi.est, estim.method, n, phi, phi.se, lam, lam.se
  names(df) = c(
    'lam.est',
    'phi.est',
    'estim.method',
    'n',
    'phi.true',
    'phi.se',
    'lam.true',
    'lam.se'
  )

  df[, 1:2] = all[, 1:2]
  df[, 3] = 'particle'
  df[, 4] = n
  df[, 5] = ARParm
  df[, 6] = all[, 4]
  df[, 7] = MargParm
  df[, 8] = all[, 3]

  return(df)
}
