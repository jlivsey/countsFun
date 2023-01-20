LGCSimulation = function(nsim           = nsim,
                         SampleSize     = NULL,
                         Regressor      = NULL,
                         EstMethod      = "PFR",
                         CountDist      = NULL,
                         MargParm       = NULL,
                         ARParm         = NULL,
                         MAParm         = NULL,
                         ParticleNumber = 10,
                         epsilon        = 0.5,
                         initialParam   = NULL,
                         OptMethod      = "bobyqa",
                         ParamScheme    = NULL)
  {


  # require necessary libraries.
  require(foreach)
  require(doParallel)

  # Simulation scheme details
  no_cores      = detectCores() - 1                  # Select the number of cores
  ##-------------------------------------------------------------------------------------------------#

  # generate all the realizations and save in a list
  l <- list()
  for (i in 1:nsim) {
    set.seed(i)
    l[[i]] = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm)
  }

  ARMAorder = c(length(ARParm), length(MAParm))

  # initiate and register the cluster
  cl <- makeCluster(no_cores)

  #clusterSetRNGStream(cl, 1001) #make the bootstrapping exactly the same as above to equate computation time
  registerDoParallel(cl)

  # run foreach
  all = foreach(index = 1:nsim,
             .combine = rbind,
            .packages = c("ltsa", "optimx"))  %dopar%  {
                  lgc(DependentVar   = l[[index]],
                      Regressor      = Regressor,
                      EstMethod      = EstMethod,
                      CountDist      = CountDist,
                      ARMAorder      = ARMAorder,
                      ParticleNumber = ParticleNumber,
                      epsilon        = epsilon,
                      initialParam   = initialParam,
                      OptMethod      = OptMethod
                  )
                }

  stopCluster(cl)


  # # Prepare results for the plot.
  # df = data.frame(matrix(ncol = 8, nrow = nsim))
  #
  # #Create columns lam.est, phi.est, estim.method, n, phi, phi.se, lam, lam.se
  # names(df) = c(
  #   'lam.est',
  #   'phi.est',
  #   'estim.method',
  #   'n',
  #   'phi.true',
  #   'phi.se',
  #   'lam.true',
  #   'lam.se'
  # )
  #
  # df[, 1:2] = all[, 1:2]
  # df[, 3] = 'particle'
  # df[, 4] = n
  # df[, 5] = ARParm
  # df[, 6] = all[, 4]
  # df[, 7] = MargParm
  # df[, 8] = all[, 3]

  return(all)
}
