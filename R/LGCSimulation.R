LGCSimulation = function(nsim           = NULL,
                         SampleSize     = NULL,
                         Regressor      = NULL,
                         EstMethod      = NULL,
                         CountDist      = NULL,
                         MargParm       = NULL,
                         ARParm         = NULL,
                         MAParm         = NULL,
                         ParticleNumber = NULL,
                         epsilon        = NULL,
                         initialParam   = NULL,
                         OptMethod      = NULL,
                         Optimization   = NULL,
                         OutputType     = NULL,
                         ParamScheme    = NULL,
                         no_cores       = NULL)
  {


  # require necessary libraries.
  require(foreach)
  require(doParallel)

  # Simulation scheme details

  if(is.null(no_cores)) no_cores = detectCores() - 1                  # Select the number of cores
  ##-------------------------------------------------------------------------------------------------#

  # generate all the realizations and save in a list
  l <- list()
  for (i in 1:nsim) {
    set.seed(i)
    l[[i]] = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm)
  }

  TrueParam = c(MargParm, ARParm, MAParm)
  ARMAorder = c(length(ARParm), length(MAParm))

  # initiate and register the cluster
  cl <- makeCluster(no_cores)

  #clusterSetRNGStream(cl, 1001) #make the bootstrapping exactly the same as above to equate computation time
  registerDoParallel(cl)

  # run foreach
  all = foreach(index = 1:nsim,
             .combine = rbind,
            .packages = c("ltsa", "optimx", 'tictoc', 'countsFun'))  %dopar%  {
                  lgc(DependentVar   = l[[index]],
                      Regressor      = Regressor,
                      EstMethod      = EstMethod,
                      CountDist      = CountDist,
                      ARMAorder      = ARMAorder,
                      ParticleNumber = ParticleNumber,
                      epsilon        = epsilon,
                      initialParam   = initialParam,
                      TrueParam      = TrueParam,
                      Optimization   = Optimization,
                      OptMethod      = OptMethod,
                      OutputType     = OutputType,
                      ParamScheme    = ParamScheme
                  )
                }

  stopCluster(cl)

  return(all)
}


















