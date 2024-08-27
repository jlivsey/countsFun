#---------------------------------------------------------------------------------------------------------#
# PURPOSE: Main Wrapper for the lgc package.
#
#
#
#
# Authors: Stefanos Kechagias, Jiajie Kong, James Livsey, Robert Lund, Vladas Pipiras
#---------------------------------------------------------------------------------------------------------#

# Final wrapper function
lgc = function(DependentVar   = NULL,
               Regressor      = NULL,
               EstMethod      = "PFR",
               CountDist      = NULL,
               ARMAModel      = NULL,
               ParticleNumber = 5,
               epsilon        = 0.5,
               initialParam   = NULL,
               TrueParam      = NULL,
               Task           = 'Evaluation',
               SampleSize     = NULL,
               nsim           = NULL,
               no_cores       = 1,
               OptMethod      = "L-BFGS-B",
               OutputType     = "list",
               ParamScheme    = NULL,
               maxdiff        = 10^(-8),
               ntrials        = NULL,...){

  # parse all the parameters and the data into a list called mod
  mod = ModelScheme(DependentVar   = DependentVar,
                    Regressor      = Regressor,
                    EstMethod      = EstMethod,
                    ARMAModel      = ARMAModel,
                    CountDist      = CountDist,
                    ParticleNumber = ParticleNumber,
                    epsilon        = epsilon,
                    initialParam   = initialParam,
                    TrueParam      = TrueParam,
                    Task           = Task,
                    SampleSize     = SampleSize,
                    OptMethod      = OptMethod,
                    OutputType     = OutputType,
                    ParamScheme    = ParamScheme,
                    maxdiff        = maxdiff,
                    ntrials        = ntrials)

  # if(Task=='Synthesis'){
  #   out
  # }


  # compute initial parameters if they haven't been provided
  if (is.null(mod$initialParam)){
    mod$initialParam = InitialEstimates(mod)
  }

  # if simulation task has been chosen simulate the data and compute initial estimates
  # check me how fast is this?
  if(Task=='Simulation'){
    AllSimulatedSeries <- vector(mode='list', length=nsim)
    AllInitialParam    <- vector(mode='list', length=nsim)
    for (i in 1:nsim) {
      set.seed(i)
      AllSimulatedSeries[[i]] = mod$DependentVar =sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm)
      AllInitialParam[[i]]    = InitialEstimates(mod)
    }

    # renew the task (after we simulated we need to fit)
    mod$Task = 'Optimization'

    # we need some cores in order to fit the data - we ll use all but one
    if(is.null(no_cores)) no_cores = detectCores() - 1

    # initiate and register the cluster
    cl <- makeCluster(no_cores)

    #clusterSetRNGStream(cl, 1001) #make the bootstrapping exactly the same as above to equate computation time
    registerDoParallel(cl)

    # run foreach
    # fix me: need to be very careful here with packages - and run tests for all distributions with and without regressors
    out = foreach(index = 1:nsim,
                .combine = rbind,
                .packages = c("ltsa", "optimx", 'tictoc', 'countsFun'))  %dopar%  {
                  mod$DependentVar =  AllSimulatedSeries[[index]]
                  theta  = mod$initialParam = AllInitialParam[[index]]
                  FitMultiplePF_Res(theta,mod)
                }

    stopCluster(cl)
  }

  if(Task %in% c('Evaluation', 'Optimization')){
      theta  = mod$initialParam
      out = FitMultiplePF_Res(theta, mod)
  }

# create an lgc object and save the initial estimate
class(out) = "lgc"

return(out)

}

