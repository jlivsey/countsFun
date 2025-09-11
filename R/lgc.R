#---------------------------------------------------------------------------------------------------------#
# PURPOSE: Main Wrapper for the lgc package.
#
#
#
#
# Authors: Stefanos Kechagias, Jiajie Kong, James Livsey, Robert Lund, Vladas Pipiras
#---------------------------------------------------------------------------------------------------------#

# Final wrapper function
lgc = function(formula        = NULL,
               data           = NULL,
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

  # parse the regression formula
  parsed_formula <- parse_formula(formula)

  # retrieve the Dependent variable
  DependentVar = data[parsed_formula$DependentVar]

  # retrieve the Regressors variable
  if(is.null(parsed_formula$Regressor)){
    Regressor = NULL
  } else{
    Regressor =   data[parsed_formula$Regressor]
  }


  # retrieve intercept
  Intercept = parsed_formula$intercept

  # add a column of ones in the Regressors if Intercept is present
  if (!is.null(Regressor) && Intercept){
    Regressor = cbind(rep(1,dim(data)[1]),Regressor)
    names(Regressor)[1] = "Intercept"
  }

  # parse all the parameters and the data into a list called mod
  mod = ModelScheme(DependentVar   = DependentVar,
                    Regressor      = Regressor,
                    Intercept      = Intercept,
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

    # retrieve the parameters
    Parms = RetrieveParameters(TrueParam,mod)

    AllSimulatedSeries <- vector(mode='list', length=nsim)
    AllInitialParam    <- vector(mode='list', length=nsim)
    for (i in 1:nsim) {
      set.seed(i)
      AllSimulatedSeries[[i]] = mod$DependentVar =sim_lgc(SampleSize, mod$CountDist, Parms$MargParms, Parms$AR, Parms$MA,
                                                          mod$Regressor,mod$Intercept)
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

    # .packages = c("ltsa", "optimx", 'tictoc', 'countsFun', 'itsmr',
    #               'doParallel','numDeriv','VGAM','iZID','extraDistr','devtools',
    #               'parallel','MASS','mixtools', 'optextras')

    out = foreach(index = 1:nsim,
                .packages = c("optimx", 'countsFun'),.export= c("FitMultiplePF_Res_New"))  %dopar%  {
                  mod$DependentVar =  AllSimulatedSeries[[index]]
                  theta  = mod$initialParam = AllInitialParam[[index]]
                  FitMultiplePF_Res_New(theta,mod)
                }

    stopCluster(cl)
  }

  if(Task %in% c('Evaluation', 'Optimization')){
      theta  = mod$initialParam
      FitResults = FitMultiplePF_Res_New(theta, mod)

      # gather the input information and the Fit Results in one output structure
      out = PrepareOutput(mod, FitResults)
  }

  # create an lgc object and save the initial estimate
  # check me: also for
  if(Task %in% c('Evaluation', 'Optimization')){
    class(out) <- "lgc"
  }

  if(Task=='Simulation'){
    for (i in seq_along(out)) {
      class(out[[i]]) <- "lgc"
    }
  }
return(out)

}

