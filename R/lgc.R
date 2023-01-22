
# Final wrapper function
lgc = function(DependentVar, Regressor=NULL, EstMethod="PFR", CountDist=NULL, ARMAorder=c(0,0),
                   ParticleNumber=400, epsilon = 0.5, initialParam = NULL, TrueParam = NULL,
                   Optimization = TRUE, OptMethod = "bobyqa", OutputType="list", ParamScheme = NULL){

  # parse all the parameters and the data into a list called mod
  mod = ModelScheme(DependentVar, Regressor, EstMethod, ARMAorder, CountDist,ParticleNumber, epsilon,
                    initialParam, TrueParam, Optimization, OptMethod, OutputType, ParamScheme)

  # fix me: we need to do implement error checking
  if(mod$error) stop(mod$errorMsg)

  # compute initial parameters if they haven't been provided and save them in mod
  if (is.null(initialParam)){
    theta  = InitEst =  InitialEstimates(mod)
    mod$initialParam = InitEst
  }else{
    theta  = InitEst = mod$initialParam
  }

  # fit the model using PF
  if(EstMethod=="PFR"){
    out = FitMultiplePF_Res(theta, mod)
    return(out)
  }

  # create an lgc object and save the initial estimate
  class(out) = "lgc"
  if(OutputType == "list"){
    out$InitialEstim = InitEst
  }
  return(out)

}

