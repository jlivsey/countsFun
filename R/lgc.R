
# Final wrapper function
lgc = function(DependentVar, Regressor=NULL, EstMethod="PFR", CountDist=NULL, ARMAorder=c(0,0),
                   ParticleNumber=400, epsilon = 0.5, initialParam = NULL, TrueParam = NULL,
                   Optimization = TRUE, OptMethod = "bobyqa", OutputType="list", ParamScheme = NULL){


  mod = ModelScheme(DependentVar, Regressor, EstMethod, ARMAorder, CountDist,ParticleNumber, epsilon,
                    initialParam, TrueParam, Optimization, OptMethod, OutputType, ParamScheme)

  # stop if there was an error in model specification
  if(mod$error) stop(mod$errorMsg)

  # fix me: I need a function that computes initial parameters
  if (is.null(initialParam)){
    theta  = InitEst =  InitialEstimates(mod)
  }else{
    theta  = InitEst = mod$initialParam
  }

  # fit the model using PF
  if(EstMethod=="PFR"){
    out = FitMultiplePF_Res(theta, mod)
    return(out)
  }

  class(out) = "lgc"
  if(OutputType == "list"){
    out$InitialEstim = InitEst
  }
  return(out)

}

