
# Final wrapper function
countC = function(data, Regressor=NULL, CountDist=NULL, EstMethod="PFR", ARMAorder=c(0,0),
                  nHC=30, MaxCdf=500, ParticleNumber=400, epsilon = 0.5, initialParam = NULL,
                  OptMethod = "bobyqa", maxit=100){

  # check validity of arguments
  # rc = CheckInputSpecs(data, Regressor, CountDist, EstMethod, ARMAorder,
  #                            nHC, MaxCdf, ParticleNumber, epsilon, initialParam,OptMethod )
  #
  # if(rc[[1]]) stop(rc[[2]])

  # parse everything in modelInfo
  # FIX ME:check where is maxit used?
  mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon, initialParam, EstMethod, maxit)

  # stop if there was an error in model specification
  if(mod$error) stop(mod$errorMsg)

  # fix me: I need a function that computes initial parameters
  if (is.null(initialParam)){
    theta  = InitialEstimates(mod)
  }else{
    theta  = mod$initialParam
  }

  # fit the model using GL
  if(EstMethod=="GL"){
    out = FitGaussianLogLik(theta, mod)
  }

  # fit the model using PF
  if(EstMethod=="PFR"){
    out = FitMultiplePF_Res(theta, mod)
  }


  # fit the model using PF
  if(EstMethod=="IYW"){
    out = FitMultiplePF_Res(theta, mod)
  }


  return(out)

}

