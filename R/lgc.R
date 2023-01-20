
# Final wrapper function
lgc = function(DependentVar, Regressor=NULL, EstMethod="PFR", CountDist=NULL, ARMAorder=c(0,0),
                   ParticleNumber=400, epsilon = 0.5, initialParam = NULL,
                  OptMethod = "bobyqa"){


  mod = ModelScheme(DependentVar, Regressor, EstMethod, ARMAorder, CountDist,ParticleNumber, epsilon, initialParam)

  # stop if there was an error in model specification
  if(mod$error) stop(mod$errorMsg)

  # fix me: I need a function that computes initial parameters
  # if (is.null(initialParam)){
  #   theta  = InitEst =  InitialEstimates(mod)
  # }else{
  #   theta  = InitEst = mod$initialParam
  # }
  theta  = InitEst = mod$initialParam= c(10,0.75)

  # fit the model using PF
  if(EstMethod=="PFR"){
    out = FitMultiplePF_Res(theta, mod)
  }


  class(out) = "lgc"
  out$InitialEstim = InitEst

  return(out)

}

