
# Final wrapper function
countC = function(data, Regressor=NULL, CountDist="Poisson", EstMethod="PFR", ARMAorder=c(0,0),
                  nHC=30, MaxCdf=500, ParticleNumber=400, epsilon = 0.5,initialParam = NULL,
                  LB=0.01, UB=Inf, OptMethod = "bobyqa"){

  # fix me: I need a function that checks validity of arguments here

  # retrieve model specifications
  mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )


  # fix me: I need a function that returns default constraints


  # fix me: I need a function that computes initial parameters

  # fit the model using GL
  if(EstMethod=="GL"){
    out = FitMultiplePF_Res(theta, data, Regressor, mod, LB, UB, OptMethod)
  }

  # fit the model using PF
  if(EstMethod=="PF"){
    out = FitMultiplePF_Res(theta, data, Regressor, mod, LB, UB, OptMethod)
  }



  return(out)

}

