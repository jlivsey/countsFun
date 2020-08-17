
# Final wrapper function
countC = function(data, Regressor=NULL, CountDist=NULL, EstMethod="PFR", ARMAorder=c(0,0),
                  nHC=30, MaxCdf=500, ParticleNumber=400, epsilon = 0.5, initialParam = NULL,
                  OptMethod = "bobyqa"){

  # check validity of arguments
  # rc = CheckInputSpecs(data, Regressor, CountDist, EstMethod, ARMAorder,
  #                            nHC, MaxCdf, ParticleNumber, epsilon, initialParam,OptMethod )
  #
  # if(rc[[1]]) stop(rc[[2]])

  # retrieve model specifications
  mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon, initialParam)

  # stop if there was an error in model specification
  if(mod$error) stop(mod$errorMsg)

  # fix me: I need a function that computes initial parameters
  if (is.null(initialParam)){
    theta  = InitialEstimates(data, Regressor, mod)
  }else{
    theta = initialParam
  }

  # fit the model using GL
  if(EstMethod=="GL"){
    out = FitGaussianLogLik(theta, data, Regressor, mod, OptMethod)
  }

  # fit the model using PF
  if(EstMethod=="PFR"){
    out = FitMultiplePF_Res(theta, data, Regressor, mod, OptMethod)
  }


  # fit the model using PF
  if(EstMethod=="IYW"){
    out = FitMultiplePF_Res(theta, data, Regressor, mod, OptMethod)
  }


  return(out)

}


# check validity of input arguments
CheckInputSpecs = function(data, Regressor, CountDist, EstMethod, ARMAorder,
                           nHC, MaxCdf, ParticleNumber, epsilon, initialParam,OptMethod ){

  rc = 0

  # Truncation of cdf
  if (MaxCdf<0) rc = 1

  # check distributions
  if ( !(EstMethod %in%  c("PFR", "GL", "IYW")))  rc = 2

  # check Particle number
  if (EstMethod=="PFR" && (epsilon > 1 || epsilon<0)) rc = 3

  # check Particle number
  if (EstMethod=="PFR" && ParticleNumber<1) rc = 4

  # check distributions
  if ( !(CountDist %in%  c("Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial" ))) rc = 5

  # check ARMAorder
  if (prod(ARMAorder)<0 || length(ARMAorder)!= 2) rc = 6

  # Mixed ARMA model
  if (mod$ARMAorder[1]>0 && mod$ARMAorder[2]>0) rc = 7

  # check data
  if (is.null(data) ||  length(data)<3) rc = 8

  errors = list(
    e1 = 'Please specify a nonnegative MaxCdf.',
    e2 = 'The argument EstMethod must take one of the following values: "GL", IYW","PFR".',
    e3 = 'Please select a value between 0 and 1 for epsilon.',
    e4 = 'Please select a nonegative value for the argument ParticleNumber.',
    e5 = 'The argument CountDist must take one of the following values:
         "Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial".',
    e6 = 'The argument ARMAorder must have length 2 and can not take negative values.',
    e7 = 'Please specify a pure AR or a pure MA model. ARMA(p,q) models with p>0 and q>0 have not yet been implemented.',
    e8 = 'Data must be non null with sample size greater than 3.'
  )


  out = list(
  rc  = rc,
  e   = errors[[rc]]
  )
  return(out)

}




# compute initial estimates
InitialEstimates = function(data, Regressor, mod){

  est  = rep(NA, mod$nMargParms+sum(mod$ARMAorder))

  #-----Poisson case
  if(mod$CountDist=="Poisson"){
    if(mod$nreg==0){
      est[1] = mean(data)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm(data~Regressor, family = poisson)[1]$coef)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }


  if(mod$CountDist=="Negative Binomial"){
    if(mod$nreg==0){
          xbar = mean(data)
          sSquare = var(data)

          # Method of Moments for negBin
          rEst = xbar^2/(sSquare - xbar)
          pEst = 1 - xbar/sSquare
          est[1:2] = c(rEst, pEst)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(data~Regressor)[1]$coef)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }



  if(mod$CountDist=="Mixed Poisson"){
    if(mod$nreg==0){
      # pmle for marginal parameters
      MixPois_PMLE <- pmle.pois(x,2)

      pEst  = MixPois_PMLE[[1]][1]
      l1Est = MixPois_PMLE[[2]][1]
      l2Est = MixPois_PMLE[[2]][2]


      # correct estimates if they are outside the feasible region
      if(pEst<LB[1]){pEst = 1.1*mod$LB[1]}
      if(pEst>UB[1]){pEst = 0.9*mod$UB[1]}

      if(l1Est<LB[2]){l1Est = 1.1*mod$LB[2]}
      if(l2Est<LB[3]){l2Est = 1.1*mod$LB[3]}

      est[1:3] = c(l1Est, l1Est, pEst)

      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(data~Regressor)[1]$coef)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }



  #-----Generalized Poisson case
  if(mod$CountDist=="Generalized Poisson"){
    if(mod$nreg==0){
      xbar = mean(data)
      sSquare = var(data)

      # Method of Moments for negBin
      rEst = xbar^2/(sSquare - xbar)
      pEst = 1 - xbar/sSquare
      est[1:2] = c(rEst, pEst)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(data~Regressor)[1]$coef)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }





  return(est)
}



