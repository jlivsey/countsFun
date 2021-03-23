

# Compute PIT values
PITvalues = function(H, predDist){
  PITvalues = rep(0,H)

  predd1 = predDist[1,]
  predd2 = predDist[2,]
  Tr = length(predd1)

  for (h in 1:H){
    id1 = (predd1 < h/H)*(h/H < predd2)
    id2 = (h/H >= predd2)
    tmp1 = (h/H-predd1)/(predd2-predd1)
    tmp1[!id1] = 0
    tmp2 = rep(0,Tr)
    tmp2[id2] = 1
    PITvalues[h] = mean(tmp1+tmp2)
  }
PITvalues = c(0,PITvalues)
return(diff(PITvalues))
}

# compute predictive distribution AR1--new   file no  resampling
PDvaluesAR1 = function(theta, data, Regressor, ARMAorder, ParticleNumber,CountDist){

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments

  # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1:(1+nreg),
                           "Negative Binomial"   = 1:(2+nreg),
                           "Mixed Poisson"       = 1:(3+nreg),
                           "Generalized Poisson" = 1:(2+nreg),
                           "Binomial"            = 1:(2+nreg))

  if(nreg<1){
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = ppois,
                   "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
                   "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = pGpois,
                   "Binomial"            = pbinom
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = dGpois,
                   "Binomial"            = dbinom
    )
  }else{
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
    )
    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ dnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ dGpois  (x, ConstTheta, DynamTheta)}
    )
  }


  # retrieve marginal distribution parameters
  MargParms  = theta[MargParmIndices]
  nMargParms = length(MargParms)
  nparms     = length(theta)

  # check if the number of parameters matches the model setting
  if(nMargParms + sum(ARMAorder)!=nparms) stop('The length of theta does not match the model specification.')


  # Add the regressor to the parameters--only works for Negative Binomial
  if(CountDist == "Negative Binomial" && nreg>0){
    beta = MargParms[1:(nreg+1)]
    k = MargParms[nreg+2]
    m = exp(Regressor%*%beta)
    r = 1/k
    p = k*m/(1+k*m)
    ConstMargParm = r
    DynamMargParm = p
  }

  if(CountDist == "Generalized Poisson" && nreg>0){
    beta   = MargParms[1:(nreg+1)]
    ConstMargParm  = MargParms[nreg+2]
    DynamMargParm = exp(Regressor%*%beta)
  }


  # retrieve ARMA parameters
  AR = NULL
  if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]

  MA = NULL
  if(ARMAorder[2]>0) MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]


  #set.seed(1)
  xt = data
  T1 = length(xt)
  N  = ParticleNumber # number of particles
  preddist = matrix(0,2,T1-1) # to collect the values of predictive distribution of interest

  # Compute integral limits
  if(nreg==0){
    a = rep( qnorm(mycdf(data[1]-1,t(MargParms)),0,1), N)
    b = rep( qnorm(mycdf(data[1],t(MargParms)),0,1), N)
  }else{
    a = rep( qnorm(mycdf(data[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
    b = rep( qnorm(mycdf(data[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
  }

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
  zhat = AR*zprev

  wprev = rep(1,N)
  rt = sqrt(1-AR^2)

  for (t in 2:T1){
    temp = rep(0,(xt[t]+1))
    for (x in 0:xt[t]){
      # Compute integral limits
      if(nreg==0){
        a = (qnorm(mycdf(x-1,t(MargParms)),0,1) - zhat)/rt
        b = (qnorm(mycdf(x,t(MargParms)),0,1) - zhat)/rt
      }else{
        a = (qnorm(mycdf(x-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
        b = (qnorm(mycdf(x,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
      }
      temp[x+1] = mean(wprev*(pnorm(b,0,1) - pnorm(a,0,1)))/mean(wprev)
    }


    err  = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
    znew = AR*zprev + rt*err
    zhat = AR*znew

    wprev = wprev*(pnorm(b,0,1) - pnorm(a,0,1))

    if (xt[t]==0){
      preddist[,t-1] = c(0,temp[1])
    }else{
      preddist[,t-1] = cumsum(temp)[xt[t]:(xt[t]+1)]
    }
  }
  return(preddist)
}

PDvaluesARp_NoRes = function(theta, data, Regressor, ARMAorder, ParticleNumber,CountDist, epsilon){

  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments

  # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1:(1+nreg),
                           "Negative Binomial"   = 1:(2+nreg),
                           "Mixed Poisson"       = 1:(3+nreg),
                           "Generalized Poisson" = 1:(2+nreg),
                           "Binomial"            = 1:(2+nreg))

  if(nreg<1){
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = ppois,
                   "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
                   "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = pGpois,
                   "Binomial"            = pbinom
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = dGpois,
                   "Binomial"            = dbinom
    )
  }else{
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
    )
    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ dnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ dGpois  (x, ConstTheta, DynamTheta)}
    )
  }


  # retrieve marginal distribution parameters
  MargParms  = theta[MargParmIndices]
  nMargParms = length(MargParms)
  nparms     = length(theta)

  # check if the number of parameters matches the model setting
  if(nMargParms + sum(ARMAorder)!=nparms) stop('The length of theta does not match the model specification.')


  # Add the regressor to the parameters--only works for Negative Binomial
  if(CountDist == "Negative Binomial" && nreg>0){
    beta = MargParms[1:(nreg+1)]
    k = MargParms[nreg+2]
    m = exp(Regressor%*%beta)
    r = 1/k
    p = k*m/(1+k*m)
    ConstMargParm = r
    DynamMargParm = p
  }

  if(CountDist == "Generalized Poisson" && nreg>0){
    beta   = MargParms[1:(nreg+1)]
    ConstMargParm  = MargParms[nreg+2]
    DynamMargParm = exp(Regressor%*%beta)
  }


  # retrieve ARMA parameters
  AR = NULL
  if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]

  MA = NULL
  if(ARMAorder[2]>0) MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]


  #set.seed(1)
  xt = data
  T1 = length(xt)
  N  = ParticleNumber # number of particles
  preddist = matrix(0,2,T1-ARMAorder[1]) # to collect the values of predictive distribution of interest
  wgh = matrix(0,T1,N)        # to collect all particle weights

  # allocate memory for zprev
  ZprevAll = matrix(0,ARMAorder[1],N)

  # Compute integral limits
  if(nreg==0){
    a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
    b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)
  }else{
    a = rep( qnorm(mycdf(xt[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
    b = rep( qnorm(mycdf(xt[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
  }

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

  # save the currrent normal variables
  ZprevAll[1,] = zprev

  # initial estimate of first AR coefficient as Gamma(1)/Gamma(0) and corresponding error
  phit = TacvfAR(AR)[2]/TacvfAR(AR)[1]
  rt = as.numeric(sqrt(1-phit^2))

  # particle filter weights
  wprev = rep(1,N)
  wgh[1,] = wprev


  if (ARMAorder[1]>=2){
    for (t in 2:ARMAorder[1]){

      # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
      if (t==2) {
        zhat = ZprevAll[1:(t-1),]*phit
      } else{
        zhat = colSums(ZprevAll[1:(t-1),]*phit)
      }

      # compute random errors from truncated normal
      err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

      # compute the new Z and add it to the previous ones
      znew = rbind(zhat + rt*err, ZprevAll[1:(t-1),])
      ZprevAll[1:t,] = znew

      # recompute weights
      wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
      wprev = wgh[t,]

      # use YW equation to compute estimates of phi and of the errors
      Gt = toeplitz(TacvfAR(AR)[1:t])
      gt = TacvfAR(AR)[2:(t+1)]
      phit = as.numeric(solve(Gt) %*% gt)
      rt =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/TacvfAR(AR)[1]))

    }
  }



  for(t in (ARMAorder[1]+1):T1){

    # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
    if(ARMAorder[1]>1){# colsums doesnt work for 1-dimensional matrix
      zhat = colSums(ZprevAll*phit)
    }else{
      zhat=ZprevAll*phit
    }

    # compute a,b and temp
    temp = rep(0,(xt[t]+1))
    for (x in 0:xt[t]){
      # Compute integral limits
      if(nreg==0){
        a = as.numeric(qnorm(mycdf(x-1,MargParms),0,1) - zhat)/rt
        b = as.numeric(qnorm(mycdf(x,  MargParms),0,1) - zhat)/rt
      }else{
        a = as.numeric(qnorm(mycdf(x-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
        b = as.numeric(qnorm(mycdf(x,  ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
      }
      temp[x+1] = mean(pnorm(b,0,1) - pnorm(a,0,1))
    }

    # draw errors from truncated normal
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

    # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
    znew = zhat + rt*err


    # recompute weights
    wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
    wprev = wgh[t,]

    # save particles
    ZprevAll = rbind(znew, ZprevAll[1:(ARMAorder[1]-1),])

    # compute predictive distribution
    if (xt[t]==0){
      preddist[,t-ARMAorder[1]] = c(0,temp[1])
    }else{
      preddist[,t-ARMAorder[1]] = cumsum(temp)[xt[t]:(xt[t]+1)]
    }



  }
  return(preddist)
}



# compute predictive distribution AR(p)--new file resampling
PDvaluesARp = function(theta, data, Regressor, ARMAorder, ParticleNumber,CountDist, epsilon){

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments

  # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1:(1+nreg),
                           "Negative Binomial"   = 1:(2+nreg),
                           "Mixed Poisson"       = 1:(3+nreg),
                           "Generalized Poisson" = 1:(2+nreg),
                           "Binomial"            = 1:(2+nreg))

  if(nreg<1){
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = ppois,
                   "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
                   "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = pGpois,
                   "Binomial"            = pbinom
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = dGpois,
                   "Binomial"            = dbinom
    )
  }else{
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
    )
    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ dnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ dGpois  (x, ConstTheta, DynamTheta)}
    )
  }


  # retrieve marginal distribution parameters
  MargParms  = theta[MargParmIndices]
  nMargParms = length(MargParms)
  nparms     = length(theta)

  # check if the number of parameters matches the model setting
  if(nMargParms + sum(ARMAorder)!=nparms) stop('The length of theta does not match the model specification.')


  # Add the regressor to the parameters--only works for Negative Binomial
  if(CountDist == "Negative Binomial" && nreg>0){
    beta = MargParms[1:(nreg+1)]
    k = MargParms[nreg+2]
    m = exp(Regressor%*%beta)
    r = 1/k
    p = k*m/(1+k*m)
    ConstMargParm = r
    DynamMargParm = p
  }

  if(CountDist == "Generalized Poisson" && nreg>0){
    beta   = MargParms[1:(nreg+1)]
    ConstMargParm  = MargParms[nreg+2]
    DynamMargParm = exp(Regressor%*%beta)
  }


  # retrieve ARMA parameters
  AR = NULL
  if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]

  MA = NULL
  if(ARMAorder[2]>0) MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]


  #set.seed(1)
  xt = data
  T1 = length(xt)
  N  = ParticleNumber # number of particles
  preddist = matrix(0,2,T1-ARMAorder[1]) # to collect the values of predictive distribution of interest
  wgh = matrix(0,T1,N)        # to collect all particle weights

  # allocate memory for zprev
  ZprevAll = matrix(0,ARMAorder[1],N)

  # Compute integral limits
  if(nreg==0){
    a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
    b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)
  }else{
    a = rep( qnorm(mycdf(xt[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
    b = rep( qnorm(mycdf(xt[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
  }

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

  # save the currrent normal variables
  ZprevAll[1,] = zprev

  # initial estimate of first AR coefficient as Gamma(1)/Gamma(0) and corresponding error
  phit = TacvfAR(AR)[2]/TacvfAR(AR)[1]
  rt = as.numeric(sqrt(1-phit^2))

  # particle filter weights
  wprev = rep(1,N)
  wgh[1,] = wprev


  if (ARMAorder[1]>=2){
    for (t in 2:ARMAorder[1]){

      # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
      if (t==2) {
        zhat = ZprevAll[1:(t-1),]*phit
      } else{
        zhat = colSums(ZprevAll[1:(t-1),]*phit)
      }

      # compute random errors from truncated normal
      err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

      # compute the new Z and add it to the previous ones
      znew = rbind(zhat + rt*err, ZprevAll[1:(t-1),])
      ZprevAll[1:t,] = znew

      # recompute weights
      wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
      wprev = wgh[t,]

      # use YW equation to compute estimates of phi and of the errors
      Gt = toeplitz(TacvfAR(AR)[1:t])
      gt = TacvfAR(AR)[2:(t+1)]
      phit = as.numeric(solve(Gt) %*% gt)
      rt =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/TacvfAR(AR)[1]))

    }
  }



  for(t in (ARMAorder[1]+1):T1){

    # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
    if(ARMAorder[1]>1){# colsums doesnt work for 1-dimensional matrix
      zhat = colSums(ZprevAll*phit)
    }else{
      zhat=ZprevAll*phit
    }


    # compute a,b and temp
    temp = rep(0,(xt[t]+1))
    for (x in 0:xt[t]){
      # Compute integral limits
      if(nreg==0){
        a = as.numeric(qnorm(mycdf(x-1,MargParms),0,1) - zhat)/rt
        b = as.numeric(qnorm(mycdf(x,  MargParms),0,1) - zhat)/rt
      }else{
        a = as.numeric(qnorm(mycdf(x-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
        b = as.numeric(qnorm(mycdf(x,  ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
      }
      temp[x+1] = mean(pnorm(b,0,1) - pnorm(a,0,1))
    }

    # draw errors from truncated normal
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

    # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
    znew = zhat + rt*err

    # compute unnormalized weights
    wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)

    # normalized weights
    wghn = wgh[t,]/sum(wgh[t,])

    old_state1 <- get_rand_state()

    # sample indices from multinomial distribution-see Step 4 of SISR in paper
    ESS = 1/sum(wghn^2)

    if(ESS<epsilon*N){
      ind = rmultinom(1,N,wghn)
      # sample particles
      znew = rep(znew,ind)
    }
    set_rand_state(old_state1)

    # save particles
    if (ARMAorder[1]>1){
      ZprevAll = rbind(znew, ZprevAll[1:(ARMAorder[1]-1),])
    }else {
      ZprevAll[1,]=znew
    }

    if (xt[t]==0){
      preddist[,t-ARMAorder[1]] = c(0,temp[1])
    }else{
      preddist[,t-ARMAorder[1]] = cumsum(temp)[xt[t]:(xt[t]+1)]
    }



  }
  return(preddist)
}


# compute residuals
ComputeResiduals = function(theta, Regressor, ARMAorder, CountDist){
  #--------------------------------------------------------------------------#
  # PURPOSE:  Compute the residuals in relation (66)
  #
  # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
  #           details. A first version of the paper can be found at:
  #           https://arxiv.org/abs/1811.00203
  #
  # INPUTS:
  #
  # OUTPUT:
  #    loglik:
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    July 2020
  #--------------------------------------------------------------------------#



  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments

  # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1:(1+nreg),
                           "Negative Binomial"   = 1:(2+nreg),
                           "Mixed Poisson"       = 1:(3+nreg),
                           "Generalized Poisson" = 1:(2+nreg),
                           "Binomial"            = 1:(2+nreg))
  if(nreg<1){
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = ppois,
                   "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
                   "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = pGpois,
                   "Binomial"            = pbinom
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = dGpois,
                   "Binomial"            = dbinom
    )
  }else{
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = function(x, ConstTheta, DynamTheta){             ppois   (x, DynamTheta)},
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
    )
    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = function(x, ConstTheta, DynamTheta){             dpois   (x, DynamTheta)},
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ dnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ dGpois  (x, ConstTheta, DynamTheta)}
    )
  }



  # retrieve marginal distribution parameters
  MargParms  = theta[MargParmIndices]
  nMargParms = length(MargParms)
  nparms     = length(theta)

  # check if the number of parameters matches the model setting
  if(nMargParms + sum(ARMAorder)!=nparms) stop('The length of theta does not match the model specification.')


  # Add the regressor to the parameters--only works for Negative Binomial
  if(CountDist == "Negative Binomial" && nreg>0){
    beta = MargParms[1:(nreg+1)]
    k = MargParms[nreg+2]
    m = exp(Regressor%*%beta)
    r = 1/k
    p = k*m/(1+k*m)
    ConstMargParm = r
    DynamMargParm = p
  }

  if(CountDist == "Generalized Poisson" && nreg>0){
    beta   = MargParms[1:(nreg+1)]
    ConstMargParm  = MargParms[nreg+2]
    DynamMargParm = exp(Regressor%*%beta)
  }

  if(CountDist == "Poisson" && nreg>0){
    beta           = MargParms[1:(nreg+1)]
    ConstMargParm  = NULL
    DynamMargParm  = exp(Regressor%*%beta)
  }


  # retrieve ARMA parameters
  AR = NULL
  if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]

  MA = NULL
  if(ARMAorder[2]>0) MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]

  xt = data

  # sample size
  n = length(xt)

  # allocate memory
  Zhat <- rep(-99,n)

  # pGpois funcion doesnt allow for vetor lambda so a for-loop is needed, however the change shouldnt
  # be too hard. note the formula subtracts 1 from the counts Xt so it may lead to negative numbers
  # thats why there is and if else below.
  for (i in 1:n){
    k <- data[i]
    if (k != 0) {
      if(nreg==0){
        # Compute limits
        a = qnorm(mycdf(xt[i]-1,t(MargParms)),0,1)
        b = qnorm(mycdf(xt[i],t(MargParms)),0,1)
      }else{
        a = qnorm(mycdf(xt[i]-1, ConstMargParm, DynamMargParm[i]) ,0,1)
        b = qnorm(mycdf(xt[i], ConstMargParm, DynamMargParm[i]) ,0,1)
      }

      Zhat[i] <- (exp(-a^2/2)-exp(-b^2/2))/sqrt(2*pi)/(mycdf(k, ConstMargParm, DynamMargParm[i])-mycdf(k-1,ConstMargParm, DynamMargParm[i]))
    }else{
      if(nreg==0){
        # Compute integral limits
        b = qnorm(mycdf(xt[i],t(MargParms)),0,1)
      }else{
        b = qnorm(mycdf(xt[i], ConstMargParm, DynamMargParm[i]) ,0,1)
      }
      Zhat[i] <- -exp(-b^2/2)/sqrt(2*pi)/mycdf(0, ConstMargParm, DynamMargParm[i])
    }
  }

  # apply AR filter--Fix me allow for MA as well
  residual = data.frame(filter(Zhat,c(1,-AR))[1:(n-ARMAorder[1])])

  names(residual) = "residual"
  return(residual)
}





# Compute Predictive Distribution of AR1 Poisson--Old file
PDvaluesPoissonAR1Old = function(theta, phi, data,   ParticleNumber){
  cdf = function(x, theta){ ppois(x, lambda=theta[1]) }
  pdf = function(x, theta){ dpois(x, lambda=theta[1]) }


  #set.seed(1)
  xt = data
  T1 = length(xt)
  N = ParticleNumber# number of particles
  preddist = matrix(0,2,T1-1) # to collect the values of predictive distribution of interest

  a = qnorm(cdf(xt[1]-1,theta),0,1)
  b = qnorm(cdf(xt[1],theta),0,1)
  a = rep(a,N)
  b = rep(b,N)
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
  zhat = phi*zprev

  wprev = rep(1,N)

  for (t in 2:T1){
    temp = rep(0,(xt[t]+1))
    for (x in 0:xt[t]){
      rt = sqrt(1-phi^2)
      a = (qnorm(cdf(x-1,theta),0,1) - zhat)/rt
      b = (qnorm(cdf(x,theta),0,1) - zhat)/rt
      temp[x+1] = mean(wprev*(pnorm(b,0,1) - pnorm(a,0,1)))/mean(wprev)
    }
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
    znew = phi*zprev + rt*err
    zhat = phi*znew

    wprev = wprev*(pnorm(b,0,1) - pnorm(a,0,1))

    if (xt[t]==0){
      preddist[,t-1] = c(0,temp[1])
    }else{
      preddist[,t-1] = cumsum(temp)[xt[t]:(xt[t]+1)]
    }
  }
  return(preddist)
}


# Compute Predictive Distribution of MA2 Poisson with Regressors -- Old file
PDvaluesMA2OldReg = function(theta, tht, data, Regressor, ARMAorder, ParticleNumber){
  #set.seed(1)
  cdf = function(x, theta){ ppois(x, lambda=theta[1]) }
  pdf = function(x, theta){ dpois(x, lambda=theta[1]) }

  xt = data
  T1 = length(xt)

  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)
  theta0 = theta[1:(nreg+1)]
  theta1 = exp(Regressor%*%theta0)

  N = ParticleNumber # number of particles
  preddist = matrix(0,2,T1-1) # to collect the values of predictive distribution of interest


  gam1 = tht[1]*(1+tht[2])/(1+tht[1]^2+tht[2]^2)
  gam2 = tht[2]/(1+tht[1]^2+tht[2]^2)

  a = qnorm(cdf(xt[1]-1,theta1[1,]),0,1)
  b = qnorm(cdf(xt[1],theta1[1,]),0,1)
  a = rep(a,N)
  b = rep(b,N)
  zprev1 = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
  tht11 = gam1
  zhat1 = tht11*zprev1

  rt1 = ( 1 - tht11^2 )^(1/2)
  wprev = rep(1,N)

  temp = rep(0,(xt[2]+1))
  for (x in 0:xt[2]){
    a2 = (qnorm(cdf(x-1,theta1[2,]),0,1) - zhat1)/rt1
    b2 = (qnorm(cdf(x,theta1[2,]),0,1) - zhat1)/rt1
    temp[x+1] = mean(wprev*(pnorm(b2,0,1) - pnorm(a2,0,1)))/mean(wprev)
  }
  err2 = qnorm(runif(length(a2),0,1)*(pnorm(b2,0,1)-pnorm(a2,0,1))+pnorm(a2,0,1),0,1)
  znew2 = zhat1 + rt1*err2
  znew_all = cbind(znew2,zprev1)

  tht22 = gam2
  tht21 = (gam1 - tht11*tht22)/rt1^2
  tht_all = cbind(rep(tht21,N),rep(tht22,N))
  zhat_all = cbind(zhat1,rep(0,N))
  zhat2 = rowSums(tht_all*(znew_all-zhat_all))

  rt2 = (1 - tht21^2*rt1^2 - tht22^2)^(1/2)

  wprev = wprev*(pnorm(b2,0,1) - pnorm(a2,0,1))

  rt_all = c(rt2,rt1)
  zhat_all = cbind(zhat2,zhat1)

  if (xt[2]==0){
    preddist[,2-1] = c(0,temp[1])
  }else{
    preddist[,2-1] = cumsum(temp)[xt[2]:(xt[2]+1)]
  }

  for (t in 3:T1)
  {

    temp = rep(0,(xt[t]+1))
    for (x in 0:xt[t]){
      a = (qnorm(cdf(x-1,theta1[t,]),0,1) - zhat_all[,1])/rt_all[1]
      b = (qnorm(cdf(x,theta1[t,]),0,1) - zhat_all[,1])/rt_all[1]
      temp[x+1] = mean(wprev*(pnorm(b,0,1) - pnorm(a,0,1)))/mean(wprev)
    }
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
    znew = zhat_all[,1] + rt_all[1]*err
    znew_all = cbind(znew,znew_all[,1])

    thtt2 = gam2/rt_all[2]^2
    thtt1 = (gam1 - tht_all[1]*thtt2*rt_all[2]^2)/rt_all[1]^2
    tht_all = cbind(rep(thtt1,N),rep(thtt2,N))
    zhat2 = rowSums(tht_all*(znew_all-zhat_all))

    rt = (1 - thtt1^2*rt_all[1]^2 - thtt2^2*rt_all[2]^2)^(1/2)

    wprev = wprev*(pnorm(b,0,1) - pnorm(a,0,1))

    rt_all = c(rt,rt_all[1])
    zhat_all = cbind(zhat2,zhat_all[,1])

    if (xt[t]==0){
      preddist[,t-1] = c(0,temp[1])
    }else{
      preddist[,t-1] = cumsum(temp)[xt[t]:(xt[t]+1)]
    }

  }

  return(preddist)

}




# two functions to help me plot the qqplot
# nsim <- function(n, m = 0, s = 1) {
#   z <- rnorm(n)
#   m + s * ((z - mean(z)) / sd(z))
# }

nboot <- function(x, R) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  do.call(rbind,
          lapply(1 : R,
                 function(i) {
                   xx <- sort(nsim(n, m, s))
                   p <- seq_along(x) / n - 0.5 / n
                   data.frame(x = xx, p = p, sim = i)
                 }))
}







#---------NON RESAMPLING FUNCTIONS need to be updated---------#
FitMultiplePF = function(initialParam, data, CountDist, nfit, ParticleSchemes){
  # Let nparts by the length of the vector ParticleSchemes.
  # This function fits maximizes the PF likelihood, nfit manys times for nparts many choices
  # of particle numbers, thus yielding a total of nfit*nparts many estimates.

  # how many choices for the number of particles
  nparts = length(ParticleSchemes)
  nparms = length(initialParam)

  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst = matrix(0,nrow=nfit*nparts,ncol=nparms)
  se =  matrix(NA,nrow=nfit*nparts,ncol=nparms)
  loglik = rep(0,nfit*nparts)


  # Each realization will be fitted nfit*nparts many times
  for (j in 1:nfit){
    set.seed(j)
    # for each fit repeat for different number of particles
    for (k in 1:nparts){
      # number of particles to be used
      ParticleNumber = ParticleSchemes[k]

      # remove the ParticleNumber from the likelihood function arguments
      myfun = function(theta,data)LikSISGenDist_ARp_Res(theta, data, ParticleNumber, CountDist)

      # run optimization for our model
      optim.output <- optim(par = initialParam, fn = myfun,
                            data=data,
                            hessian=TRUE, method = "BFGS")

      # save estimates, loglik value and diagonal hessian
      ParmEst[nfit*(k-1)+j,]  = optim.output$par
      loglik[nfit*(k-1) +j]   = optim.output$value
      se[nfit*(k-1)+j,]       = sqrt(abs(diag(solve(optim.output$hessian))))
    }
  }

  All = cbind(ParmEst, se, loglik)
  return(All)
}


ParticleFilterMA1 = function(theta, data, ARMAorder, ParticleNumber, CountDist){
  #------------------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
  #           of the a specified count time series model with an underlying MA(1)
  #           dependence structure.
  #
  # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
  #           details. A first version of the paper can be found at:
  #           https://arxiv.org/abs/1811.00203
  #           2. This function is very similar to LikSISGenDist_ARp but here
  #           I have a resampling step.
  #
  # INPUTS:
  #    theta:            parameter vector
  #    data:             data
  #    ParticleNumber:   number of particles to be used.
  #    CountDist:        count marginal distribution
  #    epsilon           resampling when ESS<epsilon*N
  #
  # OUTPUT:
  #    loglik:           approximate log-likelihood
  #
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    April 2020
  #------------------------------------------------------------------------------------#

  #print(theta)
  # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments
  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))
  # retrieve indices of marginal distribution parameters
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1,
                           "Negative Binomial"   = 1:2,
                           "Mixed Poisson"       = 1:3,
                           "Generalized Poisson" = 1:2,
                           "Binomial"            = 1:2)

  # retrieve marginal cdf
  mycdf = switch(CountDist,
                 "Poisson"                       = ppois,
                 "Negative Binomial"             = function(x, theta){ pnbinom(q = x, size = theta[1], prob = 1-theta[2]) },
                 "Mixed Poisson"                 = pMixedPoisson,
                 "Generalized Poisson"           = pGenPoisson,
                 "Binomial"                      = pbinom
  )

  # retrieve marginal pdf
  mypdf = switch(CountDist,
                 "Poisson"                       = dpois,
                 "Negative Binomial"             = function(x, theta){ dnbinom(x, size = theta[1], prob = 1-theta[2]) },
                 "Mixed Poisson"                 = dMixedPoisson,
                 "Generalized Poisson"           = dGenPoisson,
                 "Binomial"                      = dbinom
  )

  # retrieve marginal distribution parameters
  MargParms = theta[MargParmIndices]
  nMargParms = length(MargParms) # num param in MargParms


  # retrieve MA parameters
  tht = theta[nMargParms + 1]
  xt = data
  T1 = length(xt)
  N = ParticleNumber          # number of particles
  wgh = matrix(0,T1,N)        # to collect all particle weights



  # Compute integral limits
  a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
  b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

  # Exact form for one-step ahead prediction in MA(1) ase
  rt0 = 1+tht^2
  zhat = tht*zprev/rt0

  # particle filter weights
  wprev = rep(1,N)
  wgh[1,] = wprev
  nloglik = 0 # initialize likelihood

  for (t in 2:T1){
    rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
    rt = sqrt(rt0/(1+tht^2))

    # compute limits of truncated normal distribution
    a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - zhat)/rt
    b = as.numeric(qnorm(mycdf(xt[t],MargParms),0,1) -   zhat)/rt

    # draw errors from truncated normal
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

    # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
    znew = zhat + rt*err

    zhat = tht*(znew-zhat)/rt0

    wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
    wprev = wgh[t,]
  }

  # likelihood approximation
  lik = mypdf(xt[1],MargParms)*mean(na.omit(wgh[T1,]))

  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  nloglik = -log(lik)
  #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  if (nloglik==Inf | is.na(nloglik)){
    out = 10^6
  }else{
    out = nloglik
  }

  return(out)
}

ParticleFilterMA1New = function(theta, data, ARMAorder, ParticleNumber, CountDist){
  #------------------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
  #           of the a specified count time series model with an underlying MA(1)
  #           dependence structure.
  #
  # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
  #           details. A first version of the paper can be found at:
  #           https://arxiv.org/abs/1811.00203
  #           2. This function is very similar to LikSISGenDist_ARp but here
  #           I have a resampling step.
  #
  # INPUTS:
  #    theta:            parameter vector
  #    data:             data
  #    ParticleNumber:   number of particles to be used.
  #    CountDist:        count marginal distribution
  #    epsilon           resampling when ESS<epsilon*N
  #
  # OUTPUT:
  #    loglik:           approximate log-likelihood
  #
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    April 2020
  #------------------------------------------------------------------------------------#

  #print(theta)
  # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments

  # retrieve indices of marginal distribution parameters
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1,
                           "Negative Binomial"   = 1:2,
                           "Mixed Poisson"       = 1:3,
                           "Generalized Poisson" = 1:2,
                           "Binomial"            = 1:2)

  # retrieve marginal cdf
  mycdf = switch(CountDist,
                 "Poisson"                       = ppois,
                 "Negative Binomial"             = function(x, theta){ pnbinom(q = x, size = theta[1], prob = 1-theta[2]) },
                 "Mixed Poisson"                 = pMixedPoisson,
                 "Generalized Poisson"           = pGenPoisson,
                 "Binomial"                      = pbinom
  )

  # retrieve marginal pdf
  mypdf = switch(CountDist,
                 "Poisson"                       = dpois,
                 "Negative Binomial"             = function(x, theta){ dnbinom(x, size = theta[1], prob = 1-theta[2]) },
                 "Mixed Poisson"                 = dMixedPoisson,
                 "Generalized Poisson"           = dGenPoisson,
                 "Binomial"                      = dbinom
  )

  # retrieve marginal distribution parameters
  MargParms = theta[MargParmIndices]
  nMargParms = length(MargParms) # num param in MargParms


  # retrieve MA parameters
  tht = theta[nMargParms + 1]
  xt = data
  T1 = length(xt)
  N = ParticleNumber          # number of particles
  wgh = matrix(0,T1,N)        # to collect all particle weights



  # Compute integral limits
  a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
  b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

  # Exact form for one-step ahead prediction in MA(1) ase
  rt0 = 1+tht^2
  zhat = tht*zprev/rt0

  # particle filter weights
  wprev = rep(1,N)
  wgh[1,] = wprev
  nloglik = 0 # initialize likelihood

  for (t in 2:T1){
    rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
    rt = sqrt(rt0/(1+tht^2))

    # compute limits of truncated normal distribution
    a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - zhat)/rt
    b = as.numeric(qnorm(mycdf(xt[t],MargParms),0,1) -   zhat)/rt

    # draw errors from truncated normal
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

    # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
    znew = zhat + rt*err

    zhat = tht*(znew-zhat)/rt0

    wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
    wprev = wgh[t,]
  }



  # break if I got NA weight
  if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
    nloglik = 10^6
    return(nloglik)
  }

  # likelihood approximation
  lik = mypdf(xt[1],MargParms)*mean(na.omit(wgh[T1,]))

  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  nloglik = -log(lik)
  #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  out =nloglik


  return(out)
}

#---------new wrapper to fit PF likelihood---------#
FitMultiplePFMA1New = function(x0, X, CountDist, Particles, LB, UB, ARMAorder, epsilon, UseDEOptim){
  #====================================================================================#
  # PURPOSE       Fit the Particle Filter log-likelihood. This function maximizes
  #               the PF likelihood, nfit manys times for nparts many choices of
  #               particle numbers, thus yielding a total of nfit*nparts many estimates
  #
  # INPUT
  #   x0          initial parameters
  #   X           count series
  #   CountDist   prescribed count distribution
  #   Particles   vector with different choices for number of particles
  #   LB          parameter lower bounds
  #   UB          parameter upper bounds
  #   ARMAorder   order of the udnerlying ARMA model
  #   epsilon     resampling when ESS<epsilon*N
  #   UseDEOptim  flag, if=1 then use Deoptim for optimization
  #
  #   MaxCdf      cdf will be computed up to this number (for light tails cdf=1 fast)
  #   nHC         number of HC to be computed
  #
  # OUTPUT
  #   All         parameter estimates, standard errors, likelihood value
  #
  # NOTES         I may comment out se in cases where maximum is achieved at the boundary
  #
  # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date          April 2020
  # Version       3.6.3
  #====================================================================================#

  # how many choices for the number of particles
  nparts = length(Particles)
  nparms = length(x0)
  nfit = 1
  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst = matrix(0,nrow=nfit*nparts,ncol=nparms)
  se =  matrix(NA,nrow=nfit*nparts,ncol=nparms)
  loglik = rep(0,nfit*nparts)

  n = length(X)

  # Each realization will be fitted nfit*nparts many times
  for (j in 1:nfit){
    set.seed(j)
    # for each fit repeat for different number of particles
    for (k in 1:nparts){
      # number of particles to be used
      ParticleNumber = Particles[k]
      if(!UseDEOptim){
        if (n<400){
          # run optimization for our model
          optim.output <- optim(par            = x0,
                                fn             = ParticleFilterMA1,
                                data           = X,
                                ARMAorder      = ARMAorder,
                                ParticleNumber = ParticleNumber,
                                CountDist      = CountDist,
                                lower          = LB,
                                upper          = UB,
                                hessian        = TRUE,
                                method         = "L-BFGS-B")

        }else{
          # run optimization for our model
          optim.output <- optim(par            = x0,
                                fn             = likSISRMA1,
                                data           = X,
                                ARMAorder      = ARMAorder,
                                ParticleNumber = ParticleNumber,
                                CountDist      = CountDist,
                                epsilon        = epsilon,
                                lower          = LB,
                                upper          = UB,
                                hessian        = TRUE,
                                method         = "L-BFGS-B")
          }
        }else{
          if(n<400){
            optim.output<- DEoptim::DEoptim(fn             = likSISRMA1,
                                            lower          = LB,
                                            upper          = UB,
                                            data           = X,
                                            ARMAorder      = ARMAorder,
                                            ParticleNumber = ParticleNumber,
                                            CountDist      = CountDist,
                                            control        = DEoptim::DEoptim.control(trace = 10, itermax = 200, steptol = 50, reltol = 1e-5))


          }else{
            optim.output<- DEoptim::DEoptim(fn             = ParticleFilterMA1_Res,
                                            lower          = LB,
                                            upper          = UB,
                                            data           = X,
                                            ARMAorder      = ARMAorder,
                                            ParticleNumber = ParticleNumber,
                                            CountDist      = CountDist,
                                            epsilon        = epsilon,
                                            control        = DEoptim::DEoptim.control(trace = 10, itermax = 200, steptol = 50, reltol = 1e-5))
          }
        }
      if(UseDEOptim){
        ParmEst[nfit*(k-1)+j,]  = optim.output$optim$bestmem
        loglik[nfit*(k-1) +j]   = optim.output$optim$bestval
      }
      else{
        # save estimates, loglik value and diagonal hessian
        ParmEst[nfit*(k-1)+j,]  = optim.output$par
        loglik[nfit*(k-1) +j]   = optim.output$value
        se[nfit*(k-1)+j,]       = sqrt(abs(diag(solve(optim.output$hessian))))
      }
    }
  }

  All = cbind(ParmEst, se, loglik)
  return(All)
}







# #---------new wrapper to fit PF likelihood---------#
# FitMultiplePFRes = function(x0, X, CountDist, Particles, LB, UB, ARMAorder, epsilon){
#   #====================================================================================#
#   # PURPOSE       Fit the Particle Filter log-likelihood. This function maximizes
#   #               the PF likelihood, nfit manys times for nparts many choices of
#   #               particle numbers, thus yielding a total of nfit*nparts many estimates
#   #
#   # INPUT
#   #   x0          initial parameters
#   #   X           count series
#   #   CountDist   prescribed count distribution
#   #   Particles   vector with different choices for number of particles
#   #   LB          parameter lower bounds
#   #   UB          parameter upper bounds
#   #   ARMAorder   order of the udnerlying ARMA model
#   #   epsilon     resampling when ESS<epsilon*N
#   #   UseDEOptim  flag, if=1 then use Deoptim for optimization
#   #
#   #   MaxCdf      cdf will be computed up to this number (for light tails cdf=1 fast)
#   #   nHC         number of HC to be computed
#   #
#   # OUTPUT
#   #   All         parameter estimates, standard errors, likelihood value
#   #
#   # NOTES         I may comment out se in cases where maximum is achieved at the boundary
#   #
#   # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date          April 2020
#   # Version       3.6.3
#   #====================================================================================#
#
#   # how many choices for the number of particles
#   nparts = length(Particles)
#   nparms = length(x0)
#   nfit = 1
#   # allocate memory to save parameter estimates, hessian values, and loglik values
#   ParmEst = matrix(0,nrow=nfit*nparts,ncol=nparms)
#   se =  matrix(NA,nrow=nfit*nparts,ncol=nparms)
#   loglik = rep(0,nfit*nparts)
#   convcode = rep(0,nfit*nparts)
#   kkt1 = rep(0,nfit*nparts)
#   kkt2 = rep(0,nfit*nparts)
#
#   n = length(X)
#
#   # Each realization will be fitted nfit*nparts many times
#   for (j in 1:nfit){
#     set.seed(j)
#     # for each fit repeat for different number of particles
#     for (k in 1:nparts){
#       # number of particles to be used
#       ParticleNumber = Particles[k]
#       # run optimization for our model
#       optim.output <- optimx(par           = x0,
#                              fn             = ParticleFilterRes,
#                              data           = X,
#                              ARMAorder      = ARMAorder,
#                              ParticleNumber = ParticleNumber,
#                              CountDist      = CountDist,
#                              epsilon        = epsilon,
#                              lower          = LB,
#                              upper          = UB,
#                              hessian        = TRUE,
#                              method         = "L-BFGS-B")
#
#       # save estimates, loglik value and diagonal hessian
#       ParmEst[nfit*(k-1)+j,]  = c(optim.output$p1,optim.output$p2,optim.output$p3,optim.output$p4)
#       loglik[nfit*(k-1) +j]   = optim.output$value
#       convcode[nfit*(k-1) +j] = optim.output$convcode
#       kkt1[nfit*(k-1) +j]     = optim.output$kkt1
#       kkt2[nfit*(k-1) +j]     = optim.output$kkt2
#
#       # compute hessian
#       H = gHgen(par            = ParmEst[nfit*(k-1)+j,],
#                 fn             = ParticleFilterRes,
#                 data           = X,
#                 ARMAorder      = ARMAorder,
#                 CountDist      = CountDist,
#                 ParticleNumber = ParticleNumber,
#                 epsilon        = epsilon
#       )
#
#
#       # save standard errors from Hessian
#       if(H$hessOK && det(H$Hn)>10^(-8)){
#         se[nfit*(k-1)+j,]   = sqrt(abs(diag(solve(H$Hn))))
#       }else{
#         se[nfit*(k-1)+j,] = rep(NA, nparms)
#       }
#
#
#     }
#   }
#
#   All = cbind(ParmEst, se, loglik, convcode, kkt1, kkt2)
#   return(All)
# }
#
#
# #---------new wrapper to fit PF likelihood---------#
# FitMultiplePFResReg = function(x0, X, Regressor, CountDist, Particles, LB, UB, ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod){
#   #====================================================================================#
#   # PURPOSE       Fit the Particle Filter log-likelihood. This function maximizes
#   #               the PF likelihood, nfit manys times for nparts many choices of
#   #               particle numbers, thus yielding a total of nfit*nparts many estimates.
#   #               Here I am also considering a regressor
#   #
#   # INPUT
#   #   x0          initial parameters
#   #   X           count series
#   #   CountDist   prescribed count distribution
#   #   Particles   vector with different choices for number of particles
#   #   LB          parameter lower bounds
#   #   UB          parameter upper bounds
#   #   ARMAorder   order of the udnerlying ARMA model
#   #   epsilon     resampling when ESS<epsilon*N
#   #   UseDEOptim  flag, if=1 then use Deoptim for optimization
#   #
#   #   MaxCdf      cdf will be computed up to this number (for light tails cdf=1 fast)
#   #   nHC         number of HC to be computed
#   #
#   # OUTPUT
#   #   All         parameter estimates, standard errors, likelihood value
#   #
#   # NOTES         I may comment out se in cases where maximum is achieved at the boundary
#   #
#   # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date          April 2020
#   # Version       3.6.3
#   #====================================================================================#
#
#   # how many choices for the number of particles
#   nparts = length(Particles)
#   nparms = length(x0)
#   nfit = 1
#   # allocate memory to save parameter estimates, hessian values, and loglik values
#   ParmEst = matrix(0,nrow=nfit*nparts,ncol=nparms)
#   se =  matrix(NA,nrow=nfit*nparts,ncol=nparms)
#   loglik = rep(0,nfit*nparts)
#   convcode = rep(0,nfit*nparts)
#   kkt1 = rep(0,nfit*nparts)
#   kkt2 = rep(0,nfit*nparts)
#
#   n = length(X)
#
#   # Each realization will be fitted nfit*nparts many times
#   for (j in 1:nfit){
#     set.seed(j)
#     # for each fit repeat for different number of particles
#     for (k in 1:nparts){
#       # number of particles to be used
#       ParticleNumber = Particles[k]
#       # run optimization for our model --no ARMA model allowed
#       if(ARMAorder[1]>0){
#         optim.output <- optimx(par           = x0,
#                                fn             = ParticleFilterAR_Res,
#                                data           = X,
#                                Regressor      = Regressor,
#                                ARMAorder      = ARMAorder,
#                                ParticleNumber = ParticleNumber,
#                                CountDist      = CountDist,
#                                epsilon        = epsilon,
#                                lower          = LB,
#                                upper          = UB,
#                                hessian        = TRUE,
#                                method         = OptMethod)
#       }else{
#         optim.output <- optimx(par            = x0,
#                                fn             = ParticleFilterMA_Res,
#                                data           = X,
#                                Regressor      = Regressor,
#                                ARMAorder      = ARMAorder,
#                                ParticleNumber = ParticleNumber,
#                                CountDist      = CountDist,
#                                epsilon        = epsilon,
#                                lower          = LB,
#                                upper          = UB,
#                                hessian        = TRUE,
#                                method         = OptMethod)
#
#       }
#
#       # save estimates, loglik value and diagonal hessian
#       ParmEst[nfit*(k-1)+j,]  = as.numeric(optim.output[1:nparms])
#       loglik[nfit*(k-1) +j]   = optim.output$value
#       convcode[nfit*(k-1) +j] = optim.output$convcode
#       kkt1[nfit*(k-1) +j]     = optim.output$kkt1
#       kkt2[nfit*(k-1) +j]     = optim.output$kkt2
#
#       if(ARMAorder[1]>0){
#         # compute hessian
#         H = gHgen(par          = ParmEst[nfit*(k-1)+j,],
#                   fn             = ParticleFilterAR_Res,
#                   data           = X,
#                   Regressor      = Regressor,
#                   ARMAorder      = ARMAorder,
#                   CountDist      = CountDist,
#                   ParticleNumber = ParticleNumber,
#                   epsilon        = epsilon
#         )
#       }else{
#         H = gHgen(par            = ParmEst[nfit*(k-1)+j,],
#                   fn             = ParticleFilterMA_Res,
#                   data           = X,
#                   Regressor      = Regressor,
#                   ARMAorder      = ARMAorder,
#                   CountDist      = CountDist,
#                   ParticleNumber = ParticleNumber,
#                   epsilon        = epsilon
#         )
#       }
#
#       if(H$hessOK && det(H$Hn)>10^(-8)){
#         se[nfit*(k-1)+j,]   = sqrt(abs(diag(solve(H$Hn))))
#       }else{
#         se[nfit*(k-1)+j,] = rep(NA, nparms)
#       }
#
#     }
#   }
#   # Compute model selection criteria
#   Criteria = ComputeCriteria(loglik, nparms, n, Particles)
#
#
#   # get the names of the final output
#   parmnames = colnames(optim.output)
#   mynames = c(parmnames[1:nparms],paste("se", parmnames[1:nparms], sep="_"), "loglik", "AIC", "BIC","AICc", "status", "kkt1", "kkt2")
#
#
#   All = matrix(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)
#   colnames(All) = mynames
#
#   return(All)
# }
#
# #---------new wrapper to fit PF likelihood---------#
# FitMultiplePFMA1Res = function(x0, X, CountDist, Particles, LB, UB, ARMAorder, epsilon){
#   #====================================================================================#
#   # PURPOSE       Fit the Particle Filter log-likelihood. This function maximizes
#   #               the PF likelihood, nfit manys times for nparts many choices of
#   #               particle numbers, thus yielding a total of nfit*nparts many estimates
#   #
#   # INPUT
#   #   x0          initial parameters
#   #   X           count series
#   #   CountDist   prescribed count distribution
#   #   Particles   vector with different choices for number of particles
#   #   LB          parameter lower bounds
#   #   UB          parameter upper bounds
#   #   ARMAorder   order of the udnerlying ARMA model
#   #   epsilon     resampling when ESS<epsilon*N
#   #   UseDEOptim  flag, if=1 then use Deoptim for optimization
#   #
#   #   MaxCdf      cdf will be computed up to this number (for light tails cdf=1 fast)
#   #   nHC         number of HC to be computed
#   #
#   # OUTPUT
#   #   All         parameter estimates, standard errors, likelihood value
#   #
#   # NOTES         I may comment out se in cases where maximum is achieved at the boundary
#   #
#   # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date          April 2020
#   # Version       3.6.3
#   #====================================================================================#
#
#   # how many choices for the number of particles
#   nparts = length(Particles)
#   nparms = length(x0)
#   nfit = 1
#   # allocate memory to save parameter estimates, hessian values, and loglik values
#   ParmEst = matrix(0,nrow=nfit*nparts,ncol=nparms)
#   se =  matrix(NA,nrow=nfit*nparts,ncol=nparms)
#   loglik = rep(0,nfit*nparts)
#   convcode = rep(0,nfit*nparts)
#   kkt1 = rep(0,nfit*nparts)
#   kkt2 = rep(0,nfit*nparts)
#
#   n = length(X)
#
#   # Each realization will be fitted nfit*nparts many times
#   for (j in 1:nfit){
#     set.seed(j)
#     # for each fit repeat for different number of particles
#     for (k in 1:nparts){
#       # number of particles to be used
#       ParticleNumber = Particles[k]
#       # run optimization for our model
#
#       optim.output <- optimx(par           = x0,
#                              fn             = ParticleFilterMA1_Res,
#                              data           = X,
#                              ARMAorder      = ARMAorder,
#                              ParticleNumber = ParticleNumber,
#                              CountDist      = CountDist,
#                              epsilon        = epsilon,
#                              lower          = LB,
#                              upper          = UB,
#                              hessian        = TRUE,
#                              method         = "L-BFGS-B")
#
#       # save estimates, loglik value and diagonal hessian
#       ParmEst[nfit*(k-1)+j,]  = c(optim.output$p1,optim.output$p2,optim.output$p3)
#       loglik[nfit*(k-1) +j]   = optim.output$value
#       convcode[nfit*(k-1) +j] = optim.output$convcode
#       kkt1[nfit*(k-1) +j]     = optim.output$kkt1
#       kkt2[nfit*(k-1) +j]     = optim.output$kkt2
#
#       # compute hessian
#       H = gHgen(par            = ParmEst[nfit*(k-1)+j,],
#                 fn             = ParticleFilterMA1_Res,
#                 data           = X,
#                 ARMAorder      = ARMAorder,
#                 CountDist      = CountDist,
#                 ParticleNumber = ParticleNumber,
#                 epsilon        = epsilon
#       )
#
#       if(H$hessOK){
#         se[nfit*(k-1)+j,]   =  sqrt(abs(diag(solve(H$Hn))))
#       }
#     }
#   }
#
#   All = cbind(ParmEst, se, loglik, convcode, kkt1, kkt2)
#   return(All)
# }








# PF likelihood with resampling
# ParticleFilterRes = function(theta, data, ARMAorder, ParticleNumber, CountDist, epsilon){
#   #--------------------------------------------------------------------------#
#   # REDUNANT--Initially I used this function for AR(p) ParticleFilter resampling
#   #           without regressors. The function ParticleFilterAR_Res  considers
#   #           general AR(p) with and without regressors.
#   # PURPOSE:  Use particle filtering with resampling
#   #           to approximate the likelihood of the
#   #           a specified count time series model with an underlying AR(p)
#   #           dependence structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paer can be found at:
#   #           https://arxiv.org/abs/1811.00203
#   #           2. This function is very similar to LikSISGenDist_ARp but here
#   #           I have a resampling step.
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             data
#   #    ParticleNumber:   number of particles to be used.
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    November 2019
#   #--------------------------------------------------------------------------#
#
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments
#
#   # retrieve indices of marginal distribution parameters
#   MargParmIndices = switch(CountDist,
#                            "Poisson"             = 1,
#                            "Negative Binomial"   = 1:2,
#                            "Mixed Poisson"       = 1:3,
#                            "Generalized Poisson" = 1:2,
#                            "Binomial"            = 1:2)
#
#   # retrieve marginal cdf
#   mycdf = switch(CountDist,
#                  "Poisson"             = ppois,
#                  "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
#                  "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
#                  "Generalized Poisson" = pGenPoisson,
#                  "Binomial"            = pbinom
#   )
#
#   # retrieve marginal pdf
#   mypdf = switch(CountDist,
#                  "Poisson"             = dpois,
#                  "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
#                  "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
#                  "Generalized Poisson" = dGenPoisson,
#                  "Binomial"            = dbinom
#   )
#
#   # retrieve marginal distribution parameters
#   MargParms  = theta[MargParmIndices]
#   nMargParms = length(MargParms)
#   nparms     = length(theta)
#
#
#   # retrieve ARMA parameters
#   if(ARMAorder[1]>0){
#     AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
#   }else{
#     AR = NULL
#   }
#
#   if(ARMAorder[2]>0){
#     MA = theta[ (length(theta) - ARMAorder[2]) : length(theta)]
#   }else{
#     MA = NULL
#   }
#
#   # retrieve AR order
#   ARorder = ARMAorder[1]
#
#   if (prod(abs(polyroot(c(1,-AR))) > 1)){ # check if the ar model is causal
#
#     xt = data
#     T1 = length(xt)
#     N = ParticleNumber          # number of particles
#     prt = matrix(0,N,T1)        # to collect all particles
#     wgh = matrix(0,T1,N)        # to collect all particle weights
#
#     # allocate memory for zprev
#     ZprevAll = matrix(0,ARMAorder[1],N)
#
#     # Compute integral limits
#     a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
#     b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)
#
#     # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#     zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#     # save the currrent normal variables
#     ZprevAll[1,] = zprev
#
#     # initial estimate of first AR coefficient as Gamma(1)/Gamma(0) and corresponding error
#     phit = TacvfAR(AR)[2]/TacvfAR(AR)[1]
#     rt = as.numeric(sqrt(1-phit^2))
#
#     # particle filter weights
#     wprev = rep(1,N)
#     wgh[1,] = wprev
#     nloglik = 0 # initialize likelihood
#     #t0 = proc.time()
#     # First p steps:
#
#     if (ARorder>=2){
#       for (t in 2:ARorder){
#
#         # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
#         if (t==2) {
#           ZpreviousTimesPhi = ZprevAll[1:(t-1),]*phit
#         } else{
#           ZpreviousTimesPhi = colSums(ZprevAll[1:(t-1),]*phit)
#         }
#
#         # Recompute integral limits
#         a = (qnorm(mycdf(xt[t]-1,MargParms),0,1) - ZpreviousTimesPhi)/rt
#         b = (qnorm(mycdf(xt[t],MargParms),0,1) - ZpreviousTimesPhi)/rt
#
#         # compute random errors from truncated normal
#         err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#         # compute the new Z and add it to the previous ones
#         znew = rbind(ZpreviousTimesPhi + rt*err, ZprevAll[1:(t-1),])
#         ZprevAll[1:t,] = znew
#
#         # recompute weights
#         wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
#         wprev = wgh[t,]
#
#         # use YW equation to compute estimates of phi and of the erros
#         Gt = toeplitz(TacvfAR(AR)[1:t])
#         gt = TacvfAR(AR)[2:(t+1)]
#         phit = as.numeric(solve(Gt) %*% gt)
#         rt =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/TacvfAR(AR)[1]))
#
#       }
#     }
#
#
#     # From p to T1 I dont need to estimate phi anymore
#     for (t in (ARorder+1):T1){
#       # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
#       if(ARorder>1){# colsums doesnt work for 1-dimensional matrix
#         ZpreviousTimesPhi = colSums(ZprevAll*AR)
#       }else{
#         ZpreviousTimesPhi=ZprevAll*AR
#       }
#       # compute limits of truncated normal distribution
#       a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - ZpreviousTimesPhi)/rt
#       b = as.numeric(qnorm(mycdf(xt[t],MargParms),0,1) -   ZpreviousTimesPhi)/rt
#
#       # draw errors from truncated normal
#       err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#       # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
#       znew = ZpreviousTimesPhi + rt*err
#
#       # Resampling Step--here the function differs from LikSISGenDist_ARp
#
#       # compute unnormalized weights
#       wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)
#
#       # break if I got NA weight
#       if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
#         nloglik = 10^8
#         break
#       }
#
#       # normalized weights
#       wghn = wgh[t,]/sum(wgh[t,])
#
#       old_state1 <- get_rand_state()
#
#       # sample indices from multinomial distribution-see Step 4 of SISR in paper
#       ESS = sum(1/wghn^2)
#       if(ESS<epsilon*N){
#         ind = rmultinom(1,N,wghn)
#         # sample particles
#         znew = rep(znew,ind)
#
#         # use low variance resampling
#         #znew = lowVarianceRS(znew, wghn, N)
#       }
#       set_rand_state(old_state1)
#
#
#       # save particles
#       if (ARorder>1){
#         ZprevAll = rbind(znew, ZprevAll[1:(ARorder-1),])
#       }else {
#         ZprevAll[1,]=znew
#       }
#       # update likelihood
#       nloglik = nloglik - log(mean(wgh[t,]))
#     }
#
#     # likelihood approximation
#     nloglik = nloglik - log(mypdf(xt[1],MargParms))
#
#
#     # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#     nloglik = nloglik
#     #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#
#     out =nloglik
#
#   }else{
#     nloglik = 10^8
#   }
#
#   return(out)
# }


#
# # PF likelihood with resampling for MA(1)
# ParticleFilterMA1_Res = function(theta, data, ARMAorder, ParticleNumber, CountDist, epsilon){
#   #------------------------------------------------------------------------------------#
#   # REDUNANT--Initially I used this function for MA(1) ParticleFilter resampling
#   #           without regressors. The function ParticleFilterMA_Res_Reg  considers
#   #           general MA(q) with and without regressors.
#   #
#   # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
#   #           of the a specified count time series model with an underlying MA(1)
#   #           dependence structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paper can be found at:
#   #           https://arxiv.org/abs/1811.00203
#   #           2. This function is very similar to LikSISGenDist_ARp but here
#   #           I have a resampling step.
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             data
#   #    ParticleNumber:   number of particles to be used.
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    April 2020
#   #------------------------------------------------------------------------------------#
#
# old_state <- get_rand_state()
# on.exit(set_rand_state(old_state))
#  #print(theta)
#   # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments
#
#   # retrieve indices of marginal distribution parameters
#   MargParmIndices = switch(CountDist,
#                            "Poisson"             = 1,
#                            "Negative Binomial"   = 1:2,
#                            "Mixed Poisson"       = 1:3,
#                            "Generalized Poisson" = 1:2,
#                            "Binomial"            = 1:2)
#
#   # retrieve marginal cdf
#   mycdf = switch(CountDist,
#                  "Poisson"                       = ppois,
#                  "Negative Binomial"             = function(x, theta){ pnbinom(x,  theta[1], 1-theta[2]) },
#                  "Mixed Poisson"                 = pMixedPoisson,
#                  "Generalized Poisson"           = pGenPoisson,
#                  "Binomial"                      = pbinom
#   )
#
#   # retrieve marginal pdf
#   mypdf = switch(CountDist,
#                  "Poisson"                       = dpois,
#                  "Negative Binomial"             = function(x, theta){ dnbinom(x,  theta[1],  1-theta[2]) },
#                  "Mixed Poisson"                 = dMixedPoisson,
#                  "Generalized Poisson"           = dGenPoisson,
#                  "Binomial"                      = dbinom
#   )
#
#   # retrieve marginal distribution parameters
#   MargParms = theta[MargParmIndices]
#   nMargParms = length(MargParms) # num param in MargParms
#
#
#   # retrieve MA parameters
#     tht = theta[nMargParms + 1]
#     xt = data
#     T1 = length(xt)
#     N = ParticleNumber          # number of particles
#     wgh = matrix(0,T1,N)        # to collect all particle weights
#
#
#     # Compute integral limits
#     a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
#     b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)
#
#     # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#     zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#     # Exact form for one-step ahead prediction in MA(1) ase
#     rt0 = 1+tht^2
#     zhat = tht*zprev/rt0
#
#     # particle filter weights
#     wprev = rep(1,N)
#     wgh[1,] = wprev
#     nloglik = 0 # initialize likelihood
#
#     for (t in 2:T1){
#       rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
#       rt = sqrt(rt0/(1+tht^2))
#
#       # compute limits of truncated normal distribution
#       a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - zhat)/rt
#       b = as.numeric(qnorm(mycdf(xt[t],MargParms),0,1) -   zhat)/rt
#
#       # draw errors from truncated normal
#       err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#       # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
#       znew = zhat + rt*err
#
#       # compute unnormalized weights
#       wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)
#
#       # break if I got NA weight
#       if (any(is.na(wgh[t,]))| sum(wgh[t,])<10^(-8) ){
#         nloglik = 10^8
#         break
#       }
#
#       # normalized weights
#       wghn = wgh[t,]/sum(wgh[t,])
#
#       old_state1 <- get_rand_state()
#
#       # Resampling: sample indices from multinomial distribution-see Step 4 of SISR in paper
#       ESS = 1/sum(wghn^2)
#
#       if(ESS<epsilon*N){
#         ind = rmultinom(1,N,wghn)
#         # sample particles
#         znew = rep(znew,ind)
#
#         # use low variance resampling
#         #znew = lowVarianceRS(znew, wghn, N)
#       }
#       zhat = tht*(znew-zhat)/rt0
#       set_rand_state(old_state1)
#
#       # update likelihood
#       nloglik = nloglik - log(mean(wgh[t,]))
#     }
#
#     # likelihood approximation
#     nloglik = nloglik - log(mypdf(xt[1],MargParms))
#
#
#     # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#     nloglik = nloglik
#     #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#
#     out =nloglik
#
#     if (out==Inf | is.na(out)){
#       out = 10^8
#     }
#
#   return(out)
# }



# PF likelihood with resampling for MA(1) with regressors
# ParticleFilterMA1_Res_Reg = function(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon){
#   #------------------------------------------------------------------------------------#
#   # REDUNANT--Initially I used this functyion for MA(1) ParticleFilter resampling
#   #           with regressors. The function ParticleFilterMA_Res_Reg  considers
#   #           general MA(q).
#   # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
#   #           of the a specified count time series model with an underlying MA(1)
#   #           dependence structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paper can be found at:
#   #           https://arxiv.org/abs/1811.00203
#   #           2. This function is very similar to LikSISGenDist_ARp but here
#   #           I have a resampling step.
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             data
#   #    ParticleNumber:   number of particles to be used.
#   #    Regressor:        independent variable
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    July 2020
#   #------------------------------------------------------------------------------------#
#
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # number of regressors assuming there is an intercept
#   nreg = dim(Regressor)[2]-1
#   # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments
#
#   # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
#   MargParmIndices = switch(CountDist,
#                            "Poisson"             = 1:(1+nreg),
#                            "Negative Binomial"   = 1:(2+nreg),
#                            "Mixed Poisson"       = 1:(3+nreg),
#                            "Generalized Poisson" = 1:(2+nreg),
#                            "Binomial"            = 1:(2+nreg))
#   if(nreg<1){
#     # retrieve marginal cdf
#     mycdf = switch(CountDist,
#                    "Poisson"             = ppois,
#                    "Negative Binomial"   = function(x, theta){ pnbinom (q = x, size = theta[1], prob = 1-theta[2])},
#                    "Mixed Poisson"       = function(x, theta){ pmixpois(x, p = theta[1], lam1 = theta[2], lam2 = theta[3])},
#                    "Generalized Poisson" = pGpois,
#                    "Binomial"            = pbinom
#     )
#
#     # retrieve marginal pdf
#     mypdf = switch(CountDist,
#                    "Poisson"             = dpois,
#                    "Negative Binomial"   = function(x, theta){ dnbinom (x, size = theta[1], prob = 1-theta[2]) },
#                    "Mixed Poisson"       = function(x, theta){ dmixpois(x, p = theta[1], lam1 = theta[2], lam2 = theta[3])},
#                    "Generalized Poisson" = dGpois,
#                    "Binomial"            = dbinom
#     )
#   }else{
#     # retrieve marginal cdf
#     mycdf = switch(CountDist,
#                    "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
#                    "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
#     )
#     # retrieve marginal pdf
#     mypdf = switch(CountDist,
#                    "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ dnbinom (x, ConstTheta, 1-DynamTheta)},
#                    "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ dGpois  (x, ConstTheta, DynamTheta)}
#     )
#   }
#
#
#
#   # retrieve marginal distribution parameters
#   MargParms  = theta[MargParmIndices]
#   nMargParms = length(MargParms)
#   nparms     = length(theta)
#
#
#   # Add the regressor to the parameters--only works for Negative Binomial
#   if(CountDist == "Negative Binomial"){
#     beta = MargParms[1:(nreg+1)]
#     k = MargParms[nreg+2]
#     m = exp(Regressor%*%beta)
#     r = 1/k
#     p = k*m/(1+k*m)
#     ConstMargParm = r
#     DynamMargParm = p
#   }
#
#   if(CountDist == "Generalized Poisson"){
#     beta   = MargParms[1:(nreg+1)]
#     ConstMargParm  = MargParms[nreg+2]
#     DynamMargParm = exp(Regressor%*%beta)
#   }
#
#
#   # retrieve ARMA parameters
#   AR = NULL
#   if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]
#
#   MA = NULL
#   if(ARMAorder[2]>0) MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]
#
#   if (checkPoly(AR,MA)[1]=="Causal" && checkPoly(AR,MA)[2]=="Invertible"){
#     # retrieve MA parameters
#
#     tht = MA[1]
#     xt = data
#     T1 = length(xt)
#     N = ParticleNumber          # number of particles
#     wgh = matrix(0,T1,N)        # to collect all particle weights
#
#     if(nreg==0){
#       # Compute integral limits
#       a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
#       b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)
#     }else{
#       a = rep( qnorm(mycdf(xt[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
#       b = rep( qnorm(mycdf(xt[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
#     }
#
#     # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#     zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#     # Exact form for one-step ahead prediction in MA(1) case
#     if (is.null(MA) && is.null(AR)){
#       rt0 = 1
#       zhat = 0
#     }else{
#       rt0 = 1+tht^2
#       zhat = tht*zprev/rt0
#     }
#
#
#     # particle filter weights
#     wprev = rep(1,N)
#     wgh[1,] = wprev
#     nloglik = 0 # initialize likelihood
#     for (t in 2:T1){
#
#       if (is.null(MA) && is.null(AR)){
#         rt0 = 1
#         rt  = 1
#       }else{
#         rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
#         rt = sqrt(rt0/(1+tht^2))
#       }
#
#       if(nreg==0){
#         # compute limits of truncated normal distribution
#         a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - zhat)/rt
#         b = as.numeric(qnorm(mycdf(xt[t],MargParms),0,1) -   zhat)/rt
#       }else{
#         a = as.numeric(qnorm(mycdf(xt[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#         b = as.numeric(qnorm(mycdf(xt[t],ConstMargParm, DynamMargParm[t]),0,1) -   zhat)/rt
#       }
#
#       # draw errors from truncated normal
#       err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#       # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
#       znew = zhat + rt*err
#
#       # compute unnormalized weights
#       wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)
#
#       # break if I got NA weight
#       if (any(is.na(wgh[t,]))| sum(wgh[t,])<10^(-8) ){
#         nloglik = 10^8
#         break
#       }
#
#       # normalized weights
#       wghn = wgh[t,]/sum(wgh[t,])
#
#       old_state1 <- get_rand_state()
#
#       # Resampling: sample indices from multinomial distribution-see Step 4 of SISR in paper
#       ESS = 1/sum(wghn^2)
#
#       if(ESS<epsilon*N){
#         ind = rmultinom(1,N,wghn)
#         # sample particles
#         znew = rep(znew,ind)
#
#         # use low variance resampling
#         #znew = lowVarianceRS(znew, wghn, N)
#       }
#
#       if (is.null(MA) && is.null(AR)){
#         zhat = 0
#       }else{
#         zhat = tht*(znew-zhat)/rt0
#       }
#       set_rand_state(old_state1)
#
#       # update likelihood
#       nloglik = nloglik - log(mean(wgh[t,]))
#     }
#
#     # likelihood approximation
#     if(nreg<1){
#       nloglik = nloglik - log(mypdf(xt[1],MargParms))
#     }else{
#       nloglik = nloglik - log(mypdf(xt[1], ConstMargParm, DynamMargParm[1]))
#       }
#
#     # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#     nloglik = nloglik
#     #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#
#     out =nloglik
#
#     if (out==Inf | is.na(out)){
#       out = 10^8
#       }
#     }else{
#       out = 10^9
#       }
#
#   return(out)
# }




# z.rest = function(a,b){
#   # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
#   qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
# }
#
# likSISR = function(theta, data){
#   cdf = function(x, theta1){ pnbinom(q = x, size = theta1[1], prob = theta1[2]) }
#   pdf = function(x, theta1){ dnbinom(x = x, size = theta1[1], prob = theta1[2]) }
#   #set.seed(1)
#   theta1.idx = 1:2
#   theta2.idx = 3
#   theta1 = theta[theta1.idx]
#   n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
#   theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
#   tht = theta[theta2.idx]
#   xt = data
#   T1 = length(xt)
#   N = 1000 # number of particles
#   prt = matrix(0,N,T1) # to collect all particles
#   wgh = matrix(0,N,T1) # to collect all particle weights
#
#   a = qnorm(cdf(xt[1]-1,theta1),0,1)
#   b = qnorm(cdf(xt[1],theta1),0,1)
#   a = rep(a,N)
#   b = rep(b,N)
#   zprev = z.rest(a,b)
#   rt0 = 1+tht^2
#   zhat = tht*zprev/rt0
#   prt[,1] = zhat
#
#
#   nloglik <- 0
#   for (t in 2:T1)
#   {
#     rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
#     rt = sqrt(rt0/(1+tht^2))
#     a = (qnorm(cdf(xt[t]-1,theta1),0,1) - zhat)/rt
#     b = (qnorm(cdf(xt[t],theta1),0,1) - zhat)/rt
#     err = z.rest(a,b)
#     znew =   + rt*err
#     wgh <- pnorm(b,0,1) - pnorm(a,0,1)
#     if (any(is.na(wgh))) # see apf code below for explanation
#     {
#       #nloglik <- NaN
#       nloglik <- Inf
#       break
#     }
#     if (sum(wgh)==0)
#     {
#       #nloglik <- NaN
#       nloglik <- Inf
#       break
#     }
#     wghn <- wgh/sum(wgh)
#     ind <- rmultinom(1, N, wghn)
#     znew <- rep(znew,ind)
#
#     zhat = tht*(znew-zhat)/rt0
#     prt[,t] = zhat
#
#     nloglik <- nloglik -2*log(mean(wgh))
#   }
#
#   nloglik <- nloglik - 2*log(pdf(xt[1],theta1))
#
#
#   out = if (is.na(nloglik)) Inf else nloglik
#   return(out)
# }




# Add some functions that I will need in particle filtering approximation of
# likelihood. See file LikSIS_ARpGenDist.R


# # Add some functions that I will need
# get_rand_state <- function() {
#   # Using `get0()` here to have `NULL` output in case object doesn't exist.
#   # Also using `inherits = FALSE` to get value exactly from global environment
#   # and not from one of its parent.
#   get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
# }
#
# set_rand_state <- function(state) {
#   # Assigning `NULL` state might lead to unwanted consequences
#   if (!is.null(state)) {
#     assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
#   }
# }
#
#
# # innovations algorithm code
# innovations.algorithm <- function(acvf,n.max=length(acvf)-1){
#   # Found this onlinbe need to check it
#   # http://faculty.washington.edu/dbp/s519/R-code/innovations-algorithm.R
#   thetas <- vector(mode="list",length=n.max)
#   vs <- rep(acvf[1],n.max+1)
#   for(n in 1:n.max)
#   {
#     thetas[[n]] <- rep(0,n)
#     thetas[[n]][n] <- acvf[n+1]/vs[1]
#     if(n>1)
#     {
#       for(k in 1:(n-1))
#       {
#         js <- 0:(k-1)
#         thetas[[n]][n-k] <- (acvf[n-k+1] - sum(thetas[[k]][k-js]*thetas[[n]][n-js]*vs[js+1]))/vs[k+1]
#       }
#     }
#     js <- 0:(n-1)
#     vs[n+1] <- vs[n+1] - sum(thetas[[n]][n-js]^2*vs[js+1])
#   }
#   structure(list(vs=vs,thetas=thetas))
# }


# PF likelihood with resampling for AR(p)
# ParticleFilter_Res_AR = function(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon){
#   #--------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling
#   #           to approximate the likelihood of the
#   #           a specified count time series model with an underlying AR(p)
#   #           dependence structure. A singloe dummy regression is added here.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paer can be found at:
#   #           https://arxiv.org/abs/1811.00203
#   #           2. This function is very similar to LikSISGenDist_ARp but here
#   #           I have a resampling step.
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             data
#   #    ParticleNumber:   number of particles to be used.
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    July  2020
#   #--------------------------------------------------------------------------#
#
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # number of regressors assuming there is an intercept
#   nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)
#
#   # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments
#
#   # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
#   MargParmIndices = switch(CountDist,
#                            "Poisson"             = 1:(1+nreg),
#                            "Negative Binomial"   = 1:(2+nreg),
#                            "Mixed Poisson"       = 1:(3+nreg),
#                            "Generalized Poisson" = 1:(2+nreg),
#                            "Binomial"            = 1:(2+nreg))
#   if(nreg<1){
#     # retrieve marginal cdf
#     mycdf = switch(CountDist,
#                  "Poisson"             = ppois,
#                  "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
#                  "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
#                  "Generalized Poisson" = pGpois,
#                  "Binomial"            = pbinom
#     )
#
#     # retrieve marginal pdf
#     mypdf = switch(CountDist,
#                  "Poisson"             = dpois,
#                  "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
#                  "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
#                  "Generalized Poisson" = dGpois,
#                  "Binomial"            = dbinom
#     )
#     }else{
#       # retrieve marginal cdf
#       mycdf = switch(CountDist,
#                      "Poisson"             = function(x, ConstTheta, DynamTheta){             ppois   (x, DynamTheta)},
#                      "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
#                      "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
#       )
#       # retrieve marginal pdf
#       mypdf = switch(CountDist,
#                      "Poisson"             = function(x, ConstTheta, DynamTheta){             dpois   (x, DynamTheta)},
#                      "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ dnbinom (x, ConstTheta, 1-DynamTheta)},
#                      "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ dGpois  (x, ConstTheta, DynamTheta)}
#       )
#     }
#
#
#
#   # retrieve marginal distribution parameters
#   MargParms  = theta[MargParmIndices]
#   nMargParms = length(MargParms)
#   nparms     = length(theta)
#
#   # check if the number of parameters matches the model setting
#   if(nMargParms + sum(ARMAorder)!=nparms) stop('The length of theta does not match the model specification.')
#
#
#   # Add the regressor to the parameters--only works for Negative Binomial
#   if(CountDist == "Negative Binomial" && nreg>0){
#     beta = MargParms[1:(nreg+1)]
#     k = MargParms[nreg+2]
#     m = exp(Regressor%*%beta)
#     r = 1/k
#     p = k*m/(1+k*m)
#     ConstMargParm = r
#     DynamMargParm = p
#   }
#
#   if(CountDist == "Generalized Poisson" && nreg>0){
#     beta   = MargParms[1:(nreg+1)]
#     ConstMargParm  = MargParms[nreg+2]
#     DynamMargParm = exp(Regressor%*%beta)
#   }
#
#   if(CountDist == "Poisson" && nreg>0){
#     beta           = MargParms[1:(nreg+1)]
#     ConstMargParm  = NULL
#     DynamMargParm  = exp(Regressor%*%beta)
#   }
#
#
#   # retrieve ARMA parameters
#   AR = NULL
#   if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]
#
#   MA = NULL
#   if(ARMAorder[2]>0) MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]
#
#
#   if (checkPoly(AR,MA)[1]=="Causal" && checkPoly(AR,MA)[2]=="Invertible"){
#
#     xt = data
#     T1 = length(xt)
#     N = ParticleNumber          # number of particles
#     prt = matrix(0,N,T1)        # to collect all particles
#     wgh = matrix(0,T1,N)        # to collect all particle weights
#
#     # allocate memory for zprev
#     ZprevAll = matrix(0,ARMAorder[1],N)
#
#     if(nreg==0){
#       # Compute integral limits
#       a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
#       b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)
#     }else{
#       a = rep( qnorm(mycdf(xt[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
#       b = rep( qnorm(mycdf(xt[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
#     }
#     # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#     zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#     # save the currrent normal variables
#     ZprevAll[1,] = zprev
#
#     # initial estimate of first AR coefficient as Gamma(1)/Gamma(0) and corresponding error
#     phit = TacvfAR(AR)[2]/TacvfAR(AR)[1]
#     rt = as.numeric(sqrt(1-phit^2))
#
#     # particle filter weights
#     wprev = rep(1,N)
#     wgh[1,] = wprev
#     nloglik = 0 # initialize likelihood
#
#     # First p steps:
#     if (ARMAorder[1]>=2){
#       for (t in 2:ARMAorder[1]){
#
#         # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
#         if (t==2) {
#           zhat = ZprevAll[1:(t-1),]*phit
#         } else{
#           zhat = colSums(ZprevAll[1:(t-1),]*phit)
#         }
#
#         # Recompute integral limits
#         if(nreg==0){
#           a = (qnorm(mycdf(xt[t]-1,t(MargParms)),0,1) - zhat)/rt
#           b = (qnorm(mycdf(xt[t],t(MargParms)),0,1) - zhat)/rt
#         }else{
#           a = (qnorm(mycdf(xt[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#           b = (qnorm(mycdf(xt[t],ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#         }
#
#         # compute random errors from truncated normal
#         err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#         # compute the new Z and add it to the previous ones
#         znew = rbind(zhat + rt*err, ZprevAll[1:(t-1),])
#         ZprevAll[1:t,] = znew
#
#         # recompute weights
#         wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
#         wprev = wgh[t,]
#
#         # use YW equation to compute estimates of phi and of the erros
#         Gt = toeplitz(TacvfAR(AR)[1:t])
#         gt = TacvfAR(AR)[2:(t+1)]
#         phit = as.numeric(solve(Gt) %*% gt)
#         rt =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/TacvfAR(AR)[1]))
#
#       }
#     }
#
#
#     # From p to T1 I dont need to estimate phi anymore
#     for (t in (ARMAorder[1]+1):T1){
#
#       # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
#       if(ARMAorder[1]>1){# colsums doesnt work for 1-dimensional matrix
#         zhat = colSums(ZprevAll*AR)
#       }else{
#         zhat=ZprevAll*AR
#       }
#
#       # compute limits of truncated normal distribution
#       if(nreg==0){
#         a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - zhat)/rt
#         b = as.numeric(qnorm(mycdf(xt[t],  MargParms),0,1) - zhat)/rt
#       }else{
#         a = as.numeric(qnorm(mycdf(xt[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#         b = as.numeric(qnorm(mycdf(xt[t],  ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#       }
#
#       # draw errors from truncated normal
#       err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#       # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
#       znew = zhat + rt*err
#
#       # Resampling Step--here the function differs from LikSISGenDist_ARp
#
#       # compute unnormalized weights
#       wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)
#
#       # break if I got NA weight
#       if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
#         nloglik = 10^8
#         break
#       }
#
#       # normalized weights
#       wghn = wgh[t,]/sum(wgh[t,])
#
#       old_state1 <- get_rand_state()
#
#       # sample indices from multinomial distribution-see Step 4 of SISR in paper
#       ESS = 1/sum(wghn^2)
#       if(ESS<epsilon*N){
#         ind = rmultinom(1,N,wghn)
#         # sample particles
#         znew = rep(znew,ind)
#
#         # use low variance resampling
#         #znew = lowVarianceRS(znew, wghn, N)
#       }
#       set_rand_state(old_state1)
#
#
#       # save particles
#       if (ARMAorder[1]>1){
#         ZprevAll = rbind(znew, ZprevAll[1:(ARMAorder[1]-1),])
#       }else {
#         ZprevAll[1,]=znew
#       }
#       # update likelihood
#       nloglik = nloglik - log(mean(wgh[t,]))
#     }
#
#     # likelihood approximation
#     if(nreg==0){
#       nloglik = nloglik - log(mypdf(xt[1],MargParms))
#     }else{
#       nloglik = nloglik - log(mypdf(xt[1], ConstMargParm, DynamMargParm[1]))
#     }
#
#
#
#     # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#     nloglik = nloglik
#     #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#
#     out = nloglik
#
#   }else{
#     out = 10^9
#   }
#
#   return(out)
# }


# # PF likelihood with resampling for MA(q)
# ParticleFilter_Res_MA = function(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon){
#   #------------------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
#   #           of the a specified count time series model with an underlying MA(1)
#   #           dependence structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paper can be found at:
#   #           https://arxiv.org/abs/1811.00203
#   #           2. This function is very similar to LikSISGenDist_ARp but here
#   #           I have a resampling step.
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             data
#   #    ParticleNumber:   number of particles to be used.
#   #    Regressor:        independent variable
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    July 2020
#   #------------------------------------------------------------------------------------#
#
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # number of regressors assuming there is an intercept
#   nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)
#
#   # sample size and number of particles
#   T1  = length(data)
#   N   = ParticleNumber
#
#   # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments
#
#   # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
#   MargParmIndices = switch(CountDist,
#                            "Poisson"             = 1:(1+nreg),
#                            "Negative Binomial"   = 1:(2+nreg),
#                            "Mixed Poisson"       = 1:(3+nreg),
#                            "Generalized Poisson" = 1:(2+nreg),
#                            "Binomial"            = 1:(2+nreg))
#   if(nreg<1){
#     # retrieve marginal cdf
#     mycdf = switch(CountDist,
#                    "Poisson"             = ppois,
#                    "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
#                    "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
#                    "Generalized Poisson" = pGpois,
#                    "Binomial"            = pbinom
#     )
#
#     # retrieve marginal pdf
#     mypdf = switch(CountDist,
#                    "Poisson"             = dpois,
#                    "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
#                    "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
#                    "Generalized Poisson" = dGpois,
#                    "Binomial"            = dbinom
#     )
#   }else{
#     # retrieve marginal cdf
#     mycdf = switch(CountDist,
#                    "Poisson"             = function(x, ConstTheta, DynamTheta){             ppois   (x, DynamTheta)},
#                    "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
#                    "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
#     )
#     # retrieve marginal pdf
#     mypdf = switch(CountDist,
#                    "Poisson"             = function(x, ConstTheta, DynamTheta){             dpois   (x, DynamTheta)},
#                    "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ dnbinom (x, ConstTheta, 1-DynamTheta)},
#                    "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ dGpois  (x, ConstTheta, DynamTheta)}
#     )
#   }
#
#
#
#   # retrieve marginal distribution parameters
#   MargParms  = theta[MargParmIndices]
#   nMargParms = length(MargParms)
#   nparms     = length(theta)
#
#   # check if the number of parameters matches the model setting
#   if(nMargParms + sum(ARMAorder)!=nparms) stop('The length of theta does not match the model specification.')
#
#   # Add the regressor to the parameters--only works for Negative Binomial
#   if(CountDist == "Negative Binomial" && nreg>0){
#     beta          = MargParms[1:(nreg+1)]
#     k             = MargParms[nreg+2]
#     m             = exp(Regressor%*%beta)
#     r             = 1/k
#     p             = k*m/(1+k*m)
#     ConstMargParm = r
#     DynamMargParm = p
#   }
#
#   if(CountDist == "Generalized Poisson" && nreg>0){
#     beta           = MargParms[1:(nreg+1)]
#     ConstMargParm  = MargParms[nreg+2]
#     DynamMargParm  = exp(Regressor%*%beta)
#   }
#
#   if(CountDist == "Poisson" && nreg>0){
#     beta           = MargParms[1:(nreg+1)]
#     ConstMargParm  = NULL
#     DynamMargParm  = exp(Regressor%*%beta)
#   }
#
#
#   # retrieve AR parameters
#   AR = NULL
#   if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]
#
#   # retrieve MA parameters
#   MA = NULL
#   if(ARMAorder[2]>0) MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]
#
#   #--------------------focus on stable models----------------------------
#   if (checkPoly(AR,MA)[1]=="Causal" && checkPoly(AR,MA)[2]=="Invertible"){
#
#     # allocate matrix to collect all particle weights
#     wgh = matrix(0,length(data),N)
#
#     # Compute integral limits
#     if(nreg==0){
#       a = rep( qnorm(mycdf(data[1]-1,t(MargParms)),0,1), N)
#       b = rep( qnorm(mycdf(data[1],t(MargParms)),0,1), N)
#     }else{
#       a = rep( qnorm(mycdf(data[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
#       b = rep( qnorm(mycdf(data[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
#     }
#
#     # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#     zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#
#     # run innovations Algorithm for MA models that are not WN
#     if(ARMAorder[2]>0) Inn = matrix(0,N,ARMAorder[2])       # I will save here the q many innovations (Z - Zhat) --see (5.3.9) BD book
#     if (is.null(MA) && is.null(AR)){
#       v0   = 1
#       zhat = 0
#     }else{
#       MA.acvf <- as.vector(tacvfARMA(theta = MA, maxLag=T1))
#       ia = innovations.algorithm(MA.acvf)
#       Theta = ia$thetas
#       # first stage of Innovations
#       v0    = ia$vs[1]
#       zhat = -Theta[[1]][1]*zprev
#       Inn[,ARMAorder[2]] = zprev
#     }
#
#     # particle filter weights
#     wprev   = rep(1,N)
#     wgh[1,] = wprev
#     nloglik = 0 # initialize likelihood
#
#     for (t in 2:T1){
#
#       # update innovations quantities if not White noise
#       if (is.null(MA) && is.null(AR)){
#         vt=1
#       }else{
#         vt0 = ia$vs[t]
#         vt  = sqrt(vt0/v0)
#         }
#
#       # roll the old Innovations to earlier columns
#       if(ARMAorder[2]>1) Inn[,1:(ARMAorder[2]-1)] = Inn[,2:(ARMAorder[2])]
#
#       # compute limits of truncated normal distribution
#       if(nreg==0){
#         a = as.numeric(qnorm(mycdf(data[t]-1,MargParms),0,1) - zhat)/vt
#         b = as.numeric(qnorm(mycdf(data[t],MargParms),0,1) -   zhat)/vt
#       }else{
#         a = as.numeric(qnorm(mycdf(data[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/vt
#         b = as.numeric(qnorm(mycdf(data[t],ConstMargParm, DynamMargParm[t]),0,1) -   zhat)/vt
#       }
#
#       # draw errors from truncated normal
#       err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#       # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
#       znew = zhat + vt*err
#
#       # compute new innovation
#       Inn[,ARMAorder[2]] = (znew-zhat)
#
#       # compute unnormalized weights
#       wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)
#
#       # break if I got NA weight
#       if (any(is.na(wgh[t,]))| sum(wgh[t,])<10^(-8) ){
#         nloglik = 10^8
#         break
#       }
#
#       # normalized weights
#       wghn = wgh[t,]/sum(wgh[t,])
#
#       old_state1 <- get_rand_state()
#
#       # Resampling: sample indices from multinomial distribution-see Step 4 of SISR in paper
#       ESS = 1/sum(wghn^2)
#
#       if(ESS<epsilon*N){
#         ind = rmultinom(1,N,wghn)
#         # sample particles
#         znew = rep(znew,ind)
#
#         # use low variance resampling
#         #znew = lowVarianceRS(znew, wghn, N)
#       }
#
#       # update zhat--fix me can probably be vectorized
#       if (is.null(MA) && is.null(AR)){
#         zhat = 0
#       }else{
#         S = 0
#         for(j in 1:min(t,ARMAorder[2])){
#           S = S-Theta[[t]][j]*Inn[,ARMAorder[2]-j+1]
#         }
#         zhat = S
#       }
#
#       set_rand_state(old_state1)
#
#       # update likelihood
#       nloglik = nloglik - log(mean(wgh[t,]))
#     }
#
#     # likelihood approximation
#     if(nreg<1){
#       nloglik = nloglik - log(mypdf(data[1],MargParms))
#     }else{
#       nloglik = nloglik - log(mypdf(data[1], ConstMargParm, DynamMargParm[1]))
#     }
#
#     # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#     nloglik = nloglik
#     #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#
#     out =nloglik
#
#     if (out==Inf | is.na(out)){
#       out = 10^8
#     }
#   }else{
#     out = 10^9
#   }
#
#   return(out)
# }


# # PF likelihood with resampling
# ParticleFilter_Res = function(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon){
#   #--------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling
#   #           to approximate the likelihood of the
#   #           a specified count time series model with an underlying AR(p)
#   #           dependence structure or MA(q) structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paer can be found at:
#   #           https://arxiv.org/abs/1811.00203
#
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             dependent variable
#   #    Regressor:        independent variables
#   #    ParticleNumber:   number of particles to be used in likelihood approximation
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    July 2020
#   #--------------------------------------------------------------------------#
#
#   # check Particle number
#   if(epsilon > 1 || epsilon<0) stop('Please select a value between 0 and 1 for epsilon.')
#
#   # check Particle number
#   if(ParticleNumber<1) stop('Please select a nonegative value for the argument ParticleNumber.')
#
#   # check distributions
#   if ( !(CountDist %in%  c("Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial" )))
#     stop('The argument CountDist must take one of the following values:
#          "Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial".')
#
#   # check ARMAorder
#   if(prod(ARMAorder)<0 || length(ARMAorder)!= 2) stop('The argument ARMAorder must have length 2 and can not take negative values.')
#
#   # Mixed ARMA model
#   if(ARMAorder[1]>0 && ARMAorder[2]>0) stop('Please specify a pure AR or a pure MA model. ARMA(p,q) models with p>0 and q>0 have not yet been implemented.')
#
#   # Pure AR model
#   if(ARMAorder[1]>0 && ARMAorder[2]==0) loglik = ParticleFilter_Res_AR(theta, data, Regressor, ARMAorder,
#                                                                       ParticleNumber, CountDist, epsilon)
#   # Pure MA model or White noise
#   if(ARMAorder[1]==0&& ARMAorder[2]>=0) loglik = ParticleFilter_Res_MA(theta, data, Regressor, ARMAorder,
#                                                                       ParticleNumber, CountDist, epsilon)
#
#   return(loglik)
# }



# # Optimization wrapper to fit PF likelihood with resamplinbg
# FitMultiplePF_Res = function(theta, data, Regressor, ARMAorder, Particles, CountDist, epsilon, LB, UB, OptMethod){
#   #====================================================================================#
#   # PURPOSE       Fit the Particle Filter log-likelihood with resampling.
#   #               This function maximizes the PF likelihood, nfit manys times for nparts
#   #               many choices of particle numbers, thus yielding a total of nfit*nparts
#   #               many estimates.
#   #
#   # INPUT
#   #   theta       initial parameters
#   #   data        count series
#   #   CountDist   prescribed count distribution
#   #   Particles   vector with different choices for number of particles
#   #   ARMAorder   order of the udnerlying ARMA model
#   #   epsilon     resampling when ESS<epsilon*N
#   #   LB          parameter lower bounds
#   #   UB          parameter upper bounds
#   #   OptMethod
#   #
#   #
#   # OUTPUT
#   #   All         parameter estimates, standard errors, likelihood value, AIC, etc
#   #
#   # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date          April 2020
#   # Version       3.6.3
#   #====================================================================================#
#
#   # retrieve parameter, sample size etc
#   nparts = length(Particles)
#   nparms = length(theta)
#   nfit   = 1
#   n      = length(data)
#
#   # allocate memory to save parameter estimates, hessian values, and loglik values
#   ParmEst  = matrix(0,nrow=nfit*nparts,ncol=nparms)
#   se       =  matrix(NA,nrow=nfit*nparts,ncol=nparms)
#   loglik   = rep(0,nfit*nparts)
#   convcode = rep(0,nfit*nparts)
#   kkt1     = rep(0,nfit*nparts)
#   kkt2     = rep(0,nfit*nparts)
#
#
#   # Each realization will be fitted nfit*nparts many times
#   for (j in 1:nfit){
#     set.seed(j)
#     # for each fit repeat for different number of particles
#     for (k in 1:nparts){
#       # number of particles to be used
#       ParticleNumber = Particles[k]
#
#       # run optimization for our model --no ARMA model allowed
#       optim.output <- optimx(par            = theta,
#                              fn             = ParticleFilter_Res,
#                              data           = data,
#                              Regressor      = Regressor,
#                              ARMAorder      = ARMAorder,
#                              ParticleNumber = ParticleNumber,
#                              CountDist      = CountDist,
#                              epsilon        = epsilon,
#                              lower          = LB,
#                              upper          = UB,
#                              hessian        = TRUE,
#                              method         = OptMethod)
#
#
#
#       # save estimates, loglik value and diagonal hessian
#       ParmEst[nfit*(k-1)+j,]  = as.numeric(optim.output[1:nparms])
#       loglik[nfit*(k-1) +j]   = optim.output$value
#       convcode[nfit*(k-1) +j] = optim.output$convcode
#       kkt1[nfit*(k-1) +j]     = optim.output$kkt1
#       kkt2[nfit*(k-1) +j]     = optim.output$kkt2
#
#
#       # compute Hessian
#       H = gHgen(par            = ParmEst[nfit*(k-1)+j,],
#                 fn             = ParticleFilter_Res,
#                 data           = data,
#                 Regressor      = Regressor,
#                 ARMAorder      = ARMAorder,
#                 CountDist      = CountDist,
#                 ParticleNumber = ParticleNumber,
#                 epsilon        = epsilon)
#
#       # save standard errors from Hessian
#       if(H$hessOK && det(H$Hn)>10^(-8)){
#         se[nfit*(k-1)+j,]   = sqrt(abs(diag(solve(H$Hn))))
#       }else{
#         se[nfit*(k-1)+j,] = rep(NA, nparms)
#       }
#
#     }
#   }
#
#   # Compute model selection criteria (assuming one fit)
#   Criteria = ComputeCriteria(loglik, nparms, n, Particles)
#
#
#   # get the names of the final output
#   parmnames = colnames(optim.output)
#   mynames = c(parmnames[1:nparms],paste("se", parmnames[1:nparms], sep="_"), "loglik", "AIC", "BIC","AICc", "status", "kkt1", "kkt2")
#
#
#   All = matrix(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)
#   colnames(All) = mynames
#
#   return(All)
# }



