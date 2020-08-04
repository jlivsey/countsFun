# particle filter for MA
ParticleFilterRes_MA = function(xt,  Regressor, MargParms, ConstMargParm, DynamMargParm,
                                ARMAorder, AR, MA, ParticleNumber, mycdf, epsilon){


    old_state <- get_rand_state()
    on.exit(set_rand_state(old_state))

    # number of regressors assuming there is an intercept
    nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

    # sample size and number of particles
    T1  = length(xt)
    N   = ParticleNumber


    # allocate matrix to collect all particle weights
    wgh = matrix(0,length(xt),N)

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


    # run innovations Algorithm for MA models that are not WN
    if(ARMAorder[2]>0) Inn = matrix(0,N,ARMAorder[2])       # I will save here the q many innovations (Z - Zhat) --see (5.3.9) BD book
    if (is.null(MA) && is.null(AR)){
      v0   = 1
      zhat = 0
    }else{
      MA.acvf <- as.vector(tacvfARMA(theta = MA, maxLag=T1))
      ia = innovations.algorithm(MA.acvf)
      Theta = ia$thetas
      # first stage of Innovations
      v0    = ia$vs[1]
      zhat = -Theta[[1]][1]*zprev
      Inn[,ARMAorder[2]] = zprev
    }

    # particle filter weights
    wprev   = rep(1,N)
    wgh[1,] = wprev
    nloglik = 0 # initialize likelihood

    for (t in 2:T1){

      # update innovations quantities if not White noise
      if (is.null(MA) && is.null(AR)){
        vt=1
      }else{
        vt0 = ia$vs[t]
        vt  = sqrt(vt0/v0)
      }

      # roll the old Innovations to earlier columns
      if(ARMAorder[2]>1) Inn[,1:(ARMAorder[2]-1)] = Inn[,2:(ARMAorder[2])]

      # compute limits of truncated normal distribution
      if(nreg==0){
        a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - zhat)/vt
        b = as.numeric(qnorm(mycdf(xt[t],MargParms),0,1) -   zhat)/vt
      }else{
        a = as.numeric(qnorm(mycdf(xt[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/vt
        b = as.numeric(qnorm(mycdf(xt[t],ConstMargParm, DynamMargParm[t]),0,1) -   zhat)/vt
      }

      # draw errors from truncated normal
      err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

      # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
      znew = zhat + vt*err

      # compute new innovation
      Inn[,ARMAorder[2]] = (znew-zhat)

      # compute unnormalized weights
      wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)

      # break if I got NA weight
      if (any(is.na(wgh[t,]))| sum(wgh[t,])<10^(-8) ){
        nloglik = 10^8
        break
      }

      # normalized weights
      wghn = wgh[t,]/sum(wgh[t,])

      old_state1 <- get_rand_state()

      # Resampling: sample indices from multinomial distribution-see Step 4 of SISR in paper
      ESS = 1/sum(wghn^2)

      if(ESS<epsilon*N){
        ind = rmultinom(1,N,wghn)
        # sample particles
        znew = rep(znew,ind)

        # use low variance resampling
        #znew = lowVarianceRS(znew, wghn, N)
      }

      # update zhat--fix me can probably be vectorized
      if (is.null(MA) && is.null(AR)){
        zhat = 0
      }else{
        S = 0
        for(j in 1:min(t,ARMAorder[2])){
          S = S-Theta[[t]][j]*Inn[,ARMAorder[2]-j+1]
        }
        zhat = S
      }

      set_rand_state(old_state1)

      # update likelihood
      nloglik = nloglik - log(mean(wgh[t,]))
    }

    # likelihood approximation
    if(nreg<1){
      nloglik = nloglik - log(mypdf(xt[1],MargParms))
    }else{
      nloglik = nloglik - log(mypdf(xt[1], ConstMargParm, DynamMargParm[1]))
    }

    # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
    nloglik = nloglik
    #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))


    if (nloglik==Inf | is.na(nloglik)){
      nloglik = 10^8
    }

  return(nloglik)
}

# particle filter for AR
ParticleFilterRes_AR = function(xt,  Regressor, MargParms, ConstMargParm, DynamMargParm,
                                ARMAorder, AR, MA, ParticleNumber, mycdf, epsilon){

    old_state <- get_rand_state()
    on.exit(set_rand_state(old_state))

    # retrieve sample size and particle number
    T1 = length(xt)
    N = ParticleNumber

    # allocate memory for particles
    wgh = matrix(0,T1,N)

    # allocate memory for previous predictors
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
    nloglik = 0 # initialize likelihood

    # First p steps:
    if (ARMAorder[1]>=2){
      for (t in 2:ARMAorder[1]){

        # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
        if (t==2) {
          zhat = ZprevAll[1:(t-1),]*phit
        } else{
          zhat = colSums(ZprevAll[1:(t-1),]*phit)
        }

        # Recompute integral limits
        if(nreg==0){
          a = (qnorm(mycdf(xt[t]-1,t(MargParms)),0,1) - zhat)/rt
          b = (qnorm(mycdf(xt[t],t(MargParms)),0,1) - zhat)/rt
        }else{
          a = (qnorm(mycdf(xt[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
          b = (qnorm(mycdf(xt[t],ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
        }

        # compute random errors from truncated normal
        err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

        # compute the new Z and add it to the previous ones
        znew = rbind(zhat + rt*err, ZprevAll[1:(t-1),])
        ZprevAll[1:t,] = znew

        # recompute weights
        wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
        wprev = wgh[t,]

        # use YW equation to compute estimates of phi and of the erros
        Gt = toeplitz(TacvfAR(AR)[1:t])
        gt = TacvfAR(AR)[2:(t+1)]
        phit = as.numeric(solve(Gt) %*% gt)
        rt =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/TacvfAR(AR)[1]))

      }
    }


    # From p to T1 I dont need to estimate phi anymore
    for (t in (ARMAorder[1]+1):T1){

      # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
      if(ARMAorder[1]>1){# colsums doesnt work for 1-dimensional matrix
        zhat = colSums(ZprevAll*AR)
      }else{
        zhat=ZprevAll*AR
      }

      # compute limits of truncated normal distribution
      if(nreg==0){
        a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - zhat)/rt
        b = as.numeric(qnorm(mycdf(xt[t],  MargParms),0,1) - zhat)/rt
      }else{
        a = as.numeric(qnorm(mycdf(xt[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
        b = as.numeric(qnorm(mycdf(xt[t],  ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
      }

      # draw errors from truncated normal
      err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

      # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
      znew = zhat + rt*err

      # Resampling Step--here the function differs from LikSISGenDist_ARp

      # compute unnormalized weights
      wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)

      # break if I got NA weight
      if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
        nloglik = 10^8
        break
      }

      # normalized weights
      wghn = wgh[t,]/sum(wgh[t,])

      old_state1 <- get_rand_state()

      # sample indices from multinomial distribution-see Step 4 of SISR in paper
      ESS = 1/sum(wghn^2)
      if(ESS<epsilon*N){
        ind = rmultinom(1,N,wghn)
        # sample particles
        znew = rep(znew,ind)

        # use low variance resampling
        #znew = lowVarianceRS(znew, wghn, N)
      }
      set_rand_state(old_state1)


      # save particles
      if (ARMAorder[1]>1){
        ZprevAll = rbind(znew, ZprevAll[1:(ARMAorder[1]-1),])
      }else {
        ZprevAll[1,]=znew
      }
      # update likelihood
      nloglik = nloglik - log(mean(wgh[t,]))
    }

    # likelihood approximation
    if(nreg==0){
      nloglik = nloglik - log(mypdf(xt[1],MargParms))
    }else{
      nloglik = nloglik - log(mypdf(xt[1], ConstMargParm, DynamMargParm[1]))
    }


    # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
    nloglik = nloglik
    #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))


  return(nloglik)
}


# PF likelihood with resampling
ParticleFilter_Res = function(xt,  Regressor, MargParms, ConstMargParm, DynamMargParm,
                              ARMAorder, AR, MA, ParticleNumber, mycdf, epsilon){


  # check Particle number
  if(epsilon > 1 || epsilon<0) stop('Please select a value between 0 and 1 for epsilon.')

  # check Particle number
  if(ParticleNumber<1) stop('Please select a nonegative value for the argument ParticleNumber.')

  # check distributions
  if ( !(CountDist %in%  c("Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial" )))
    stop('The argument CountDist must take one of the following values:
         "Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial".')

  # check ARMAorder
  if(prod(ARMAorder)<0 || length(ARMAorder)!= 2) stop('The argument ARMAorder must have length 2 and can not take negative values.')

  # Mixed ARMA model
  if(ARMAorder[1]>0 && ARMAorder[2]>0) stop('Please specify a pure AR or a pure MA model. ARMA(p,q) models with p>0 and q>0 have not yet been implemented.')

  # Pure AR model
  if(ARMAorder[1]>0 && ARMAorder[2]==0) loglik = ParticleFilterRes_MA(xt,  Regressor, MargParms, ConstMargParm, DynamMargParm,
                                                                                 ARMAorder, AR, MA, ParticleNumber, mycdf, epsilon)
  # Pure MA model or White noise
  if(ARMAorder[1]==0&& ARMAorder[2]>=0) loglik = ParticleFilter_Res_MA(theta, xt, Regressor, ARMAorder,
                                                                       ParticleNumber, CountDist, epsilon)

  return(loglik)
}


# wrapper function to do everything
FitCountModel = function (xt,  Regressor,  CountDist,  ARMAorder, method, initial,
                          ParticleNumber, epsilon, MaxCdf, nHC){

  #==========================================================================================#
  # PURPOSE:  Model count time series with specificied marginal distribution and a dependence
  #           structure.
  #
  # NOTES:    1. We are implementing the methods of the paper "Count Series Modeling with
  #           Gaussian Copulas" that was resurbmitted on August of 2020. An earlier version of
  #           this work appeared udner the title "Latent Gaussian Count
  #           Time Series Modeling"  at:  https://arxiv.org/abs/1811.00203
  #
  # INPUTS:
  #    xt:               count time series data
  #    Regressor:        independent variables
  #    CountDist:        count marginal distribution
  #    ARMAorder:        order of underlying ARMA model
  #    method:           estimation method
  #    initial:          initial parameter
  #    ParticleNumber:   number of particles to be used.
  #    epsilon:          parameter affecting resampling in particle filter method
  #    MaxCdf:           trucnation parameter
  #    nHC:              number of Hermitte coefficients
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    August  2020
  #--------------------------------------------------------------------------#



  # retrieve number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # retrieve sample size
  n = length(xt)


  # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1:(1+nreg),
                           "Negative Binomial"   = 1:(2+nreg),
                           "Mixed Poisson"       = 1:(3+nreg),
                           "Generalized Poisson" = 1:(2+nreg),
                           "Binomial"            = 1:(2+nreg))

  # retrieve marginal cdf and pdf if there is no regressor
  if(nreg<1){
    mycdf = switch(CountDist,
                   "Poisson"             = ppois,
                   "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
                   "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = pGpois,
                   "Binomial"            = pbinom
    )
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = dGpois,
                   "Binomial"            = dbinom
    )
  }else{
    # retrieve marginal cdf and pdf if there is at least one regressor
    mycdf = switch(CountDist,
                   "Poisson"             = function(x, ConstTheta, DynamTheta){             ppois   (x, DynamTheta)},
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
    )
    mypdf = switch(CountDist,
                   "Poisson"             = function(x, ConstTheta, DynamTheta){             dpois   (x, DynamTheta)},
                   "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ dnbinom (x, ConstTheta, 1-DynamTheta)},
                   "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ dGpois  (x, ConstTheta, DynamTheta)}
    )
  }

  # retrieve marginal distribution parameters and number of parameters
  MargParms  = theta[MargParmIndices]
  nMargParms = length(MargParms)
  nparms     = length(theta)

  # retrieve regrssor parameters
  if(nreg>0) beta = MargParms[1:(nreg+1)]

  # check if the number of parameters matches the model setting
  if(nMargParms + sum(ARMAorder)!=nparms) stop('The length of theta does not match the model specification.')

  # If there is a regressor link it to the parameters

  if(CountDist == "Negative Binomial" && nreg>0){
    k = MargParms[nreg+2]
    m = exp(Regressor%*%beta)
    ConstMargParm = 1/k
    DynamMargParm = k*m/(1+k*m)
  }

  if(CountDist == "Generalized Poisson" && nreg>0){
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


  # work only with stable ARMA model
  if (!(checkPoly(AR,MA)[1]=="Causal" && checkPoly(AR,MA)[2]=="Invertible")) return()


  # select objective function according to estimation method
  objFunction = switch(method,
                 "PMLE"   = GaussianMLE,
                 "PF"     = ParticleFilter_Res,
                 "IYW"    = ImpliedYW

  )













































  }




