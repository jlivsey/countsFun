#---------retrieve the model scheme
ModelScheme = function(DependentVar, Regressor, EstMethod, ARMAorder, CountDist, ParticleNumber, epsilon,
                       initialParam, TrueParam, Optimization, OptMethod, OutputType, ParamScheme){

  error = 0
  errorMsg = NULL

  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # retrieve sample size
  n = length(DependentVar)

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
                   "Generalized Poisson" = function(x,theta) { pGpois  (x, theta[1], theta[2])},
                   "Binomial"            = pbinom
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = function(x,theta) { dGpois  (x, theta[1], theta[2])},
                   "Binomial"            = dbinom



    )

    # retrieve marginal inverse cdf
    myinvcdf = switch(CountDist,
                      "Poisson"             = qpois,
                      "Negative Binomial"   = function(x, theta){ qnbinom (x, theta[1], 1-theta[2]) },
                      "Mixed Poisson"       = function(x, theta){ qmixpois(x, theta[1], theta[2], theta[3])},
                      "Generalized Poisson" = function(x, theta){ qGpois  (x, theta[1], theta[2])},
                      "Binomial"            = qbinom

    )
  }else{
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = function(x, ConstMargParm, DynamMargParm){             ppois   (x, DynamMargParm)},
                   "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ pnbinom (x, ConstMargParm, 1-DynamMargParm)},
                   "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ pGpois  (x, ConstMargParm, DynamMargParm)}
    )
    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = function(x, ConstMargParm, DynamMargParm){             dpois   (x, DynamMargParm)},
                   "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ dnbinom (x, ConstMargParm, 1-DynamMargParm)},
                   "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ dGpois  (x, ConstMargParm, DynamMargParm)}
    )

    # retrieve marginal inverse cdf
    myinvcdf = switch(CountDist,
                      "Poisson"             = function(x, ConstMargParm, DynamMargParm){             qpois   (x, DynamMargParm)},
                      "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ qnbinom (x, ConstMargParm, 1-DynamMargParm)},
                      "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ qGpois  (x, ConstMargParm, DynamMargParm)}
    )


  }

  # retrieve marginal distribution parameters
  nMargParms = length(MargParmIndices)
  nparms     = nMargParms + sum(ARMAorder)
  nAR        = ARMAorder[1]
  nMA        = ARMAorder[2]

  if(!is.null(initialParam) && length(initialParam)!=nparms) {
    error = 1
    errorMsg = "The length of the initial parameter doesn't match the model specifications."
  }

  # create the constraints
  if(CountDist =="Poisson"){
    if(nreg==0){
      LB = c(0.001, rep(-Inf, sum(ARMAorder)))
      UB = rep(Inf, sum(ARMAorder)+1)
    }else{
      LB = rep(-Inf, sum(ARMAorder)+nreg+1)
      UB = rep(Inf, sum(ARMAorder)+nreg+1)
    }
  }

  if(CountDist == "Negative Binomial"){
    if(nreg==0){
      LB = c(0.01, 0.01, rep(-Inf, sum(ARMAorder)))
      UB = c(Inf, 0.99,   rep( Inf, sum(ARMAorder)))
    }else{
      LB = c(rep(-Inf, nreg+1), 0.001, rep(-Inf, sum(ARMAorder)))
      UB = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAorder)))
    }
  }


  if(CountDist == "Generalized Poisson"){
    if(nreg==0){
      LB = c(0.001, 0.001, rep(-Inf, sum(ARMAorder)))
      UB = c(Inf, Inf,     rep( Inf, sum(ARMAorder)))
    }else{
      LB = c(rep(-Inf, nreg+1), 0.001, rep(-Inf, sum(ARMAorder)))
      UB = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAorder)))
    }
  }


  out = list(
    mycdf           = mycdf,
    mypdf           = mypdf,
    myinvcdf        = myinvcdf,
    MargParmIndices = MargParmIndices,
    initialParam    = initialParam,
    TrueParam       = TrueParam,
    # ARMAModel       = ARMAModel,
    # MargParm        = MargParm,
    # ARParm          = ARParm,
    # MAParm          = MAParm,
    nMargParms      = nMargParms,
    nAR             = nAR,
    nMA             = nMA,
    n               = n,
    nreg            = nreg,
    # MaxCdf          = MaxCdf,
    # nHC             = nHC,
    CountDist       = CountDist,
    ARMAorder       = ARMAorder,
    ParticleNumber  = ParticleNumber,
    epsilon         = epsilon,
    nparms          = nparms,
    UB              = UB,
    LB              = LB,
    error           = error,
    errorMsg        = errorMsg,
    EstMethod       = EstMethod,
    # maxit         = maxit,
    DependentVar    = DependentVar,
    Regressor       = Regressor,
    Optimization    = Optimization,
    OptMethod       = OptMethod,
    OutputType      = OutputType,
    ParamScheme     = ParamScheme
  )
  return(out)

}

# PF likelihood with resampling for AR(p)
ParticleFilter_Res_AR = function(theta, mod){
  #--------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling
  #           to approximate the likelihood of the
  #           a specified count time series model with an underlying AR(p)
  #           dependence structure. A singloe dummy regression is added here.
  #
  # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
  #           details. A first version of the paer can be found at:
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
  # DATE:    July  2020
  #--------------------------------------------------------------------------#

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))


  #-----------  Step 0: Retrieve values from the mod Structure --------------#

  # marginal parameters
  MargParms        = theta[mod$MargParmIndices]

  # regressor parameters
  if(mod$nreg>0){
    beta  = MargParms[1:(mod$nreg+1)]
    m     = exp(mod$Regressor%*%beta)
  }

  # GLM type parameters
  if(mod$CountDist == "Negative Binomial" && mod$nreg>0){
    ConstMargParm  = 1/MargParms[mod$nreg+2]
    DynamMargParm  = MargParms[mod$nreg+2]*m/(1+MargParms[mod$nreg+2]*m)
  }

  if(mod$CountDist == "Generalized Poisson" && mod$nreg>0){
    ConstMargParm  = MargParms[mod$nreg+2]
    DynamMargParm  = m
  }

  if(mod$CountDist == "Poisson" && mod$nreg>0){
    ConstMargParm  = NULL
    DynamMargParm  = m
  }

  # ARMA parameters
  AR = NULL
  if(mod$ARMAorder[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAorder[1])]

  MA = NULL
  if(mod$ARMAorder[2]>0) MA = theta[(mod$nMargParms+mod$ARMAorder[1]+1) :
                                      (mod$nMargParms + mod$ARMAorder[1] + mod$ARMAorder[2]) ]

  # check for causality
  if( CheckStability(AR,MA) ) return(10^(8))


  # sample size and number of particles
  T1      = length(mod$DependentVar)
  N       = mod$ParticleNumber

  # Initialize the negative log likelihood computation
  nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],MargParms)),
                   - log(mod$mypdf(mod$DependentVar[1], ConstMargParm, DynamMargParm[1])))

  # Compute the theoretical covariance for the AR model for current estimate
  gt      = ARMAacf(ar = AR, ma = MA)[2:(max(mod$ARMAorder)+1)]

  # Compute the best linear predictor coefficients and errors using Durbin Levinson
  DL      = DLAcfToAR(gt)
  phit    = DL[,1]
  Rt      = sqrt(DL[,3])


  # allocate memory for particle weights and the latent Gaussian Series particles
  wgh     = matrix(0,T1,N)
  Z       = matrix(0,mod$ARMAorder[1],N)

  #======================   Start the SIS algorithm   ======================#
  # Initialize the weights and the latent Gaussian series particles
  wgh[1,] = rep(1,N)
  a       = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1,t(MargParms)),0,1), N)
  b       = rep( qnorm(mod$mycdf(mod$DependentVar[1],t(MargParms)),0,1), N)
  Z[1,]   = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

  # =================== Loop from 2 to AR order===================== #
  if (mod$ARMAorder[1]>=2){
    for (t in 2: (mod$ARMAorder[1])){
      # STEP 1 in SIS: Compute the latent Gaussian predictions Zhat using Durbin Levinson
      if (t==2) {
        Zhat = Z[1:(t-1),]*phit[1:(t-1)]
      } else{
        Zhat = colSums(Z[1:(t-1),]*phit[1:(t-1)])
      }

      # STEP 2 is SIS: Update the latent Gaussian series Z and the importance weights w
      if(mod$nreg==0){
        a = (qnorm(mod$mycdf(mod$DependentVar[t]-1,t(MargParms)),0,1) - Zhat)/Rt[t]
        b = (qnorm(mod$mycdf(mod$DependentVar[t],t(MargParms)),0,1) - Zhat)/Rt[t]
      }else{
        a = (qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - Zhat)/Rt[t]
        b = (qnorm(mod$mycdf(mod$DependentVar[t],ConstMargParm, DynamMargParm[t]),0,1) - Zhat)/Rt[t]
      }

      Z[1:t,] = rbind(qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)*Rt[t] + Zhat, Z[1:(t-1),])
      wgh[t,] = wgh[t-1,]*(pnorm(b,0,1) - pnorm(a,0,1))

      # update likelihood
      nloglik = nloglik - log(mean(wgh[t,]))
      # print(t)
      # print(nloglik)
    }
  }
  # =================== Loop from AR order + 1  to T ===================== #
  # From p to T1 I don't need to estimate phi anymore
  for (t in (mod$ARMAorder[1]+1):T1){

    # STEP 1 in SIS: Compute the latent Gaussian predictions Zhat using Durbin Levinson
    if(mod$ARMAorder[1]>1){# colsums doesnt work for 1-dimensional matrix
      Zhat = colSums(Z*phit)
    }else{
      Zhat =  Z*phit
    }

    # STEP 2 is SISR: Update the latent Gaussian series Z
    if(mod$nreg==0){
      a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,t(MargParms)),0,1)) - Zhat)/Rt[mod$ARMAorder[1]]
      b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],t(MargParms)),0,1)) - Zhat)/Rt[mod$ARMAorder[1]]
    }else{
      a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1)) - Zhat)/Rt[mod$ARMAorder[1]]
      b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],ConstMargParm, DynamMargParm[t]),0,1)) - Zhat)/Rt[mod$ARMAorder[1]]
    }

    Znew = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)*Rt[mod$ARMAorder[1]] + Zhat

    # compute unnormalized weights
    # wgh[t,] = wgh[t-1,]*(pnorm(b,0,1) - pnorm(a,0,1))
    wgh[t,] = (pnorm(b,0,1) - pnorm(a,0,1))
    # break if I got NA weight
    if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
      message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
      return(10^8)
    }

    # compute normalized weights
    wghn = wgh[t,]/sum(wgh[t,])

    # STEP 3 is SISR: Resample
    old_state1 <- get_rand_state()
    ESS = 1/sum(wghn^2)
    if(ESS<mod$epsilon*N){
      ind = rmultinom(1,N,wghn)
      # sample particles
      Znew = rep(Znew,ind)
    }
    set_rand_state(old_state1)


    # save particles
    if (mod$ARMAorder[1]>1){
      Z = rbind(Znew, Z[1:(mod$ARMAorder[1]-1),])
    }else {
      Z[1,]=Znew
    }
    # update likelihood
    nloglik = nloglik - log(mean(wgh[t,]))
    # print(t)
    # print(nloglik)
  }

  return(nloglik)
}

# PF likelihood with resampling for MA(q)
ParticleFilter_Res_MA = function(theta, mod){
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
  #    Regressor:        independent variable
  #    CountDist:        count marginal distribution
  #    epsilon           resampling when ESS<epsilon*N
  #
  # OUTPUT:
  #    loglik:           approximate log-likelihood
  #
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    July 2020
  #------------------------------------------------------------------------------------#

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # retrieve marginal parameters
  MargParms        = theta[mod$MargParmIndices]

  # retrieve regressor parameters
  if(mod$nreg>0){
    beta  = MargParms[1:(mod$nreg+1)]
    m     = exp(mod$Regressor%*%beta)
  }

  # retrieve GLM type  parameters
  if(mod$CountDist == "Negative Binomial" && mod$nreg>0){
    ConstMargParm  = 1/MargParms[mod$nreg+2]
    DynamMargParm  = MargParms[mod$nreg+2]*m/(1+MargParms[mod$nreg+2]*m)
  }

  if(mod$CountDist == "Generalized Poisson" && mod$nreg>0){
    ConstMargParm  = MargParms[mod$nreg+2]
    DynamMargParm  = m
  }

  if(mod$CountDist == "Poisson" && mod$nreg>0){
    ConstMargParm  = NULL
    DynamMargParm  = m
  }

  # retrieve ARMA parameters
  AR = NULL
  if(mod$ARMAorder[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAorder[1])]

  MA = NULL
  if(mod$ARMAorder[2]>0) MA = theta[(mod$nMargParms+mod$ARMAorder[1]+1) : (mod$nMargParms + mod$ARMAorder[1] + mod$ARMAorder[2]) ]

  # check for causality
  if( CheckStability(AR,MA) ) return(10^(-6))


  T1 = length(mod$DependentVar)
  N = mod$ParticleNumber          # number of particles


  # allocate matrix to collect all particle weights
  wgh = matrix(0,length(mod$DependentVar),N)

  # Compute integral limits
  if(mod$nreg==0){
    a = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1,t(MargParms)),0,1), N)
    b = rep( qnorm(mod$mycdf(mod$DependentVar[1],t(MargParms)),0,1), N)
  }else{
    a = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
    b = rep( qnorm(mod$mycdf(mod$DependentVar[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
  }

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)


  # run innovations Algorithm for MA models that are not WN
  if(mod$ARMAorder[2]>0) Inn = matrix(0,N,mod$ARMAorder[2])   # I will save here the q many innovations (Z - Zhat) --see (5.3.9) BD book
  if (is.null(MA) && is.null(AR)){
    v0   = 1
    zhat = 0
  }else{
    # FIX ME: Check if the code below is correct in terms of the ARMAacf
    MA.acvf <- as.vector(ARMAacf(ma = MA, lag.max=T1))
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
    if(mod$ARMAorder[2]>1) Inn[,1:(mod$ARMAorder[2]-1)] = Inn[,2:(mod$ARMAorder[2])]

    # compute limits of truncated normal distribution
    if(mod$nreg==0){
      a = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t]-1,MargParms),0,1) - zhat)/vt
      b = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t],MargParms),0,1) -   zhat)/vt
    }else{
      a = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/vt
      b = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t],ConstMargParm, DynamMargParm[t]),0,1) -   zhat)/vt
    }

    # draw errors from truncated normal
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

    # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
    znew = zhat + vt*err

    # compute new innovation
    Inn[,mod$ARMAorder[2]] = (znew-zhat)

    # compute unnormalized weights
    wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)

    # break if I got NA weight
    if (any(is.na(wgh[t,]))| sum(wgh[t,])<10^(-8) ){
      nloglik = 10^8
      break
    }

    # normalized weights
    wghn = wgh[t,]/sum(wgh[t,])


    # Resampling: sample indices from multinomial distribution-see Step 4 of SISR in paper
    ESS = 1/sum(wghn^2)
    old_state1 <- get_rand_state()
    if(ESS<mod$epsilon*N){
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
      for(j in 1:min(t,mod$ARMAorder[2])){
        S = S-Theta[[t]][j]*Inn[,mod$ARMAorder[2]-j+1]
      }
      zhat = S
    }

    set_rand_state(old_state1)

    # update likelihood
    nloglik = nloglik - log(mean(wgh[t,]))
  }

  # likelihood approximation
  if(mod$nreg<1){
    nloglik = nloglik - log(mod$mypdf(mod$DependentVar[1],MargParms))
  }else{
    nloglik = nloglik - log(mod$mypdf(mod$DependentVar[1], ConstMargParm, DynamMargParm[1]))
  }

  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))


  if (nloglik==Inf | is.na(nloglik)){
    nloglik = 10^8
  }


  return(nloglik)
}


# PF likelihood with resampling
ParticleFilter_Res = function(theta, mod){
  #--------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling
  #           to approximate the likelihood of the
  #           a specified count time series model with an underlying AR(p)
  #           dependence structure or MA(q) structure.
  #
  # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
  #           details. A first version of the paer can be found at:
  #           https://arxiv.org/abs/1811.00203

  #
  # INPUTS:
  #    theta:            parameter vector
  #    data:             dependent variable
  #    Regressor:        independent variables
  #    ParticleNumber:   number of particles to be used in likelihood approximation
  #    CountDist:        count marginal distribution
  #    epsilon           resampling when ESS<epsilon*N
  #
  # OUTPUT:
  #    loglik:           approximate log-likelihood
  #
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    July 2020
  #--------------------------------------------------------------------------#


  # Pure AR model
  if(mod$ARMAorder[1]>0 && mod$ARMAorder[2]==0) loglik = ParticleFilter_Res_AR(theta, mod)
  # Pure MA model or White noise
  if(mod$ARMAorder[1]==0&& mod$ARMAorder[2]>=0) loglik = ParticleFilter_Res_MA(theta, mod)

  return(loglik)
}


# Optimization wrapper to fit PF likelihood with resampling
FitMultiplePF_Res = function(theta, mod){
  #====================================================================================#
  # PURPOSE       Fit the Particle Filter log-likelihood with resampling.
  #               This function maximizes the PF likelihood, nfit many times for nparts
  #               many choices of particle numbers, thus yielding a total of nfit*nparts
  #               many estimates.
  #
  # INPUT
  #   theta       initial parameters
  #   data        count series
  #   CountDist   prescribed count distribution
  #   Particles   vector with different choices for number of particles
  #   ARMAorder   order of the udnerlying ARMA model
  #   epsilon     resampling when ESS<epsilon*N
  #   LB          parameter lower bounds
  #   UB          parameter upper bounds
  #   OptMethod
  #
  #
  # OUTPUT
  #   All         parameter estimates, standard errors, likelihood value, AIC, etc
  #
  # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date          April 2020
  # Version       3.6.3
  #====================================================================================#


  # retrieve parameter, sample size etc
  nparts   = length(mod$ParticleNumber)
  nparms   = length(theta)
  nfit     = 1
  n        = length(mod$DependentVar)

  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst  = matrix(0,nrow=nfit*nparts,ncol=nparms)
  se       = matrix(NA,nrow=nfit*nparts,ncol=nparms)
  loglik   = rep(0,nfit*nparts)
  convcode = rep(0,nfit*nparts)
  kkt1     = rep(0,nfit*nparts)
  kkt2     = rep(0,nfit*nparts)

  # Each realization will be fitted nfit*nparts many times
  for (j in 1:nfit){
    set.seed(j)
    # for each fit repeat for different number of particles
    for (k in 1:nparts){
      # FIX ME: I need to somehow update this in mod. (number of particles to be used). I t now works only for 1 choice of ParticleNumber
      ParticleNumber = mod$ParticleNumber[k]

      if(mod$Optimization == TRUE){
      # run optimization for our model --no ARMA model allowed
      optim.output <- optimx(
        par     = theta,
        fn      = ParticleFilter_Res,
        lower   = mod$LB,
        upper   = mod$UB,
        hessian = TRUE,
        method  = mod$OptMethod,
        mod     = mod)
      }else{
        optim.output = as.data.frame(matrix(rep(NA,8+length(theta)), ncol=8+length(theta)))
        names(optim.output) = c(mapply(sprintf, rep("p%.f",length(theta)), (1:length(theta)), USE.NAMES = FALSE),
                                "value",  "fevals", "gevals", "niter", "convcode",  "kkt1", "kkt2", "xtime")

        optim.output[,1:length(theta)] = theta
        t1 = tic()
        loglikelihood = ParticleFilter_Res(theta,mod)
        t2 = tic()
        optim.output[,(length(theta)+1)] = loglikelihood
        optim.output[,(length(theta)+2)] = 1
        optim.output[,(length(theta)+3)] = 1
        optim.output[,(length(theta)+4)] = 0
        optim.output[,(length(theta)+8)] = as.numeric(t2-t1)
      }


      # save estimates, loglik value and diagonal hessian
      ParmEst[nfit*(k-1)+j,]  = as.numeric(optim.output[1:nparms])
      loglik[nfit*(k-1) +j]   = optim.output$value
      convcode[nfit*(k-1) +j] = optim.output$convcode
      kkt1[nfit*(k-1) +j]     = optim.output$kkt1
      kkt2[nfit*(k-1) +j]     = optim.output$kkt2


      # compute Hessian
      H = gHgen(par            = ParmEst[nfit*(k-1)+j,],
                fn             = ParticleFilter_Res,
                mod            = mod)

      # if I get all na for one row and one col of H remove it
      # H$Hn[rowSums(is.na(H$Hn)) != ncol(H$Hn), colSums(is.na(H$Hn)) != nrow(H$Hn)]

      if (!any(is.na(rowSums(H$Hn)))){
        # save standard errors from Hessian
        if(H$hessOK && det(H$Hn)>10^(-8)){
          se[nfit*(k-1)+j,]   = sqrt(abs(diag(solve(H$Hn))))
        }else{
          se[nfit*(k-1)+j,] = rep(NA, nparms)
        }
      }else{
        # remove the NA rows and columns from H
        Hnew = H$Hn[rowSums(is.na(H$Hn)) != ncol(H$Hn), colSums(is.na(H$Hn)) != nrow(H$Hn)]

        # find which rows are missing and which are not
        NAIndex = which(colSums(is.na(H$Hn))==nparms)
        NonNAIndex = which(colSums(is.na(H$Hn))==1)

        #repeat the previous ifelse for the reduced H matrix
        if(det(Hnew)>10^(-8)){
          se[nfit*(k-1)+j,NonNAIndex]   = sqrt(abs(diag(solve(Hnew))))
        }else{
          se[nfit*(k-1)+j,NAIndex] = rep(NA, length(NAindex))
        }

      }


    }
  }

  # Compute model selection criteria (assuming one fit)
  Criteria = Criteria.lgc(loglik, mod)


  # get the names of the final output
  parmnames     = colnames(optim.output)
  # mynames       = c(parmnames[1:nparms],paste("se", parmnames[1:nparms], sep="_"), "loglik", "AIC", "BIC","AICc", "ConvergenceStatus", "kkt1", "kkt2")
  # All           = matrix(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)
  # colnames(All) = mynames
  # All           = data.frame(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)

  if(mod$OutputType=="list"){
    #  save the results in a list
    ModelOutput = list(data.frame(matrix(rep(NA,nparms),nrow=1)),
                      data.frame(matrix(rep(NA,nparms),nrow=1)),
                      data.frame(matrix(rep(NA,4),nrow=1)),
                      data.frame(matrix(rep(NA,3),nrow=1)))
    names(ModelOutput)         = c("ParamEstimates", "StdErrors", "FitStatistics", "OptimOutput")
    ModelOutput$ParamEstimates = ParmEst
    ModelOutput$StdErrors      = se
    ModelOutput$FitStatistics  = c(loglik, Criteria)
    ModelOutput$OptimOutput    = c(convcode,kkt1,kkt2)
    ModelOutput$CountDist      = mod$CountDist
    ModelOutput$EstMethod      = mod$EstMethod
    ModelOutput$ARMAModel      = mod$ARMAorder
    ModelOutput$Optimization   = mod$Optimization
  }else{
    ModelOutput  = data.frame(matrix(ncol = 4*mod$nparms+16, nrow = 1))

    # names of the output data frame
    colnames(ModelOutput) = c(
      'CountDist','ARMAModel', 'Regressor',
      paste("True", parmnames[1:nparms], sep="_"), paste("InitialEstim", parmnames[1:nparms], sep="_"),
      parmnames[1:nparms], paste("se", parmnames[1:nparms], sep="_"),
      'EstMethod', 'SampleSize', 'ParticleNumber', 'epsilon', 'OptMethod', 'ParamScheme',
      "loglik", "AIC", "BIC", "AICc", "ConvergeStatus", "kkt1", "kkt2")

    # Start Populating the output data frame
    ModelOutput$CountDist      = mod$CountDist
    ModelOutput$ARMAModel      = paste(mod$ARMAorder, collapse=' ')
    ModelOutput$Regressor      = !is.null(mod$Regressor)


    # true Parameters
    offset = 4
    if(mod$nMargParms>0){
      ModelOutput[, offset:(offset + mod$nMargParms -1)] = mod$TrueParam[1:mod$nMargParms]
      offset = offset + mod$nMargParms
    }

    if(mod$nAR>0){
      ModelOutput[, offset:(offset + mod$nAR -1)]        = mod$TrueParam[ (mod$nMargParms+1):(mod$nMargParms+mod$nAR) ]
      offset = offset + mod$nAR
    }

    if(mod$nMA>0){
      ModelOutput[, offset:(offset + mod$nMA -1)]        = mod$TrueParam[ (mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA)]
      offset = offset + mod$nMA
     }

    # Initial Parameter Estimates
    if(mod$nMargParms>0){
      ModelOutput[, offset:(offset + mod$nMargParms -1 )] = mod$initialParam[1:mod$nMargParms ]
      offset = offset + mod$nMargParms
    }
    if(mod$nAR>0){
      ModelOutput[, offset:(offset + mod$nAR-1)]        = mod$initialParam[(mod$nMargParms+1):(mod$nMargParms+mod$nAR) ]
      offset = offset + mod$nAR
    }

    if(mod$nMA>0){
      ModelOutput[, offset:(offset + mod$nMA-1)]        = mod$initialParam[(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA) ]
      offset = offset + mod$nMA
    }

    # Parameter Estimates
    if(mod$nMargParms>0){
      ModelOutput[, offset:(offset + mod$nMargParms-1)] = ParmEst[,1:mod$nMargParms]
      offset = offset + mod$nMargParms
    }

    if(mod$nAR>0){
      ModelOutput[, offset:(offset + mod$nAR-1)]        = ParmEst[,(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]
      offset = offset + mod$nAR
    }

    if(mod$nMA>0){
      ModelOutput[, offset:(offset + mod$nMA-1)]        = ParmEst[,(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA)]
      offset = offset + mod$nMA
    }

    # Parameter Std Errors
    if(mod$nMargParms>0){
      ModelOutput[, offset:(offset + mod$nMargParms-1)] = se[,1:mod$nMargParms]
      offset = offset + mod$nMargParms
    }

    if(mod$nAR>0){
      ModelOutput[, offset:(offset + mod$nAR-1)]        = se[,(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]
      offset = offset + mod$nAR
    }

    if(mod$nMA>0){
      ModelOutput[, offset:(offset + mod$nMA-1)]        = se[,(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA)]
    }

    ModelOutput$EstMethod      = mod$EstMethod
    ModelOutput$SampleSize     = mod$n
    ModelOutput$ParticleNumber = mod$ParticleNumber
    ModelOutput$epsilon        = mod$epsilon
    ModelOutput$OptMethod      = row.names(optim.output)
    ModelOutput$ParamScheme    = mod$ParamScheme
    ModelOutput$loglik         = loglik
    ModelOutput$AIC            = Criteria[1]
    ModelOutput$BIC            = Criteria[2]
    ModelOutput$AICc           = Criteria[3]
    ModelOutput$ConvergeStatus = convcode
    ModelOutput$kkt1           = kkt1
    ModelOutput$kkt2           = kkt2

  }



  return(ModelOutput)
}

# Add some functions that I will need
get_rand_state <- function() {
  # Using `get0()` here to have `NULL` output in case object doesn't exist.
  # Also using `inherits = FALSE` to get value exactly from global environment
  # and not from one of its parent.
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

set_rand_state <- function(state) {
  # Assigning `NULL` state might lead to unwanted consequences
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
  }
}


#---------check causality and invertibility
CheckStability = function(AR,MA){
  if (is.null(AR) && is.null(MA)) return(0)

  # return 1 if model is not stable (causal and invertible)
  if(!is.null(AR) && is.null(MA)){
    rc = ifelse(any(abs( polyroot(c(1, -AR))  ) < 1), 1,0)
  }

  if(!is.null(MA) && is.null(AR)){
    rc = ifelse(any(abs( polyroot(c(1, -MA))  ) < 1),1,0)
  }

  if(!is.null(MA) && !is.null(AR)){
    rc = ifelse(checkPoly(AR,MA)[1]!="Causal" && check(poly)[2]!="Invertible", 1,0)
  }

  return(rc)
}


# compute initial estimates
InitialEstimates = function(mod){
  # require(itsmr)
  est  = rep(NA, mod$nMargParms+sum(mod$ARMAorder))
  #-----Poisson case
  if(mod$CountDist=="Poisson"){
    if(mod$nreg==0){
      est[1] = mean(mod$DependentVar)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      # GLM for the mean that depends on time
      # CHECK ME: If I fit a Poisson AR(3) in the the data example of the JASA paper, but the code below doesn't specify poisson family (it would pick up the default distribution that glm function has) then there will be a numerical error in the likelihood. Check it!
      glmPoisson            = glm(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)], family = "poisson")
      est[1:mod$nMargParms] = as.numeric(glmPoisson[1]$coef)

      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAorder))] = arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }

  #-----Neg Binomial case
  if(mod$CountDist=="Negative Binomial"){
    if(mod$nreg==0){
      xbar = mean(mod$DependentVar)
      sSquare = var(mod$DependentVar)

      # Method of Moments for negBin
      rEst = xbar^2/(sSquare - xbar)
      pEst = 1 - xbar/sSquare
      est[1:2] = c(rEst, pEst)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      # GLM for the mean that depends on time
      glmNegBin                 = glm.nb(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)])
      est[1:(mod$nMargParms-1)] = as.numeric(glmNegBin[1]$coef)
      # Mom on constant variance
      est[mod$nMargParms]       = NegBinMoM(mod$DependentVar,glmNegBin$fitted.values)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$theta
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

      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(mod$DependentVar~mod$Regressor)[1]$coef)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }



  #-----Generalized Poisson case
  if(mod$CountDist=="Generalized Poisson"){
    if(mod$nreg==0){
      xbar = mean(mod$DependentVar)
      sSquare = var(mod$DependentVar)

      # Method of Moments for negBin
      rEst = xbar^2/(sSquare - xbar)
      pEst = 1 - xbar/sSquare
      est[1:2] = c(rEst, pEst)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(mod$DependentVar~mod$Regressor)[1]$coef)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$DependentVar,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }





  return(est)
}


NegBinMoM = function(data, GLMMeanEst){
  # the GLMMeanEst is the GLM estimate of the standard log-link
  # th formula below is standard MoM for the over dispersion param in NegBin2 parametrization
  PhiMomEst = sum(GLMMeanEst^2)/(sum((data-GLMMeanEst)^2-GLMMeanEst))
  return(PhiMomEst)
}

# innovations algorithm code
innovations.algorithm <- function(acvf,n.max=length(acvf)-1){
  # Found this onlinbe need to check it
  # http://faculty.washington.edu/dbp/s519/R-code/innovations-algorithm.R
  thetas <- vector(mode="list",length=n.max)
  vs <- rep(acvf[1],n.max+1)
  for(n in 1:n.max){
    thetas[[n]] <- rep(0,n)
    thetas[[n]][n] <- acvf[n+1]/vs[1]
    if(n>1){
      for(k in 1:(n-1)){
        js <- 0:(k-1)
        thetas[[n]][n-k] <- (acvf[n-k+1] - sum(thetas[[k]][k-js]*thetas[[n]][n-js]*vs[js+1]))/vs[k+1]
      }
    }
    js <- 0:(n-1)
    vs[n+1] <- vs[n+1] - sum(thetas[[n]][n-js]^2*vs[js+1])
  }
  return(structure(list(vs=vs,thetas=thetas)))
}

# simulate from our model
sim_lgc = function(n, CountDist, MargParm, ARParm, MAParm){

  myinvcdf = switch(CountDist,
                    "Poisson"             = qpois,
                    "Negative Binomial"   = function(x, theta){ qnbinom (x, theta[1], 1-theta[2]) },
                    "Mixed Poisson"       = function(x, theta){ qmixpois(x, theta[1], theta[2], theta[3])},
                    "Generalized Poisson" = function(x, theta){ qGpois  (x, theta[1], theta[2])},
                    "Binomial"            = qbinom
  )

  # Generate latent Gaussian model
  z  =arima.sim(model = list( ar = ARParm, ma=MAParm  ), n = n)
  # Generate the count data (without regressor)
  x = myinvcdf(pnorm(z), MargParm)
  return(x)
}


#---------Compute AIC, BIC, AICc
Criteria.lgc = function(loglik, mod){
  #---------------------------------#
  # Purpose:     Compute AIC, BIC, AICc for our models
  #
  # Inputs:
  #  loglik      log likelihood value at optimum
  #  nparms      number of parametersin the model
  #       n      sample size
  #
  # NOTES:       The Gaussian likelihood we are minimizing has the form:
  #              l0 = 0.5*logdet(G) + 0.5*X'*inv(G)*X. We will need to
  #              bring this tothe classic form before we compute the
  #              criteria (see log of lrealtion (8.6.1) in Brockwell and Davis)
  #
  #              No correction is necessary in the PF case
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         July 2020
  # Version      3.6.3
  #---------------------------------#

  if(mod$EstMethod!="GL"){
    l1 = -loglik
  }else{
    l1 = -loglik  - mod$n/2*log(2*pi)
  }

  AIC = 2*mod$nparms - 2*l1
  BIC = log(mod$n)*mod$nparms - 2*l1
  AICc = AIC + (2*mod$nparms^2 + 2*mod$nparms)/(mod$n-mod$nparms-1)

  AllCriteria = c(AIC, BIC, AICc)
}

logLik.lgc = function(object){
  return(object$FitStatistics[1])
}

AIC.lgc = function(object){
  return(object$FitStatistics[2])
}

BIC <- function(object, ...) UseMethod("BIC")

BIC.lgc = function(object){
  return(object$FitStatistics[3])
}

se <- function(object, ...) UseMethod("se")

se.lgc = function(object){
  return(object$StdErrors)
}

coefficients <- function(object, ...) UseMethod("coefficients")

coefficients.lgc = function(object){
  return(object$ParamEstimates)
}

model <- function(object, ...) UseMethod("model")

model.lgc = function(object){
  if ((object$ARMAModel[1]>0) &&  (object$ARMAModel[2]>0)){
    ARMAModel = sprintf("ARMA(%.0f, %.0f)",object$ARMAModel[1], object$ARMAModel[2])
  }
  if ((object$ARMAModel[1]>0) &&  (object$ARMAModel[2]==0)){
    ARMAModel = sprintf("AR(%.0f)",object$ARMAModel[1])
  }
  if ((object$ARMAModel[1]==0) &&  (object$ARMAModel[2]>0)){
    ARMAModel = sprintf("MA(%.0f)",object$ARMAModel[2])
  }
  if ((object$ARMAModel[1]==0) &&  (object$ARMAModel[2]==0)){
    ARMAModel = "White Noise"
  }

  a = data.frame(object$CountDist, ARMAModel)
  names(a) = c("Distribution", "Model")
  return(a)
}



# PF likelihood with resampling for AR(p)
# ParticleFilter_Res_AR_old = function(theta, mod){
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
#   # retrieve marginal parameters
#   MargParms        = theta[mod$MargParmIndices]
#
#   # retrieve regressor parameters
#   if(mod$nreg>0){
#     beta  = MargParms[1:(mod$nreg+1)]
#     m     = exp(mod$Regressor%*%beta)
#   }
#
#   # retrieve GLM type  parameters
#   if(mod$CountDist == "Negative Binomial" && mod$nreg>0){
#     ConstMargParm  = 1/MargParms[mod$nreg+2]
#     DynamMargParm  = MargParms[mod$nreg+2]*m/(1+MargParms[mod$nreg+2]*m)
#   }
#
#   if(mod$CountDist == "Generalized Poisson" && mod$nreg>0){
#     ConstMargParm  = MargParms[mod$nreg+2]
#     DynamMargParm  = m
#   }
#
#   if(mod$CountDist == "Poisson" && mod$nreg>0){
#     ConstMargParm  = NULL
#     DynamMargParm  = m
#   }
#
#   # retrieve ARMA parameters
#   AR = NULL
#   if(mod$ARMAorder[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAorder[1])]
#
#   MA = NULL
#   if(mod$ARMAorder[2]>0) MA = theta[(mod$nMargParms+mod$ARMAorder[1]+1) : (mod$nMargParms + mod$ARMAorder[1] + mod$ARMAorder[2]) ]
#
#   # check for causality
#   if( CheckStability(AR,MA) ) return(10^(8))
#
#
#   T1 = length(data)
#   N = mod$ParticleNumber          # number of particles
#
#   wgh = matrix(0,T1,N)        # to collect all particle weights
#
#   # allocate memory for zprev
#   ZprevAll = matrix(0,mod$ARMAorder[1],N)
#
#   if(mod$nreg==0){
#     # Compute integral limits
#     a = rep( qnorm(mod$mycdf(data[1]-1,t(MargParms)),0,1), N)
#     b = rep( qnorm(mod$mycdf(data[1],t(MargParms)),0,1), N)
#   }else{
#     a = rep( qnorm(mod$mycdf(data[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
#     b = rep( qnorm(mod$mycdf(data[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
#   }
#   # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#   zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#   # save the currrent normal variables
#   ZprevAll[1,] = zprev
#
#   # initial estimate of first AR coefficient as Gamma(1)/Gamma(0) and corresponding error
#   # phit = TacvfAR(AR)[2]/TacvfAR(AR)[1] - this FitAR package is obsolete
#   # FIX ME: check if the code below is correct
#   phit = ARMAacf(ar = AR)[2]/ARMAacf(ar = AR)[1]
#   rt = as.numeric(sqrt(1-phit^2))
#
#   # particle filter weights
#   wprev = rep(1,N)
#   wgh[1,] = wprev
#
#   # Initiate log likelihood
#   nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],MargParms)),
#                    - log(mod$mypdf(mod$DependentVar[1], ConstMargParm, DynamMargParm[1])))
#   # print(nloglik)
#   # First p steps:
#   if (mod$ARMAorder[1]>=2){
#     for (t in 2:mod$ARMAorder[1]){
#
#       # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
#       if (t==2) {
#         zhat = ZprevAll[1:(t-1),]*phit
#       } else{
#         zhat = colSums(ZprevAll[1:(t-1),]*phit)
#       }
#
#       # Recompute integral limits
#       if(mod$nreg==0){
#         a = (qnorm(mod$mycdf(data[t]-1,t(MargParms)),0,1) - zhat)/rt
#         b = (qnorm(mod$mycdf(data[t],t(MargParms)),0,1) - zhat)/rt
#       }else{
#         a = (qnorm(mod$mycdf(data[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#         b = (qnorm(mod$mycdf(data[t],ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#       }
#
#       # compute random errors from truncated normal
#       err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#       # compute the new Z and add it to the previous ones
#       znew = rbind(zhat + rt*err, ZprevAll[1:(t-1),])
#       ZprevAll[1:t,] = znew
#
#       # recompute weights
#       wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
#       wprev = wgh[t,]
#
#       # use YW equation to compute estimates of phi and of the errors
#       # FIX ME: Check if the code below is correct in terms of the ARMAacf
#       Gt    = toeplitz(ARMAacf(AR)[1:t])
#       gt    = ARMAacf(AR)[2:(t+1)]
#       phit  = as.numeric(solve(Gt) %*% gt)
#       rt    =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/ARMAacf(AR)[1]))
#
#       # update likelihood
#       nloglik = nloglik - log(mean(wgh[t,]))
#       # print(t)
#       # print(nloglik)
#     }
#   }
#
#
#   # From p to T1 I dont need to estimate phi anymore
#   for (t in (mod$ARMAorder[1]+1):T1){
#
#     # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
#     if(mod$ARMAorder[1]>1){# colsums doesnt work for 1-dimensional matrix
#       zhat = colSums(ZprevAll*AR)
#     }else{
#       zhat=ZprevAll*AR
#     }
#
#     # compute limits of truncated normal distribution
#     if(mod$nreg==0){
#       a = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t]-1,MargParms),0,1) - zhat)/rt
#       b = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t],  MargParms),0,1) - zhat)/rt
#     }else{
#       a = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#       b = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t],  ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
#     }
#
#     # draw errors from truncated normal
#     err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#
#     # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
#     znew = zhat + rt*err
#
#     # Resampling Step--here the function differs from LikSISGenDist_ARp
#
#     # compute unnormalized weights
#     wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)
#
#     # break if I got NA weight
#     if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
#       message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
#       return(10^8)
#       # I am fitting a poisson with initialParam= NULL in the Sales data and the initial estimated
#       # mean form GLM, returns parameters that yield a, and b above that lead to Zprev = -Inf. I will comment
#       # out the break and return, also change the value from 10^8 to 10^18.
#       # return(nloglik)
#     }
#
#     # normalized weights
#     wghn = wgh[t,]/sum(wgh[t,])
#
#     old_state1 <- get_rand_state()
#
#     # sample indices from multinomial distribution-see Step 4 of SISR in paper
#     ESS = 1/sum(wghn^2)
#     if(ESS<mod$epsilon*N){
#       ind = rmultinom(1,N,wghn)
#       # sample particles
#       znew = rep(znew,ind)
#
#       # use low variance resampling
#       #znew = lowVarianceRS(znew, wghn, N)
#     }
#     set_rand_state(old_state1)
#
#
#     # save particles
#     if (mod$ARMAorder[1]>1){
#       ZprevAll = rbind(znew, ZprevAll[1:(mod$ARMAorder[1]-1),])
#     }else {
#       ZprevAll[1,]=znew
#     }
#     # update likelihood
#     nloglik = nloglik - log(mean(wgh[t,]))
#     # print(t)
#     # print(nloglik)
#   }
#
#   # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#   # nloglik = nloglik
#   #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#   # print(t)
#   # print(nloglik)
#   return(nloglik)
# }
#





