#---------retrieve the model scheme
ModelScheme = function(DependentVar, Regressor, EstMethod, ARMAModel, CountDist, ParticleNumber, epsilon,
                       initialParam, TrueParam=NULL, Task, SampleSize, OptMethod, OutputType, ParamScheme){

  error = 0
  errorMsg = NULL

  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # retrieve sample size
  n = ifelse(!is.null(DependentVar), length(DependentVar), SampleSize)


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
  nparms     = nMargParms + sum(ARMAModel)
  nAR        = ARMAModel[1]
  nMA        = ARMAModel[2]

  # if(!is.null(initialParam) && length(initialParam)!=nparms) {
  #   error = 1
  #   errorMsg = "The length of the initial parameter doesn't match the model specifications."
  # }

  # create names of parameters that will be used for output - start with marginal parameters
  if(nreg<1){
    MargParmsNames = switch(CountDist,
                            "Poisson"             = c("lambda"),
                            "Negative Binomial"   = c("r","p"),
                            "Mixed Poisson"       = c("lambda_1", "lambda_2", "p"),
                            "Generalized Poisson" = c("lambda", "a"),
                            "Binomial"            = c("n", "p")
    )
  }else{
    # fix me: check mixed poisson case
    MargParmsNames = switch(CountDist,
                            "Poisson"             = paste(rep("b_",nreg),0:nreg,sep=""),
                            "Negative Binomial"   = c(paste(rep("b_",nreg),0:nreg,sep=""), "k"),
                            "Mixed Poisson"       = c(paste(rep("b_1",nreg),0:nreg,sep=""),paste(rep("b_2",nreg),0:nreg,sep=""), "p"),
                            "Generalized Poisson" = c(paste(rep("b_",nreg),0:nreg,sep=""), "a"),
                            "Binomial"            = c("n", paste(rep("b_",nreg),0:nreg,sep=""))
    )
  }

  # create names of the ARMA parameters
  if(nAR>0) ARNames = paste("AR_",1:ARMAModel[1], sep="")
  if(nMA>0) MANames = paste("MA_",1:ARMAModel[2], sep="")

  # put all the names together
  if(nAR>0 && nMA<1) parmnames = c(MargParmsNames, ARNames)
  if(nAR<1 && nMA>0) parmnames = c(MargParmsNames, MANames)
  if(nAR>0 && nMA>0) parmnames = c(MargParmsNames, ARNames, MANames)

  # add the parmnames on theta fix me: does this affect performance?
  if(!is.null(initialParam)) names(initialParam) = parmnames

  # create the constraints
  if(CountDist =="Poisson"){
    if(nreg==0){
      LB = c(0.001, rep(-Inf, sum(ARMAModel)))
      UB = rep(Inf, sum(ARMAModel)+1)
    }else{
      LB = rep(-Inf, sum(ARMAModel)+nreg+1)
      UB = rep(Inf, sum(ARMAModel)+nreg+1)
    }
  }

  if(CountDist == "Negative Binomial"){
    if(nreg==0){
      LB = c(0.01, 0.01, rep(-Inf, sum(ARMAModel)))
      UB = c(Inf, 0.99,   rep( Inf, sum(ARMAModel)))
    }else{
      LB = c(rep(-Inf, nreg+1), 0.001, rep(-Inf, sum(ARMAModel)))
      UB = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAModel)))
    }
  }


  if(CountDist == "Generalized Poisson"){
    if(nreg==0){
      LB = c(0.001, 0.001, rep(-Inf, sum(ARMAModel)))
      UB = c(Inf, Inf,     rep( Inf, sum(ARMAModel)))
    }else{
      LB = c(rep(-Inf, nreg+1), 0.001, rep(-Inf, sum(ARMAModel)))
      UB = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAModel)))
    }
  }


  out = list(
    mycdf           = mycdf,
    mypdf           = mypdf,
    myinvcdf        = myinvcdf,
    MargParmIndices = MargParmIndices,
    initialParam    = initialParam,
    TrueParam       = TrueParam,
    parmnames       = parmnames,
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
    ARMAModel       = ARMAModel,
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
    Task            = Task,
    OptMethod       = OptMethod,
    OutputType      = OutputType,
    ParamScheme     = ParamScheme
  )
  return(out)

}

# PF likelihood with resampling for AR(p) - written more concicely
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


  # Retrieve parameters ans save them in a list

  Parms = RetrieveParameters(theta,mod)

    # check for causality
  if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))

  # Initialize the negative log likelihood computation
  nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
                   - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1])))

  # Compute the theoretical covariance for the AR model for current estimate
  gt    = ARMAacf(ar = Parms$AR, ma = Parms$MA)[2:(max(mod$ARMAModel)+1)]

  # Compute the best linear predictor coefficients and errors using Durbin Levinson
  DL    = DLAcfToAR(gt)
  phit  = DL[,1]
  Rt    = sqrt(DL[,3])


  # allocate memory for particle weights and the latent Gaussian Series particles
  w     = matrix(0,mod$n, mod$ParticleNumber)
  Z     = matrix(0,mod$ARMAModel[1],mod$ParticleNumber)

  #======================   Start the SIS algorithm   ======================#
  # Initialize the weights
  w[1,] = rep(1,mod$ParticleNumber)

  # Compute the first integral limits Limit$ a and Limit$b
  Limit = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # Initialize particles from truncated normal distribution
  Z[1,] = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # =================== Loop over t ===================== #
  for (t in 2:mod$n){
    # Compute the latent Gaussian predictions Zhat_t using Durbin Levinson
    Zhat  = ComputeZhat(Z,phit,t)

    # Compute integral limits
    Limit = ComputeLimits(mod, Parms, t, Zhat, Rt)

    # Sample truncated normal particles
    Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b,t, Zhat, Rt)

    # update weights
    w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
      return(10^8)
    }

    # Resample the particles using common random numbers
    old_state1 = get_rand_state()
    Znew = ResampleParticles(mod, w, t, Znew)
    set_rand_state(old_state1)

    # Combine current particles, with particles from previous iterations
    if (mod$ARMAModel[1]>1){
      Z = rbind(Znew, Z[1:( min(t,mod$ARMAModel[1]) -1),])
    }else {
      Z[1,]=Znew
    }

    # update log-likelihood
    nloglik = nloglik - log(mean(w[t,]))
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
  if(mod$ARMAModel[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAModel[1])]

  MA = NULL
  if(mod$ARMAModel[2]>0) MA = theta[(mod$nMargParms+mod$ARMAModel[1]+1) : (mod$nMargParms + mod$ARMAModel[1] + mod$ARMAModel[2]) ]

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
  if(mod$ARMAModel[2]>0) Inn = matrix(0,N,mod$ARMAModel[2])   # I will save here the q many innovations (Z - Zhat) --see (5.3.9) BD book
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
    if(mod$ARMAModel[2]>1) Inn[,1:(mod$ARMAModel[2]-1)] = Inn[,2:(mod$ARMAModel[2])]

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
    Inn[,mod$ARMAModel[2]] = (znew-zhat)

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
      for(j in 1:min(t,mod$ARMAModel[2])){
        S = S-Theta[[t]][j]*Inn[,mod$ARMAModel[2]-j+1]
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
  if(mod$ARMAModel[1]>0 && mod$ARMAModel[2]==0) loglik = ParticleFilter_Res_AR(theta, mod)
  # Pure MA model or White noise
  if(mod$ARMAModel[1]==0&& mod$ARMAModel[2]>=0) loglik = ParticleFilter_Res_MA(theta, mod)

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

      if(mod$Task == 'Optimization'){
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


  if(mod$OutputType=="list"){
    #  save the results in a list
    ModelOutput = list(data.frame(matrix(rep(NA,nparms),nrow=1)),
                       data.frame(matrix(rep(NA,nparms),nrow=1)),
                       data.frame(matrix(rep(NA,4),     nrow=1)),
                       data.frame(matrix(rep(NA,3),     nrow=1))
    )
    # specify output list names
    names(ModelOutput)         = c("ParamEstimates", "StdErrors", "FitStatistics", "OptimOutput")
    ModelOutput$ParamEstimates = ParmEst
    ModelOutput$StdErrors      = se
    ModelOutput$FitStatistics  = c(loglik, Criteria)
    ModelOutput$OptimOutput    = c(convcode,kkt1,kkt2)
    ModelOutput$CountDist      = mod$CountDist
    ModelOutput$EstMethod      = mod$EstMethod
    ModelOutput$ARMAModel      = paste("ARMA(",mod$ARMAModel[1],",",mod$ARMAModel[2],")",sep="")
    ModelOutput$Task           = mod$Task

    # assign names to all output elements
    colnames(ModelOutput$ParamEstimates) = mod$parmnames
    colnames(ModelOutput$StdErrors)      = paste("se(", mod$parmnames,")", sep="")
    names(ModelOutput$FitStatistics)     = c("loglik", "AIC", "BIC", "AICc")
    names(ModelOutput$OptimOutput)       = c("ConvergeStatus", "kkt1", "kkt2")

  }else{
    ModelOutput  = data.frame(matrix(ncol = 4*mod$nparms+16, nrow = 1))

    # specify output names
    if(!is.null(mod$TrueParam)){
      colnames(ModelOutput) = c(
        'CountDist','ARMAModel', 'Regressor',
        paste("True_", mod$parmnames, sep=""), paste("InitialEstim_", mod$parmnames, sep=""),
        mod$parmnames, paste("se(", mod$parmnames,")", sep=""),
        'EstMethod', 'SampleSize', 'ParticleNumber', 'epsilon', 'OptMethod', 'ParamScheme',
        "loglik", "AIC", "BIC", "AICc", "ConvergeStatus", "kkt1", "kkt2")
    }
    colnames(ModelOutput) = c(
      'CountDist','ARMAModel', 'Regressor',
      paste("InitialEstim_", mod$parmnames, sep=""),
      mod$parmnames, paste("se(", mod$parmnames,")", sep=""),
      'EstMethod', 'SampleSize', 'ParticleNumber', 'epsilon', 'OptMethod',
      "loglik", "AIC", "BIC", "AICc", "ConvergeStatus", "kkt1", "kkt2")

    # Start Populating the output data frame
    ModelOutput$CountDist      = mod$CountDist
    ModelOutput$ARMAModel      = paste("ARMA(",mod$ARMAModel[1],",",mod$ARMAModel[2],")",sep="")
    ModelOutput$Regressor      = !is.null(mod$Regressor)

    offset = 4
    # true Parameters
    if(!is.null(mod$TrueParam)){
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
    if(!is.null(mod$TrueParam)) ModelOutput$ParamScheme    = mod$ParamScheme
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
  est  = rep(NA, mod$nMargParms+sum(mod$ARMAModel))
  #-----Poisson case
  if(mod$CountDist=="Poisson"){
    if(mod$nreg==0){
      est[1] = mean(mod$DependentVar)
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$nMargParms+mod$ARMAModel[1]):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }else{
      # GLM for the mean that depends on time
      # CHECK ME: If I fit a Poisson AR(3) in the the data example of the JASA paper, but the code below doesn't specify poisson family (it would pick up the default distribution that glm function has) then there will be a numerical error in the likelihood. Check it!
      glmPoisson            = glm(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)], family = "poisson")
      est[1:mod$nMargParms] = as.numeric(glmPoisson[1]$coef)

      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$ARMAModel[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
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
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$nMargParms+mod$ARMAModel[1]):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }else{
      # GLM for the mean that depends on time
      glmNegBin                 = glm.nb(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)])
      est[1:(mod$nMargParms-1)] = as.numeric(glmNegBin[1]$coef)
      # Mom on constant variance
      est[mod$nMargParms]       = NegBinMoM(mod$DependentVar,glmNegBin$fitted.values)
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$ARMAModel[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
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

      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$nMargParms+mod$ARMAModel[1]):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(mod$DependentVar~mod$Regressor)[1]$coef)
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$ARMAModel[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
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
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$nMargParms+mod$ARMAModel[1]):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(mod$DependentVar~mod$Regressor)[1]$coef)
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$ARMAModel[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }
  }


  # add the parmnames on theta fix me: does this affect performance?
  names(est) = mod$parmnames


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
sim_lgc = function(n, CountDist, MargParm, ARParm, MAParm, Regressor=NULL){


  # Generate latent Gaussian model
  z  =arima.sim(model = list( ar = ARParm, ma=MAParm  ), n = n)

  # select the correct count model
  if(is.null(Regressor)){
    myinvcdf = switch(CountDist,
                      "Poisson"             = qpois,
                      "Negative Binomial"   = function(x, theta){ qnbinom (x, theta[1], 1-theta[2]) },
                      "Mixed Poisson"       = function(x, theta){ qmixpois(x, theta[1], theta[2], theta[3])},
                      "Generalized Poisson" = function(x, theta){ qGpois  (x, theta[1], theta[2])},
                      "Binomial"            = qbinom)

    x = myinvcdf(pnorm(z), MargParm)
  }else{
    myinvcdf = switch(CountDist,
                      "Poisson"             = function(x, ConstMargParm, DynamMargParm){             qpois   (x, DynamMargParm)},
                      "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ qnbinom (x, ConstMargParm, 1-DynamMargParm)},
                      "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ qGpois  (x, ConstMargParm, DynamMargParm)}
    )

    # number of regressors assuming there is an intercept, fix me: check the case with not intercept
    nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)
    beta  = MargParm[1:(nreg+1)]
    m     = exp(Regressor%*%beta)

    if(CountDist == "Poisson" && nreg>0){
      ConstMargParm  = NULL
      DynamMargParm  = m
    }

    if(CountDist == "Negative Binomial" && nreg>0){
      ConstMargParm  = 1/MargParms[nreg+2]
      DynamMargParm  = MargParms[nreg+2]*m/(1+MargParms[nreg+2]*m)
    }

    if(CountDist == "Generalized Poisson" && nreg>0){
      ConstMargParm  = MargParms[nreg+2]
      DynamMargParm  = m
    }

    x = myinvcdf(pnorm(z), ConstMargParm, DynamMargParm)

  }

  return(x)
}



myppois = function(x, lambda){
  # compute poisson cdf as the ratio of an incomplete gamma function over the standard gamma function
  # I will also compute the derivative of the poisson cdf wrt lambda
  X  =c(lambda,x+1)
  v1 = gammainc(X)
  v2 = gamma(x+1)

  # straight forward formula from the definition of incomplete gamma integral
  v1_d = -lambda^x*exp(-lambda)
  v2_d = 0

  z  = v1/v2
  z_d = (v1_d*v2 - v2_d*v1)/v2^2
  return(c(z,z_d))
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


ComputeLimits = function(mod, Parms, t, Zhat, Rt){
  # a and b are the arguments in the two normal cdfs in the 4th line in equation (19) in JASA paper
  Lim = list()
  index = min(t, mod$ARMAModel[1])
  if(mod$nreg==0){
    Lim$a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,t(Parms$MargParms)),0,1)) - Zhat)/(Rt[index])
    Lim$b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],t(Parms$MargParms)),0,1)) - Zhat)/Rt[index]
  }else{
    Lim$a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,Parms$ConstMargParm, Parms$DynamMargParm[t]),0,1)) - Zhat)/Rt[index]
    Lim$b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],Parms$ConstMargParm, Parms$DynamMargParm[t]),0,1)) - Zhat)/Rt[index]
  }

  return(Lim)
}

SampleTruncNormParticles = function(mod, a, b, t, Zhat, Rt){
  # relation (21) in JASA paper and the inverse transform method
  # check me: this can be improved?
  index = min(t, mod$ARMAModel[1])
  z = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)*Rt[index] + Zhat
  return(z)
}

ComputeZhat = function(Z,phit,t){
  # Given particles Z and Durbin Levinson coefficients phit compute the best linear predictor Zhat_t
  # there is probably an easier way to do this but is ok for now
  if(t <= mod$ARMAModel[1]){
    if (t==2) {
      Zhat = Z[1,]*phit[1]
    } else{
      Zhat = colSums(Z[1:(t-1),]*phit[1:(t-1)])
    }
  }
  else{

    if(mod$ARMAModel[1]>1){# colsums doesnt work for 1-dimensional matrix
      Zhat = colSums(Z*phit)
    }else{
      Zhat =  Z*phit
    }

  }


  return(Zhat)
}

ComputeWeights = function(mod, a, b, t, PreviousWeights){
  # equation (21) in JASA paper
  # update weights
  if(t<=mod$ARMAModel[1]){
    NewWeights = PreviousWeights*(pnorm(b,0,1) - pnorm(a,0,1))
  }else{ # fix me: if I add the wgh[t-1,] below as I should the weights become small?
    NewWeights = (pnorm(b,0,1) - pnorm(a,0,1))
  }

  return(NewWeights)
}

ResampleParticles = function(mod, wgh, t, Znew){

  # relation (26) in JASA paper and following step
  # compute normalized weights
  wghn = wgh[t,]/sum(wgh[t,])

  # effective sample size
  ESS = 1/sum(wghn^2)

  if(ESS<mod$epsilon*mod$ParticleNumber){
    ind = rmultinom(1,mod$ParticleNumber,wghn)
    Znew = rep(Znew,ind)
  }
return(Znew)
}


RetrieveParameters = function(theta,mod){

  Parms =  vector(mode = "list", length = (mod$nparms+2))


  names(Parms) = c("MargParms", "ConstMargParm", "DynamMargParm", "AR", "MA")

  # marginal parameters
  Parms$MargParms      = theta[mod$MargParmIndices]

  # regressor parameters
  if(mod$nreg>0){
    beta  = Parms$MargParms[1:(mod$nreg+1)]
    m     = exp(mod$Regressor%*%beta)
  }

  # GLM type parameters
  if(mod$CountDist == "Negative Binomial" && mod$nreg>0){
    Parms$ConstMargParm  = 1/Parms$MargParms[mod$nreg+2]
    Parms$DynamMargParm  = Parms$MargParms[mod$nreg+2]*m/(1+Parms$MargParms[mod$nreg+2]*m)
  }

  if(mod$CountDist == "Generalized Poisson" && mod$nreg>0){
    Parms$ConstMargParm  = Parms$MargParms[mod$nreg+2]
    Parms$DynamMargParm  = m
  }

  if(mod$CountDist == "Poisson" && mod$nreg>0){
    Parms$ConstMargParm  = NULL
    Parms$DynamMargParm  = m
  }


  # Parms$AR = NULL
  if(mod$ARMAModel[1]>0) Parms$AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAModel[1])]

  # Parms$MA = NULL
  if(mod$ARMAModel[2]>0) Parms$MA = theta[(mod$nMargParms+mod$ARMAModel[1]+1) :
                                      (mod$nMargParms + mod$ARMAModel[1] + mod$ARMAModel[2]) ]


  return(Parms)
}



