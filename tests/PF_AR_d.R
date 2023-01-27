# PF likelihood with resampling for AR(p)
PF_AR_d = function(theta, mod){
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
  if(mod$ARMAModel[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAModel[1])]

  MA = NULL
  if(mod$ARMAModel[2]>0) MA = theta[(mod$nMargParms+mod$ARMAModel[1]+1) :
                                      (mod$nMargParms + mod$ARMAModel[1] + mod$ARMAModel[2]) ]

  # check for causality
  if( CheckStability(AR,MA) ) return(10^(8))


  # sample size and number of particles
  T1      = length(mod$DependentVar)
  N       = mod$ParticleNumber

  # Initialize the negative log likelihood computation
  nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],MargParms)),
                   - log(mod$mypdf(mod$DependentVar[1], ConstMargParm, DynamMargParm[1])))

  # Compute the theoretical covariance for the AR model for current estimate
  gt      = ARMAacf(ar = AR, ma = MA)[2:(max(mod$ARMAModel)+1)]

  # Compute the best linear predictor coefficients and errors using Durbin Levinson
  DL      = DLAcfToAR(gt)
  phit    = DL[,1]
  Rt      = sqrt(DL[,3])


  # allocate memory for particle weights and the latent Gaussian Series particles
  wgh     = matrix(0,T1,N)
  Z       = matrix(0,mod$ARMAModel[1],N)

  #======================   Start the SIS algorithm   ======================#
  # Initialize the weights and the latent Gaussian series particles
  wgh[1,] = rep(1,N)

  if(mod$nreg==0){
    v1 = mod$mycdf(mod$DependentVar[1]-1,t(MargParms))
    v2 = qnorm(v1,0,1)
    v3 = mod$mycdf(mod$DependentVar[1],t(MargParms))
    v4 = qnorm(v3,0,1)
    a  = rep(v2, N)
    b  = rep(v4, N)
  }else{
    a       = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1,ConstMargParm, DynamMargParm[1]),0,1), N)
    b       = rep( qnorm(mod$mycdf(mod$DependentVar[1],ConstMargParm, DynamMargParm[1]),0,1), N)
  }

  Z[1,]   = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

  # =================== Loop from 2 to AR order===================== #
  if (mod$ARMAModel[1]>=2){
    for (t in 2: (mod$ARMAModel[1])){
      # STEP 1 in SIS: Compute the latent Gaussian predictions Zhat using Durbin Levinson
      if (t==2) {
        Zhat = Z[1:(t-1),]*phit[1:(t-1)]
      } else{
        Zhat = colSums(Z[1:(t-1),]*phit[1:(t-1)])
      }

      # STEP 2 is SIS: Update the latent Gaussian series Z and the importance weights w
      if(mod$nreg==0){
        v5 = mod$mycdf(mod$DependentVar[t]-1,t(MargParms))
        v6 = qnorm(v5,0,1)
        v7 = mod$mycdf(mod$DependentVar[t],t(MargParms))
        v8 = qnorm(v7,0,1)
        a = (v6 - Zhat)/Rt[t]
        b = (v8 - Zhat)/Rt[t]

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
  for (t in (mod$ARMAModel[1]+1):T1){

    # STEP 1 in SIS: Compute the latent Gaussian predictions Zhat using Durbin Levinson
    if(mod$ARMAModel[1]>1){# colsums doesnt work for 1-dimensional matrix
      Zhat = colSums(Z*phit)
    }else{
      Zhat =  Z*phit
    }

    # STEP 2 is SISR: Update the latent Gaussian series Z
    if(mod$nreg==0){
      v9 = mod$mycdf(mod$DependentVar[t]-1,t(MargParms))
      v10 = qnorm(v9,0,1)
      v11 = mod$mycdf(mod$DependentVar[t],t(MargParms))
      v12 = qnorm(v11,0,1)

      a = as.numeric(v10 - Zhat)/Rt[mod$ARMAModel[1]]
      b = as.numeric(v12 - Zhat)/Rt[mod$ARMAModel[1]]
    }else{
      a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1)) - Zhat)/Rt[mod$ARMAModel[1]]
      b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],ConstMargParm, DynamMargParm[t]),0,1)) - Zhat)/Rt[mod$ARMAModel[1]]
    }

    Znew = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)*Rt[mod$ARMAModel[1]] + Zhat

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
    if (mod$ARMAModel[1]>1){
      Z = rbind(Znew, Z[1:(mod$ARMAModel[1]-1),])
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
