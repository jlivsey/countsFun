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
  Z_D     = matrix(0,mod$ARMAModel[1],N)

  #======================   Start the SIS algorithm   ======================#
  # Initialize the weights and the latent Gaussian series particles
  wgh[1,] = rep(1,N)
  wgh_d[1,] = rep(1,N)

  if(mod$nreg==0){
    # compute count cdf and its derivative wrt the marginal paramterer in one step
    F1   = myppois(mod$DependentVar[1]-1,t(MargParms))
    v1   = F1[1]
    v1_d = F1[2]

    v2   = qnorm(v1,0,1)
    v2_d = 1/dnorm(qnorm(v1,0,1),0,1)*v1_d

    # fix me: there is a more efficient way to compute F2 from F1 and the pmf, but for now it
    # is more readable to do it  this way
    F2   = myppois(mod$DependentVar[1],t(MargParms))
    v3   = F2[1]
    v3_d = F2[2]

    v4   = qnorm(v3,0,1)
    v4_d = 1/dnorm(qnorm(v3,0,1),0,1)*v3_d

    a    = v5   = rep(v2, N)
    a_d  = v5_d = rep(v2_d,N)

    b    = v6   = rep(v4, N)
    b_d  = v6_d = rep(v4_d,N)
  }else{
    a    = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1,ConstMargParm, DynamMargParm[1]),0,1), N)
    b    = rep( qnorm(mod$mycdf(mod$DependentVar[1],ConstMargParm, DynamMargParm[1]),0,1), N)
  }

  v7      = pnorm(v5,0,1)
  v7_d    = dnorm(v5,0,1)*v5_d
  v8      = pnorm(v6,0,1)
  v8_d    = dnorm(v6,0,1)*v6_d

  v9      = runif(length(a),0,1)*(v8-v7)+v7
  v9_d    = runif(length(a),0,1)*(v8_d-v7_d)+v7_d

  Z[1,]   = v10   = qnorm(v9,0,1)
  Z_d[1,] = v10_d = 1/dnorm(v10,0,1)*v9_d


  # =================== Loop from 2 to AR order===================== #
  if (mod$ARMAModel[1]>=2){
    for (t in 2: (mod$ARMAModel[1])){
      # STEP 1 in SIS: Compute the latent Gaussian predictions Zhat using Durbin Levinson
      if (t==2) {
        Zhat   = v11   = Z[1:(t-1),]*phit[1:(t-1)]
        # Chcek me: Is this correct?
        Zhat_d = v11_d = Z_d[1:(t-1),]*phit[1:(t-1)]
      } else{
        Zhat   = v11   = colSums(Z[1:(t-1),]*phit[1:(t-1)])
        Zhat_d = v11_d = colSums(Z_d[1:(t-1),]*phit[1:(t-1)])
      }

      # STEP 2 is SIS: Update the latent Gaussian series Z and the importance weights w
      if(mod$nreg==0){
        # compute count cdf and its derivative wrt the marginal paramterer in one step
        F1   = myppois(mod$DependentVar[t]-1,t(MargParms))
        v13   = F1[1]
        v13_d = F1[2]

        v14   = qnorm(v13,0,1)
        v14_d = 1/dnorm(qnorm(v13,0,1),0,1)*v3_d


        F2   = myppois(mod$DependentVar[t],t(MargParms))
        v15   = F1[1]
        v15_d = F1[2]

        v16   = qnorm(v15,0,1)
        v16_d = 1/dnorm(v16,0,1)*v15_d

        # fix me I need to add Rt here - for now suppose it is constant
        a   = v17   = (v14 - v11)/Rt[t]
        a_d = v17_d = (v14_d - v11_d)/Rt[t]
        b   = v18   = (v16 - v11)/Rt[t]
        b_d = v18_d = (v16 - v11)/Rt[t]

      }else{
        a = (qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - Zhat)/Rt[t]
        b = (qnorm(mod$mycdf(mod$DependentVar[t],ConstMargParm, DynamMargParm[t]),0,1) - Zhat)/Rt[t]
      }

      v19   = pnorm(v17,0,1)
      v19_d = dnorm(v17,0,1)*v17_d
      v20   = pnorm(v18,0,1)
      v20_d = dnorm(v18,0,1)*v18_d

      v21   = runif(length(a),0,1)*(v20-v19)+v19
      v21_d = runif(length(a),0,1)*(v20_d-v19_d)+v19_d

      v22   = qnorm(v21,0,1)
      v22_d = 1/dnorm(v22,0,1)*v21_d

      v23   = v22*Rt[t] + v11
      v23_d = v22_d*Rt[t] + v11_d

      Z[1:t,]   = rbind(v23, Z[1:(t-1),])
      Z_d[1:t,] = rbind(v23_d, Z_d[1:(t-1),])

      wgh[t,] = wgh[t-1,]*(v20 - v19)
      wgh_d[t,] = wgh_d[t-1,]*(v20 - v19) + wgh[t-1,]*(v20_d - v19_d)

      # update likelihood
      nloglik   = nloglik   - log(mean(wgh[t,]))
      nloglik_d = nloglik_d - 1/(mean(wgh[t,])) + mean(wgh_d[t,])
    }
  }
  # =================== Loop from AR order + 1  to T ===================== #
  # From p to T1 I don't need to estimate phi anymore
  for (t in (mod$ARMAModel[1]+1):T1){

    # STEP 1 in SIS: Compute the latent Gaussian predictions Zhat using Durbin Levinson
    if(mod$ARMAModel[1]>1){# colsums doesnt work for 1-dimensional matrix
      Zhat   = v24   = colSums(Z*phit)
      Zhat_d = v24_d = colSums(Z_d*phit)
    }else{
      Zhat   = v24   =  Z*phit
      Zhat_d = v24_d =  Z_d*phit
    }

    # STEP 2 is SISR: Update the latent Gaussian series Z
    if(mod$nreg==0){
      a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,t(MargParms)),0,1)) - Zhat)/Rt[mod$ARMAModel[1]]
      b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],t(MargParms)),0,1)) - Zhat)/Rt[mod$ARMAModel[1]]
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
