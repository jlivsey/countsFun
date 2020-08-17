# PF likelihood with resampling for AR(p)
ParticleFilter_Res_AR = function(theta, xt, Regressor, mod){
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

  # retrieve marginal parameters
  MargParms        = theta[mod$MargParmIndices]

  # retrieve regressor parameters
  if(mod$nreg>0){
    beta  = MargParms[1:(mod$nreg+1)]
    m     = exp(Regressor%*%beta)
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


    T1 = length(xt)
    N = mod$ParticleNumber          # number of particles

    wgh = matrix(0,T1,N)        # to collect all particle weights

    # allocate memory for zprev
    ZprevAll = matrix(0,mod$ARMAorder[1],N)

    if(mod$nreg==0){
      # Compute integral limits
      a = rep( qnorm(mod$mycdf(xt[1]-1,t(MargParms)),0,1), N)
      b = rep( qnorm(mod$mycdf(xt[1],t(MargParms)),0,1), N)
    }else{
      a = rep( qnorm(mod$mycdf(xt[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
      b = rep( qnorm(mod$mycdf(xt[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
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
    if (mod$ARMAorder[1]>=2){
      for (t in 2:mod$ARMAorder[1]){

        # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
        if (t==2) {
          zhat = ZprevAll[1:(t-1),]*phit
        } else{
          zhat = colSums(ZprevAll[1:(t-1),]*phit)
        }

        # Recompute integral limits
        if(mod$nreg==0){
          a = (qnorm(mod$mycdf(xt[t]-1,t(MargParms)),0,1) - zhat)/rt
          b = (qnorm(mod$mycdf(xt[t],t(MargParms)),0,1) - zhat)/rt
        }else{
          a = (qnorm(mod$mycdf(xt[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
          b = (qnorm(mod$mycdf(xt[t],ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
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
    for (t in (mod$ARMAorder[1]+1):T1){

      # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
      if(mod$ARMAorder[1]>1){# colsums doesnt work for 1-dimensional matrix
        zhat = colSums(ZprevAll*AR)
      }else{
        zhat=ZprevAll*AR
      }

      # compute limits of truncated normal distribution
      if(mod$nreg==0){
        a = as.numeric(qnorm(mod$mycdf(xt[t]-1,MargParms),0,1) - zhat)/rt
        b = as.numeric(qnorm(mod$mycdf(xt[t],  MargParms),0,1) - zhat)/rt
      }else{
        a = as.numeric(qnorm(mod$mycdf(xt[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
        b = as.numeric(qnorm(mod$mycdf(xt[t],  ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
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
      if(ESS<mod$epsilon*N){
        ind = rmultinom(1,N,wghn)
        # sample particles
        znew = rep(znew,ind)

        # use low variance resampling
        #znew = lowVarianceRS(znew, wghn, N)
      }
      set_rand_state(old_state1)


      # save particles
      if (mod$ARMAorder[1]>1){
        ZprevAll = rbind(znew, ZprevAll[1:(mod$ARMAorder[1]-1),])
      }else {
        ZprevAll[1,]=znew
      }
      # update likelihood
      nloglik = nloglik - log(mean(wgh[t,]))
    }

    # likelihood approximation
    if(mod$nreg==0){
      nloglik = nloglik - log(mod$mypdf(xt[1],MargParms))
    }else{
      nloglik = nloglik - log(mod$mypdf(xt[1], ConstMargParm, DynamMargParm[1]))
    }

    # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
    nloglik = nloglik
    #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  return(nloglik)
}


# PF likelihood with resampling for MA(q)
ParticleFilter_Res_MA = function(theta, xt, Regressor, mod){
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
    m     = exp(Regressor%*%beta)
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


  T1 = length(xt)
  N = mod$ParticleNumber          # number of particles


    # allocate matrix to collect all particle weights
    wgh = matrix(0,length(data),N)

    # Compute integral limits
    if(mod$nreg==0){
      a = rep( qnorm(mod$mycdf(data[1]-1,t(MargParms)),0,1), N)
      b = rep( qnorm(mod$mycdf(data[1],t(MargParms)),0,1), N)
    }else{
      a = rep( qnorm(mod$mycdf(data[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
      b = rep( qnorm(mod$mycdf(data[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
    }

    # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
    zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)


    # run innovations Algorithm for MA models that are not WN
    if(mod$ARMAorder[2]>0) Inn = matrix(0,N,mod$ARMAorder[2])       # I will save here the q many innovations (Z - Zhat) --see (5.3.9) BD book
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
      if(mod$ARMAorder[2]>1) Inn[,1:(mod$ARMAorder[2]-1)] = Inn[,2:(mod$ARMAorder[2])]

      # compute limits of truncated normal distribution
      if(mod$nreg==0){
        a = as.numeric(qnorm(mod$mycdf(data[t]-1,MargParms),0,1) - zhat)/vt
        b = as.numeric(qnorm(mod$mycdf(data[t],MargParms),0,1) -   zhat)/vt
      }else{
        a = as.numeric(qnorm(mod$mycdf(data[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/vt
        b = as.numeric(qnorm(mod$mycdf(data[t],ConstMargParm, DynamMargParm[t]),0,1) -   zhat)/vt
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
      nloglik = nloglik - log(mod$mypdf(data[1],MargParms))
    }else{
      nloglik = nloglik - log(mod$mypdf(data[1], ConstMargParm, DynamMargParm[1]))
    }

    # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
    # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))


    if (nloglik==Inf | is.na(nloglik)){
      nloglik = 10^8
    }


  return(nloglik)
}


# PF likelihood with resampling
ParticleFilter_Res = function(theta, data, Regressor, mod){
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
  if(mod$ARMAorder[1]>0 && mod$ARMAorder[2]==0) loglik = ParticleFilter_Res_AR(theta, data, Regressor, mod)
  # Pure MA model or White noise
  if(mod$ARMAorder[1]==0&& mod$ARMAorder[2]>=0) loglik = ParticleFilter_Res_MA(theta, data, Regressor, mod)
  return(loglik)
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



# Optimization wrapper to fit PF likelihood with resamplinbg
FitMultiplePF_Res = function(theta, data, Regressor, mod, OptMethod){
  #====================================================================================#
  # PURPOSE       Fit the Particle Filter log-likelihood with resampling.
  #               This function maximizes the PF likelihood, nfit manys times for nparts
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
  nparts = length(mod$ParticleNumber)
  nparms = length(theta)
  nfit   = 1
  n      = length(data)

  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst  = matrix(0,nrow=nfit*nparts,ncol=nparms)
  se       =  matrix(NA,nrow=nfit*nparts,ncol=nparms)
  loglik   = rep(0,nfit*nparts)
  convcode = rep(0,nfit*nparts)
  kkt1     = rep(0,nfit*nparts)
  kkt2     = rep(0,nfit*nparts)


  # Each realization will be fitted nfit*nparts many times
  for (j in 1:nfit){
    set.seed(j)
    # for each fit repeat for different number of particles
    for (k in 1:nparts){
      # number of particles to be used
      ParticleNumber = mod$ParticleNumber[k]

      # run optimization for our model --no ARMA model allowed
      optim.output <- optimx(par            = theta,
                             fn             = ParticleFilter_Res,
                             data           = data,
                             Regressor      = Regressor,
                             mod            = mod,
                             lower          = mod$LB,
                             upper          = mod$UB,
                             hessian        = TRUE,
                             method         = OptMethod)



      # save estimates, loglik value and diagonal hessian
      ParmEst[nfit*(k-1)+j,]  = as.numeric(optim.output[1:nparms])
      loglik[nfit*(k-1) +j]   = optim.output$value
      convcode[nfit*(k-1) +j] = optim.output$convcode
      kkt1[nfit*(k-1) +j]     = optim.output$kkt1
      kkt2[nfit*(k-1) +j]     = optim.output$kkt2


      # compute Hessian
      H = gHgen(par            = ParmEst[nfit*(k-1)+j,],
                fn             = ParticleFilter_Res,
                data           = data,
                Regressor      = Regressor,
                mod            = mod)

      # save standard errors from Hessian
      if(H$hessOK && det(H$Hn)>10^(-8)){
        se[nfit*(k-1)+j,]   = sqrt(abs(diag(solve(H$Hn))))
      }else{
        se[nfit*(k-1)+j,] = rep(NA, nparms)
      }

    }
  }

  # Compute model selection criteria (assuming one fit)
  Criteria = ComputeCriteria(loglik, nparms, n, mod$ParticleNumber)


  # get the names of the final output
  parmnames = colnames(optim.output)
  mynames = c(parmnames[1:nparms],paste("se", parmnames[1:nparms], sep="_"), "loglik", "AIC", "BIC","AICc", "status", "kkt1", "kkt2")


  All = matrix(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)
  colnames(All) = mynames

  return(All)
}


