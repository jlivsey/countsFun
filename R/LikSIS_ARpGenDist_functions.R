# Add some functions that I will need in particle filtering approximation of
# likelihood. See file LikSIS_ARpGenDist.R


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

# Generate AR series
sim_pois_ar = function(n, phi, lam){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qpois(pnorm(z), lam)
  return(x)
}

sim_MixPois_ar <- function(n,phi,lam1,lam2, prob){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x  = qMixedPoisson(pnorm(z),lam1,lam2, prob)
  return(x)
}


# mixed poisson cdf
pMixedPoisson = function(q, lam1, lam2, prob){
  # theta[1] is the mixing probability
  # theta[2], theta[3] are the lambda parameters
  prob*ppois(q, lam1) + (1-prob)*ppois(q, lam2)
}

dMixedPoisson = function(x, lam1, lam2, prob){
  prob*dpois(x, lam1) + (1-prob)*dpois(x, lam2)
}

qMixedPoisson = function(y, lam1, lam2, prob){
  yl = length(y)
  x = rep(0,yl)
  for (n in 1:yl){
    while(pMixedPoisson(x[n], lam1, lam2, prob) < y[n]){ # R qpois would use <y; this choice makes the function right-continuous; this does not really matter for our model
      x[n] = x[n]+1
    }
  }
  return(x)
}


# Generalized Poisson cdf, pdf
pGenPoisson = function(q, theta, lam){ cdf.vec <- rep(-99,length(q))
for (i in 1:length(q)) {
  if (q[i]>=0){
    cdf.vec[i] <- sum( exp(-(theta+(0:q[i])*lam))*theta*(theta+(0:q[i])*lam)^((0:q[i])-1)/factorial((0:q[i])))
  }else{cdf.vec[i] <-0}
}
return(cdf.vec)}

dGenPoisson = function(x, theta, lam){ exp(-(theta+x*lam))*theta*(theta+x*lam)^(x-1)/factorial(x) }


#---------wrapper to fit PF likelihood---------#
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



#---------------------------------------------------------------------------------------------#


# PF likelihood with resampling
ParticleFilterRes = function(theta, data, ARMAorder, ParticleNumber, CountDist, epsilon){
  #--------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling
  #           to approximate the likelihood of the
  #           a specified count time series model with an underlying AR(p)
  #           dependence structure.
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
  # DATE:    November 2019
  #--------------------------------------------------------------------------#

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

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
                 "Poisson"             = ppois,
                 "Negative Binomial"   = function(x, theta){ pnbinom(q = x, size = theta[1], prob = 1-theta[2])},
                 "Mixed Poisson"       = pMixedPoisson,
                 "Generalized Poisson" = pGenPoisson,
                 "Binomial"            = pbinom
  )

  # retrieve marginal pdf
  mypdf = switch(CountDist,
                 "Poisson"             = dpois,
                 "Negative Binomial"   = function(x, theta){ dnbinom(x, size = theta[1], prob = 1-theta[2]) },
                 "Mixed Poisson"       = dMixedPoisson,
                 "Generalized Poisson" = dGenPoisson,
                 "Binomial"            = dbinom
  )

  # retrieve marginal distribution parameters
  MargParms  = theta[MargParmIndices]
  nMargParms = length(MargParms)
  nparms     = length(theta)


  # retrieve ARMA parameters
  if(ARMAorder[1]>0){
    AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
  }else{
    AR = NULL
  }

  if(ARMAorder[2]>0){
    MA = theta[ (length(theta) - ARMAorder[2]) : length(theta)]
  }else{
    MA = NULL
  }

  # retrieve AR order
  ARorder = ARMAorder[1]

  if (prod(abs(polyroot(c(1,-AR))) > 1)){ # check if the ar model is causal

    xt = data
    T1 = length(xt)
    N = ParticleNumber          # number of particles
    prt = matrix(0,N,T1)        # to collect all particles
    wgh = matrix(0,T1,N)        # to collect all particle weights

    # allocate memory for zprev
    ZprevAll = matrix(0,ARorder,N)

    # Compute integral limits
    a = rep( qnorm(mycdf(xt[1]-1,t(MargParms)),0,1), N)
    b = rep( qnorm(mycdf(xt[1],t(MargParms)),0,1), N)

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
    #t0 = proc.time()
    # First p steps:

    if (ARorder>=2){
      for (t in 2:ARorder){

        # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
        if (t==2) {
          ZpreviousTimesPhi = ZprevAll[1:(t-1),]*phit
        } else{
          ZpreviousTimesPhi = colSums(ZprevAll[1:(t-1),]*phit)
        }

        # Recompute integral limits
        a = (qnorm(mycdf(xt[t]-1,MargParms),0,1) - ZpreviousTimesPhi)/rt
        b = (qnorm(mycdf(xt[t],MargParms),0,1) - ZpreviousTimesPhi)/rt

        # compute random errors from truncated normal
        err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

        # compute the new Z and add it to the previous ones
        znew = rbind(ZpreviousTimesPhi + rt*err, ZprevAll[1:(t-1),])
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
    for (t in (ARorder+1):T1){
      # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
      if(ARorder>1){# colsums doesnt work for 1-dimensional matrix
        ZpreviousTimesPhi = colSums(ZprevAll*AR)
      }else{
        ZpreviousTimesPhi=ZprevAll*AR
      }
      # compute limits of truncated normal distribution
      a = as.numeric(qnorm(mycdf(xt[t]-1,MargParms),0,1) - ZpreviousTimesPhi)/rt
      b = as.numeric(qnorm(mycdf(xt[t],MargParms),0,1) -   ZpreviousTimesPhi)/rt

      # draw errors from truncated normal
      err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

      # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
      znew = ZpreviousTimesPhi + rt*err

      # Resampling Step--here the function differs from LikSISGenDist_ARp

      # compute unnormalized weights
      wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)

      # break if I got NA weight
      if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
        nloglik = NA
        break
      }

      # normalized weights
      wghn = wgh[t,]/sum(wgh[t,])

      old_state1 <- get_rand_state()

      # sample indices from multinomial distribution-see Step 4 of SISR in paper
      ESS = sum(1/wghn^2)
      if(ESS<epsilon*N){
        ind = rmultinom(1,N,wghn)
        # sample particles
        znew = rep(znew,ind)

        # use low variance resampling
        #znew = lowVarianceRS(znew, wghn, N)
      }
      set_rand_state(old_state1)


      # save particles
      if (ARorder>1){
        ZprevAll = rbind(znew, ZprevAll[1:(ARorder-1),])
      }else {
        ZprevAll[1,]=znew
      }
      # update likelihood
      nloglik = nloglik - log(mean(wgh[t,]))
    }

    # likelihood approximation
    nloglik = nloglik - log(mypdf(xt[1],MargParms))


    # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
    nloglik = nloglik
    #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

    out =nloglik

  }else{
    out = NA # for noncasusal AR
  }

  return(out)
}

# PF likelihood with resampling for MA(1)
ParticleFilterMA1_Res = function(theta, data, ARMAorder, ParticleNumber, CountDist, epsilon){
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
print(theta)
old_state <- get_rand_state()
on.exit(set_rand_state(old_state))
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
        zhat = tht*(znew-zhat)/rt0

        # use low variance resampling
        #znew = lowVarianceRS(znew, wghn, N)
      }
      set_rand_state(old_state1)

      # update likelihood
      nloglik = nloglik - log(mean(wgh[t,]))
    }

    # likelihood approximation
    nloglik = nloglik - log(mypdf(xt[1],MargParms))


    # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
    nloglik = nloglik
    #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

    out =nloglik

    if (out==Inf | is.na(out)){
      out = 10^8
    }

  return(out)
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

    # break if I got NA weight
    if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
      nloglik = 10^6
      break
    }

    # update likelihood
    nloglik = nloglik - log(mean(wgh[t,]))
  }

  # likelihood approximation
  nloglik = nloglik - log(mypdf(xt[1],MargParms))


  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  nloglik = nloglik
  #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  out =nloglik


  return(out)
}



#---------new wrapper to fit PF likelihood---------#
FitMultiplePFNew = function(x0, X, CountDist, Particles, LB, UB, ARMAorder, epsilon, UseDEOptim){
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


  # Each realization will be fitted nfit*nparts many times
  for (j in 1:nfit){
    set.seed(j)
    # for each fit repeat for different number of particles
    for (k in 1:nparts){
      # number of particles to be used
      ParticleNumber = Particles[k]
      if(!UseDEOptim){
        # run optimization for our model
        optim.output <- optim(par            = x0,
                              fn             = ParticleFilterMA1_Res,
                              data           = X,
                              ARMAorder      = ARMAorder,
                              ParticleNumber = ParticleNumber,
                              CountDist      = CountDist,
                              epsilon        = epsilon,
                              lower          = LB,
                              upper          = UB,
                              hessian        = TRUE,
                              method         = "L-BFGS-B")
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


      # save estimates, loglik value and diagonal hessian
      ParmEst[nfit*(k-1)+j,]  = optim.output$par
      loglik[nfit*(k-1) +j]   = optim.output$value
      se[nfit*(k-1)+j,]       = sqrt(abs(diag(solve(optim.output$hessian))))
    }
  }

  All = cbind(ParmEst, se, loglik)
  return(All)
}

z.rest = function(a,b){
  # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
  qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
}


likSISR = function(theta, data){
  cdf = function(x, theta1){ pnbinom(q = x, size = theta1[1], prob = theta1[2]) }
  pdf = function(x, theta1){ dnbinom(x = x, size = theta1[1], prob = theta1[2]) }
  #set.seed(1)
  theta1.idx = 1:2
  theta2.idx = 3
  theta1 = theta[theta1.idx]
  n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
  theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
  tht = theta[theta2.idx]
  xt = data
  T1 = length(xt)
  N = 1000 # number of particles
  prt = matrix(0,N,T1) # to collect all particles
  wgh = matrix(0,N,T1) # to collect all particle weights

  a = qnorm(cdf(xt[1]-1,theta1),0,1)
  b = qnorm(cdf(xt[1],theta1),0,1)
  a = rep(a,N)
  b = rep(b,N)
  zprev = z.rest(a,b)
  rt0 = 1+tht^2
  zhat = tht*zprev/rt0
  prt[,1] = zhat


  nloglik <- 0
  for (t in 2:T1)
  {
    rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
    rt = sqrt(rt0/(1+tht^2))
    a = (qnorm(cdf(xt[t]-1,theta1),0,1) - zhat)/rt
    b = (qnorm(cdf(xt[t],theta1),0,1) - zhat)/rt
    err = z.rest(a,b)
    znew = zhat + rt*err
    wgh <- pnorm(b,0,1) - pnorm(a,0,1)
    if (any(is.na(wgh))) # see apf code below for explanation
    {
      #nloglik <- NaN
      nloglik <- Inf
      break
    }
    if (sum(wgh)==0)
    {
      #nloglik <- NaN
      nloglik <- Inf
      break
    }
    wghn <- wgh/sum(wgh)
    ind <- rmultinom(1, N, wghn)
    znew <- rep(znew,ind)

    zhat = tht*(znew-zhat)/rt0
    prt[,t] = zhat

    nloglik <- nloglik -2*log(mean(wgh))
  }

  nloglik <- nloglik - 2*log(pdf(xt[1],theta1))


  out = if (is.na(nloglik)) Inf else nloglik
  return(out)
}











