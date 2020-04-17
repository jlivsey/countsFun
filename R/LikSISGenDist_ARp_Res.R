#' Title
#'
#' @param theta
#' @param data
#'
#' @return
#' @export
#'
#'
LikSISGenDist_ARp_Res = function(theta, data, ParticleNumber, CountDist){
  ##########################################################################
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
  #
  # OUTPUT:
  #    loglik:           approximate log-likelihood
  #
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    November 2019
  ##########################################################################

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments

  # retrieve indices of marginal distribution parameters
  MargParmIndices = switch(CountDist,
                           "Poisson" = 1,
                 "Negative Binomial" = 1:2,
                     "Mixed Poisson" = 1:3,
               "Generalized Poisson" = 1:2,
                          "Binomial" = 1:2)

  # retrieve marginal cdf
  mycdf = switch(CountDist,
                 "Poisson" = ppois,
       "Negative Binomial" = pnbinom ,
           "Mixed Poisson" = pMixedPoisson,
     "Generalized Poisson" = pGenPoisson,
                "Binomial" = pbinom
  )

  # retrieve marginal pdf
  mypdf = switch(CountDist,
                 "Poisson" = dpois,
           "Mixed Poisson" = dMixedPoisson,
                "Binomial" = dbinom,
       "Negative Binomial" = dnbinom,
     "Generalized Poisson" = dGenPoisson
  )

  # retrieve marginal distribution parameters
  MargParms = theta[MargParmIndices]
  nMargParms = length(MargParms) # num param in MargParms

  # retrieve AR order--FIX ME--add argument and retrieve it from there
  ARorder = length(theta)-nMargParms

  # number of dependence parameters
  ARParms.idx = (nMargParms + 1):(nMargParms + ARorder)

  # retrieve AR parameters
  phi = theta[ARParms.idx]

  if (prod(abs(polyroot(c(1,-phi))) > 1)){ # check if the ar model is causal

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
    phit = TacvfAR(phi)[2]/TacvfAR(phi)[1]
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
        Gt = toeplitz(TacvfAR(phi)[1:t])
        gt = TacvfAR(phi)[2:(t+1)]
        phit = as.numeric(solve(Gt) %*% gt)
        rt =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/TacvfAR(phi)[1]))

      }
    }


    # From p to T1 I dont need to estimate phi anymore
    for (t in (ARorder+1):T1){
      # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
      if(ARorder>1){# colsums doesnt work for 1-dimensional matrix
        ZpreviousTimesPhi = colSums(ZprevAll*phi)
      }else{
        ZpreviousTimesPhi=ZprevAll*phi
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
      epsilon=0.5
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
