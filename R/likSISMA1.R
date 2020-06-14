likSISMA1 = function(theta, data,ARMAorder, CountDist){
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1,
                           "Negative Binomial"   = 1:2,
                           "Mixed Poisson"       = 1:3,
                           "Generalized Poisson" = 1:2,
                           "Binomial"            = 1:2)

  # retrieve marginal cdf
  cdf = switch(CountDist,
               "Poisson"                       = ppois,
               "Negative Binomial"             = function(x, theta){ pnbinom(q = x, size = theta[1], prob = 1-theta[2]) },
               "Mixed Poisson"                 = pMixedPoisson,
               "Generalized Poisson"           = pGenPoisson,
               "Binomial"                      = pbinom
  )

  # retrieve marginal pdf
  pdf = switch(CountDist,
               "Poisson"                       = dpois,
               "Negative Binomial"             = function(x, theta){ dnbinom(x, size = theta[1], prob = 1-theta[2]) },
               "Mixed Poisson"                 = dMixedPoisson,
               "Generalized Poisson"           = dGenPoisson,
               "Binomial"                      = dbinom
  )

  #set.seed(1)
  theta1.idx = MargParmIndices
  theta2.idx = ARMAorder[1]




    #set.seed(1)
    theta1 = theta[theta1.idx]
    n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
    theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
    tht = theta[theta2.idx]
    xt = data
    T1 = length(xt)
    N = 5 # number of particles
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

    wprev = rep(1,N)
    wgh[,1] = wprev

    for (t in 2:T1)
    {
      rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
      rt = sqrt(rt0/(1+tht^2))
      a = (qnorm(cdf(xt[t]-1,theta1),0,1) - zhat)/rt
      b = (qnorm(cdf(xt[t],theta1),0,1) - zhat)/rt
      err = z.rest(a,b)
      znew = zhat + rt*err
      zhat = tht*(znew-zhat)/rt0
      prt[,t] = zhat

      wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
      wprev = wgh[,t]
    }

    lik = pdf(xt[1],theta1)*mean(wgh[,T1])
    nloglik = (-2)*log(lik)

    out = if (is.na(nloglik)) Inf else nloglik
    return(out)

}
