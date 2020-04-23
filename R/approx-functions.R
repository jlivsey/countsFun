# ---- likelihood function ----
GaussLikApprox = function(theta, data){
  #######################################################################
  # PURPOSE    Compute Gaussian log-likelihood for Poisson AR series
  #
  # INPUT
  #   theta    parameter vector containing the marginal and AR parameters
  #   data     count series
  #
  # Output
  #   loglik   Gaussian log-likelihood
  #
  # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  # retrieve parameters and sample size
  lam = theta[1]
  phi = theta[-1]
  n = length(data)

  # assign large likelihood value if not causal and if lambda outside range
  if(any(abs( polyroot(c(1, -phi))  ) < 1) || (lam < 0 || lam > 100) ){
    return(NA) #check me
  }

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = CovarPoissonAR(n, lam, phi)


  x <- as.numeric(data)
  lamvec <- rep(lam, n)
  U <- chol(GAMMA)
  vec <- backsolve(U, x - lamvec, transpose = TRUE)

  quadform <- sum(vec^2)

  det <- prod(diag(U))^2

  out <- quadform * det^(1/n)

  return(out)
}

FitGaussianLikApprox = function(initialParam, x){
  #######################################################################
  # PURPOSE    Fit the Gaussian log-likelihood for Poisson AR series
  #
  # test
  # INPUT
  #   initialParam       parameter vector containing the marginal and AR parameters
  #   x                  count series
  #
  # Output
  #   optim.output$par   parameter estimates
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################
  optim.output <- optim(par = initialParam,
                        fn = GaussLikApprox,
                        data = x,
                        method = "BFGS",
                        hessian=FALSE)

  nparms = length(initialParam)
  ParmEst = matrix(0,nrow=1,ncol=nparms)
  se =  matrix(NA,nrow=1,ncol=nparms)
  loglik = rep(0,1)

  # save estimates, loglik and standard errors
  ParmEst[,1:nparms]   = optim.output$par
  loglik               = optim.output$value
  # se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))

  All      = cbind(ParmEst, se, loglik)
  return(All)

}

