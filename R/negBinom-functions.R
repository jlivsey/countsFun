
# ---- likelihood function ----
GaussLogLikNB = function(theta, data){
  #######################################################################
  # PURPOSE    Compute Gaussian log-likelihood for NegBin AR series
  #
  # INPUT
  #   theta    parameter vector containing the marginal and AR parameters
  #   data     count series
  #
  # Output
  #   loglik   Gaussian log-likelihood
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  # retrieve parameters and sample size
  r = theta[1]
  p = theta[2]
  phi = theta[-c(1,2)]
  n = length(data)

  #Select the mean value used to demean--sample or true?
  MeanValue = r*p/(1-p)

  # assign large likelihood value if not causal or if meanm outside range
  if(any(abs( polyroot(c(1, -phi))  ) < 1) || (MeanValue < 0 || MeanValue > 100) ){
    return(NA) #check me
  }

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = CovarNegBinAR(n, r, p, phi)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}



CovarNegBinAR = function(n,r, p, phi){
  #######################################################################
  # PURPOSE    Compute the covariance matrix of a NegBin AR series.
  #
  # INPUT
  #   r,p      Marginal parameters
  #   phi      AR parameter
  #   n        size of the matrix
  #
  # Output
  #   GAMMA    covariance matrix ofcount series
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
  HC = HermCoefNegBin(r,p)

  # ARMA autocorrelation function
  ar.acf <- ARMAacf(ar = phi, lag.max = n)

  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  gamma_x = CountACVF(h = 0:(n-1), myacf = ar.acf, g = HC)

  # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = toeplitz(gamma_x)
  return(GAMMA)
}


HermCoefNegBin <- function(r,p){
  #######################################################################
  # PURPOSE    Compute all Hermite Coefficients. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   r,p      Marginal parameters
  #   maxCoef  number of coefficients to return. Default = 20
  #
  # Output
  #   HC       All Hermite coeficients
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  h = 1:20 #check me
  HC = rep(NA, length(h)) # storage
  for(i in h) {
    HC[i] <- HermCoefNegBin_k(r, p , k = i)
  }

  return(HC)

}


HermCoefNegBin_k <- function(r,p, k){
  #######################################################################
  # PURPOSE    Compute kth Hermite Coefficient. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   r,p      Marginal parameters
  #   k        index of Hermite coefficient
  #
  # Output
  #   HC_k     kth hermite coeficient
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  # function for kth Hermite Polynomial
  her <- function(x){
    evalHermPolynomial(k-1, x)
  }

  # truncation numbe: check me
  N <- which(round(pnbinom(1:1000, r,p), 7) == 1)[1]

  # compute terms in the sum of relation (21) in
  terms <- exp((-qnorm(pnbinom(0:N, r,p, lower.tail= TRUE))^2)/2) *
    her(qnorm(pnbinom(0:N, r,p, lower.tail = TRUE)))

  # take the sum of all terms
  HC_k <- sum(terms) / (sqrt(2*pi) *  factorial(k))
  return(HC_k)
}


FitGaussianLikNB = function(initialParam, x){
  #######################################################################
  # PURPOSE    Fit the Gaussian log-likelihood for NegBin AR series
  #
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
                        fn = GaussLogLikNB,
                        data = x,
                        method = "BFGS",
                        hessian=TRUE)

  nparms = length(initialParam)
  ParmEst = matrix(0,nrow=1,ncol=nparms)
  se =  matrix(NA,nrow=1,ncol=nparms)
  loglik = rep(0,1)

  # save estimates, loglik and standard errors
  ParmEst[,1:nparms]   = optim.output$par
  loglik               = optim.output$value
  se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))

  All      = cbind(ParmEst, se, loglik)
  return(All)

}


sim_negbin_ar = function(n, phi, r,p){
  #######################################################################
  # PURPOSE   Simulate Poisson series with AR structure. See relation (1)
  #           in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   n        series length
  #   phi      AR parameter
  #   r,p      Marginal Parameters
  #
  # Output
  #   x        Poisson series
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qnbinom(pnorm(z), r,p)
  return(x)
}
