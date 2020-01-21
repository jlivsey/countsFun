# Generate AR series
sim_pois_ar = function(n, phi, lam){
  #######################################################################
  # PURPOSE   Simulate Poisson series with AR structure. See relation (1)
  #           in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   n        series length
  #   phi      AR parameter
  #   lam      Marginal Parameter
  #
  # Output
  #   x        Poisson series
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qpois(pnorm(z), lam)
  return(x)
}


evalPolynomial_scalarX <- function(coefs, x){
  #######################################################################
  # PURPOSE    compute arbitrary polynomial given coefficients and input
  #
  # INPUT
  #   coefs    vector of coefficients
  #   x        input
  #
  # Output
  #   value   returns f(x) = coefs[1] * x^0 +
  #                          coefs[2] * x^1 +
  #                                     ... +
  #                          coefs[n] * x^(n-1)
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################
  if(length(coefs) == 1){
      return(coefs) # handle scalar case
  }
  n <- length(coefs) - 1
  out <- sum(coefs * x^(0:n))
  return(out)
}
evalPolynomial <- Vectorize(evalPolynomial_scalarX, vectorize.args = "x")

#######################################################################
# PURPOSE    Evaluate k^th Hermite Polynomial at scaler input.
#           See relation (7) https://arxiv.org/pdf/1811.00203.pdf
#
# INPUT
#   k        index of Hermite coefficient
#   x        input to Hermite Polynomial
#
#
# Output
#   value   returns H_k(x). See relation (7)
#           in https://arxiv.org/pdf/1811.00203.pdf
#
# Authors    Stefanos Kechagias, James Livsey
# Date       January 2020
# Version    3.6.1
#######################################################################
evalHermPolynomial <- function(k, x){
  if(k < 0){
    return(0)
  }
  coefs <- HermPolyCoefs[[k+1]]
  out <- evalPolynomial(coefs, x)
  return(out)
}


HermCoef_k <- function(lam, k){
  #######################################################################
  # PURPOSE    Compute kth Hermite Coefficient. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   lam      Marginal parameter
  #   k        index of Hermite coefficient
  #
  # Output
  #   HC_k     kth hermite coeficient
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  # function for (k-1)st Hermite Polynomial
  her <- function(x){
    evalHermPolynomial(k-1, x)
  }

  # truncation numbe: check me
  N <- which(round(ppois(1:1000, lam), 7) == 1)[1]

  # compute terms in the sum of relation (21) in
  terms <- exp((-qnorm(ppois(0:N, lam, lower.tail= TRUE))^2)/2) *
    her(qnorm(ppois(0:N, lam, lower.tail = TRUE)))

  # take the sum of all terms
  HC_k <- sum(terms) / (sqrt(2*pi) *  factorial(k))
  return(HC_k)
}



HermCoef <- function(lam){
  #######################################################################
  # PURPOSE    Compute all Hermite Coefficients. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   lam      Marginal parameter
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
    HC[i] <- HermCoef_k(lam = lam, k = i)
  }

  return(HC)

}


CovarPoissonAR = function(n,lam,phi){
  #######################################################################
  # PURPOSE    Compute the covariance matrix of a Poisson AR series.
  #
  # INPUT
  #   lam      Marginal parameter
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
  HC = HermCoef(lam)

  # ARMA autocorrelation function
  ar.acf <- ARMAacf(ar = phi, lag.max = n)

  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  gamma_x = CountACVF(h = 0:(n-1), myacf = ar.acf, g = HC)

  # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = toeplitz(gamma_x)
  return(GAMMA)
}



CountACVF_h = function(h, myacf, g){
  #######################################################################
  # PURPOSE    Compute the autocovariance matrix of the count series.
  #            See relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  # INPUT
  #   h        acvf lag
  #   myacf    autocorrelation of Gaussian series: rho_z
  #   g        Hermitte Coefficients
  # Output
  #   gamma_x  count acvf
  #
  # Notes:     See also CountACVF function--vectorized version of CountACVF_h
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  k = length(g) #check me
  gamma_x = sum(g^2 *  factorial(1:k) * (myacf[h+1])^(1:k))
  return(gamma_x)

  }

CountACVF <- Vectorize(CountACVF_h, vectorize.args = "h")


EvalInvQuadForm = function(A, v, DataMean ){
  #######################################################################
  # Evaluate quadrative form v`*inv(A)*v where
  # A is symmetric positive definite
  # Want   QuadForm = v` * inv(A) * v = v` * w  where A*w = v
  # ==> (U`U)*w = v             Cholesky decomp of A
  # ==>  First solve U` z = v   for z,
  # then solve   U w = z   for w */
  #######################################################################

  U = chol(A)
  z  =  forwardsolve(t(U), v - DataMean)
  w  = backsolve(U,z)
  QuadForm = t(v-DataMean)%*%w

  logdetA = 2*sum(log(diag(U)))

  logLikComponents = c(logdetA, QuadForm)
  return(logLikComponents)
}


# ---- likelihood function ----
GaussLogLik = function(theta, data){
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
  # Authors    Stefanos Kechagias, James Livsey
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

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), lam)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}



FitGaussianLik = function(initialParam, x){
  #######################################################################
  # PURPOSE    Fit the Gaussian log-likelihood for Poisson AR series
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
                        fn = GaussLogLik,
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



