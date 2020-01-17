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

HermCoef_k <- function(lam, k, polys){
  #######################################################################
  # PURPOSE    Compute kth Hermite Coefficient. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   lam      Marginal parameter
  #   k        index of Hermite coefficient
  #   polys    Hermite Polynomial
  #
  # Output
  #   HC_k     kth hermite coeficient
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  # compute kth Hermite Polynomial
  her <- as.function(polys[[k+1]]) # polys[[k]] = H_{k-1}

  # truncation numbe: check me
  N <- which(round(ppois(1:1000, lam), 7) == 1)[1]

  # compute terms in the sum of relation (21) in
  terms <- exp(-qnorm(ppois(0:N, lam, lower.tail= TRUE))^2/2) *
    her(qnorm(ppois(0:N, lam, lower.tail = TRUE)))

  # take the sum of all terms
  HC_k <- sum(terms)/sqrt(2*pi)/factorial(k)
  return(HC_k)
}



HermCoef <- function(lam, polys, maxCoef = 20){
  #######################################################################
  # PURPOSE    Compute all Hermite Coefficients. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   lam      Marginal parameter
  #   polys    Hermite Polynomial
  #   maxCoef  number of coefficients to return. Default = 20
  #
  # Output
  #   HC       All Hermite coeficients
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  h = 0:maxCoef #check me
  HC = rep(NA, length(h)) # storage
  for(i in h) {
    HC[i+1] <- HermCoef_k(lam = lam, k = i, polys = polys)
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
  HC = HermCoef(lam,polys) # polys is Global

  # ARMA autocorrelation function
  ar.acf <- ARMAacf(ar = phi, lag.max = n)

  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  gamma_x = CountACVF(h = 0:(n-1), lam = lam, myacf = ar.acf, g = HC, 20)

  # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = toeplitz(gamma_x)
  return(GAMMA)
}



CountACVF_h = function(h, lam, myacf, g, max.terms = 20){
  #######################################################################
  # PURPOSE    Compute the autocovariance matrix of the count series.
  #            See relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  # INPUT
  #   h        acvf lag
  #   lam      Marginal Parameter
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

  k = 1:max.terms #check me
  gamma_x = sum(g[k]^2 *  factorial(k) * (myacf[h+1])^(k))
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
  GAMMA = CovarPoissonAR(n,lam,phi)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), lam)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}



FitGaussianLik = function(initialParam,x){
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
                        method = "BFGS")
  return(optim.output$par)

}








