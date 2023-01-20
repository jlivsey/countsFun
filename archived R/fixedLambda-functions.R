# CovarPoissonAR_fixedLamba = function(n, HC, phi){
#   #######################################################################
#   # PURPOSE    Compute the covariance matrix of a Poisson AR series.
#   #
#   # INPUT
#   #   lam      Marginal parameter
#   #   phi      AR parameter
#   #   n        size of the matrix
#   #
#   # Output
#   #   GAMMA    covariance matrix ofcount series
#   #
#   # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date       January 2020
#   # Version    3.6.1
#   #######################################################################
#
#   # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
#
#
#   # ARMA autocorrelation function
#   ar.acf <- ARMAacf(ar = phi, lag.max = n)
#
#   # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
#   gamma_x = CountACVF(h = 0:(n-1), myacf = ar.acf, g = HC)
#
#   # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
#   GAMMA = toeplitz(gamma_x)
#   return(GAMMA)
# }
#
#
# # ---- likelihood function ----
# GaussLogLik_fixedLambda = function(theta, data, lam, HC){
#   #######################################################################
#   # PURPOSE    Compute Gaussian log-likelihood for Poisson AR series
#   #
#   # INPUT
#   #   theta    parameter vector containing the marginal and AR parameters
#   #   data     count series
#   #
#   # Output
#   #   loglik   Gaussian log-likelihood
#   #
#   # Authors    Stefanos Kechagias, James Livsey
#   # Date       January 2020
#   # Version    3.6.1
#   #######################################################################
#
#   # retrieve parameters and sample size
#   phi = theta[1]
#   n = length(data)
#
#   # assign large likelihood value if not causal and if lambda outside range
#   if(any(abs( polyroot(c(1, -phi))  ) < 1) || (lam < 0 || lam > 100) ){
#     return(NA) #check me
#   }
#
#   # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
#   GAMMA = CovarPoissonAR_fixedLamba(n, HC, phi)
#
#   # Compute the logdet and the quadratic part
#   logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), lam)
#
#   # final loglikelihood value
#   out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]
#
#   # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
#   # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
#   return(out)
# }
#
#
#
# FitGaussianLik_fixedLambda = function(initialParam, x, lam, HC){
#   #######################################################################
#   # PURPOSE    Fit the Gaussian log-likelihood for Poisson AR series
#   #
#   # INPUT
#   #   initialParam       parameter vector containing the marginal and AR parameters
#   #   x                  count series
#   #
#   # Output
#   #   optim.output$par   parameter estimates
#   #
#   # Authors    Stefanos Kechagias, James Livsey
#   # Date       January 2020
#   # Version    3.6.1
#   #######################################################################
#   optim.output <- optim(par = initialParam,
#                         fn = GaussLogLik_fixedLambda,
#                         data = x,
#                         HC = HC,
#                         lam = lam,
#                         method = "BFGS",
#                         hessian=TRUE)
#
#   nparms = length(initialParam)
#   ParmEst = matrix(0,nrow=1,ncol=nparms)
#   se =  matrix(NA,nrow=1,ncol=nparms)
#   loglik = rep(0,1)
#
#   # save estimates, loglik and standard errors
#   ParmEst[,1:nparms]   = optim.output$par
#   loglik               = optim.output$value
#   se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))
#
#   All      = cbind(ParmEst, se, loglik)
#   return(All)
#
# }
#
#
