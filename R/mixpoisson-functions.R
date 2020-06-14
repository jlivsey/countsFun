# #---------inverse cdf of mixpoisson series---------#
# qmixpois = function(y, p, lam1, lam2){
#   yl = length(y)
#   x = rep(0,yl)
#   for (n in 1:yl){
#     while(pmixpois(x[n], p, lam1, lam2) <= y[n]){ # R qpois would use <y; this choice makes the function right-continuous; this does not really matter for our model
#       x[n] = x[n]+1
#     }
#   }
#   return(x)
# }
#
#
# #---------cdf of mixpoisson series---------#
# pmixpois = function(x, p, lam1, lam2,lower.tail){
#   y = p*ppois(x,lam1,lower.tail = lower.tail) + (1-p)*ppois(x,lam2,lower.tail=lower.tail)
#   return(y)
# }
#
# #---------generate mixpoisson series---------#
# rmixpois = function(n, p, lam1, lam2){
#   u = runif(n)
#   x = rpois(n,lam1)*(u<=p) + rpois(n,lam2)*(u>p)
#   return(x)
# }
#
#
# #---------simulate negbin series with prob of success p---------#
# sim_mixedPoisson = function(n, ARMAmodel,p,lam1,lam2){
#   #====================================================================================#
#   # PURPOSE       Simulate NegBin series with ARMA structure. See relation (1)
#   #               in https://arxiv.org/pdf/1811.00203.pdf
#   #
#   # INPUT
#   #   n           series length
#   #   ARMAmodel   list with ARMA parameters
#   #   r,p         Marginal Parameters
#   #
#   # Output
#   #   X           Neg bin series
#   #
#   # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date          April 2020
#   # Version       3.6.3
#   #====================================================================================#
#
#   z = arima.sim(model = list(ar=ARMAmodel[[1]], ma=ARMAmodel[[2]]), n = n); z = z/sd(z)
#   # The third argument of qnbinom is the prob of failure
#   X = qmixpois(pnorm(z), p, lam1, lam2)
#   return(X)
# }
#
#
# #---------Hermitte Coefficients for all k---------#
# HermCoefMixPois <- function(p, lam1,lam2, N, nHC){
#   #====================================================================================#
#   # PURPOSE    Compute all Hermite Coefficients. See relation (21) in
#   #            https://arxiv.org/pdf/1811.00203.pdf
#   #
#   # INPUT
#   #   r,p      Marginal parameters
#   #   maxCoef  number of coefficients to return. Default = 20
#   #   N        truncation of relation (21)
#   #   nHC      number of Hermitte coefficients
#   #
#   # Output
#   #   HC       All Hermite coeficients
#   #
#   # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date       April 2020
#   # Version    3.6.3
#   #====================================================================================#
#
#   h = 1:nHC #check me
#   HC = rep(NA, length(h)) # storage
#   for(i in h) {
#     HC[i] <- HermCoefMixPois_k(p,lam1,lam2 , k = i, N)
#   }
#
#   return(HC)
#
# }
#
#
# #---------Hermitte Coefficients for one k---------#
# HermCoefMixPois_k <- function(r, p, k, N){
#   #====================================================================================#
#   # PURPOSE    Compute kth Hermite Coefficient. See relation (21) in
#   #            https://arxiv.org/pdf/1811.00203.pdf
#   #
#   # INPUT
#   #   r,p      Marginal parameters
#   #   k        index of Hermite coefficient
#   #
#   # Output
#   #   HC_k     kth hermite coeficient
#   #
#   # Authors    Stefanos Kechagias, James Livsey
#   # Date       January 2020
#   # Version    3.6.1
#   #====================================================================================#
#
#   # function for kth Hermite Polynomial
#   her <- function(x){
#     evalHermPolynomial(k-1, x)
#   }
#
#   # # truncation numbe: check me
#   # N <- which(round(pnbinom(1:1000, r,p), 7) == 1)[1]
#   # if(is.na(N)){
#   #   N=1000
#   # }
#   # N = max(sapply(unique(p),function(x)which(round(pnbinom(1:1000, r,x), 7) == 1)[1]))
#
#   # compute terms in the sum of relation (21) in
#   terms <- exp((-qnorm(   pmixpois(0:max(N), p, lam1,lam2, lower.tail = TRUE)   )^2)/2) *
#              her(qnorm(   pmixpois(0:max(N), p, lam1,lam2, lower.tail = TRUE)   ))
#
#   terms[is.nan(terms)]=0
#
#   # take the sum of all terms
#   HC_k <- sum(terms) / (sqrt(2*pi) *  factorial(k))
#   return(HC_k)
# }
#
#
# #---------Covariance matrix---------#
# CovarMixPois = function(n, p, lam1,lam2, AR, MA, N, nHC){
#   #====================================================================================#
#   # PURPOSE     Compute the covariance matrix of a NegBin series.
#   #
#   # INPUT
#   #   r,p       Marginal parameters
#   #   AR,MA     ARMA parameters
#   #   n         size of the matrix
#   #   N         truncation for relation (21)
#   #   nHC       number of Hermitte coefficents to be computed
#   #
#   # Output
#   #   GAMMA     covariance matrix ofcount series
#   #
#   # Authors     Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date        April 2020
#   # Version     3.6.3
#   #====================================================================================#
#
#   # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
#   HC = HermCoefMixPois(p,lam1,lam2,N, nHC)
#
#   # ARMA autocorrelation function
#   if(!length(AR)){arma.acf <- ARMAacf(ma = MA, lag.max = n)}
#   if(!length(MA)){arma.acf <- ARMAacf(ar = AR, lag.max = n)}
#   if(length(AR) & length(MA)){arma.acf <- ARMAacf(ar = AR, ma = MA, lag.max = n)}
#
#   # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
#   gamma_x = CountACVF(h = 0:(n-1), myacf = arma.acf, g = HC)
#
#   # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
#   GAMMA = toeplitz(gamma_x)
#   return(GAMMA)
# }
#
#
# #---------Gaussian Likelihood function---------#
# GaussLogLikMP = function(theta, data, ARMAorder, MaxCdf, nHC){
#   #====================================================================================#
#   # PURPOSE      Compute Gaussian log-likelihood for NegBin series
#   #
#   # INPUT
#   #   theta      parameter vector containing the marginal and ARMA parameters
#   #   data       count series
#   #   ARMAorder  ordeer of ARMA model
#   #   MaxCdf     cdf will be computed up to this number (for light tails cdf=1 fast)
#   #   nHC        number of HC to be computed
#   #
#   #
#   # Output
#   #   loglik     Gaussian log-likelihood
#   #
#   # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date         April 2020
#   # Version      3.6.3
#   #====================================================================================#
#
#   # retrieve parameters and sample size
#   p = theta[1]
#   lam1 = theta[2]
#   lam2 = theta[3]
#   nparms  = length(theta)
#   nMargParms = nparms - sum(ARMAorder)
#   if(ARMAorder[1]>0){
#     AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
#   }else{
#     AR = NULL
#   }
#
#   if(ARMAorder[2]>0){
#     MA = theta[ (nMargParms+ARMAorder[1]+1) : length(theta)]
#   }else{
#     MA = NULL
#   }
#   n = length(data)
#
#   # compute truncation of relation (21)
#   N <- which(round( pmixpois(1:MaxCdf, p,lam1,lam2)   , 7) == 1)[1]
#
#
#   if(length(N)==0 |is.na(N) ){
#     N =MaxCdf
#   }
#
#   #Select the mean value used to demean--sample or true?
#   MeanValue = p*lam1 + (1-p)*lam2
#
#   # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
#   GAMMA = CovarMixPois(n, p, lam1, lam2, AR, MA, N, nHC)
#
#   # Compute the logdet and the quadratic part
#   logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)
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
# #---------wrapper to fit Gaussian Likelihood function---------#
# FitGaussianLikMP = function(x0, X, LB, UB, ARMAorder, MaxCdf, nHC){
#   #====================================================================================#
#   # PURPOSE       Fit the Gaussian log-likelihood for Mixed Poisson series
#   #
#   # INPUT
#   #   x0          initial parameters
#   #   X           count series
#   #   LB          parameter lower bounds
#   #   UB          parameter upper bounds
#   #   ARMAorder   order of the udnerlying ARMA model
#   #   MaxCdf      cdf will be computed up to this number (for light tails cdf=1 fast)
#   #   nHC         number of HC to be computed
#   #
#   # OUTPUT
#   #   All         parameter estimates, standard errors, likelihood value
#   #
#   # NOTES         I may comment out se in cases where maximum is achieved at the boundary
#   #
#   # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date          April 2020
#   # Version       3.6.3
#   #====================================================================================#
#   optim.output <- optim(par       = x0,
#                         fn        = GaussLogLikMP,
#                         data      = X,
#                         ARMAorder = ARMAorder,
#                         MaxCdf    = MaxCdf,
#                         nHC       = nHC,
#                         method    = "L-BFGS-B",
#                         hessian   = TRUE,
#                         lower     = LB,
#                         upper     = UB
#   )
#
#   nparms  = length(x0)
#   ParmEst = matrix(0,nrow=1,ncol=nparms)
#   se      = matrix(NA,nrow=1,ncol=nparms)
#   loglik  = rep(0,1)
#
#   # save estimates, loglik and standard errors
#   ParmEst[,1:nparms]   = optim.output$par
#   loglik               = optim.output$value
#   #se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))
#
#   All      = cbind(ParmEst, se, loglik)
#   return(All)
#
# }
#
#
