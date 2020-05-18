# .........................................................
# Mixed Poisson distribution functions
# .........................................................


#' Generate random mixed-Poisson draws
#'
#' @param n number of observations to generate
#' @param p mixing probability
#' @param lam1 mean of first component
#' @param lam2 mean of second component
#'
#' @return vector of randomly generated values
#' @export
#'


rmixpois = function(n, p, lam1, lam2){
  u = runif(n)
  x = rpois(n,lam1)*(u<=p) + rpois(n,lam2)*(u>p)
  return(x)
}

#' density of mixed-Poisson distribution
#'
#' @param x evaluation point (or vector)
#' @param p mixing probability
#' @param lam1 mean of first component
#' @param lam2 mean of second component
#'
#' @return density at x with specified parameters
#' @export
#'

dmixpois = function(x, p, lam1, lam2){
  y = p*dpois(x,lam1) + (1-p)*dpois(x,lam2)
  return(y)
}

#' CDF of mixed-Poisson distribution
#'
#' @param x evaluation point (or vector)
#' @param p mixing probability
#' @param lam1 mean of first component
#' @param lam2 mean of second component
#'
#' @return CDF at x with specified parameters
#' @export
#'

pmixpois = function(x, p, lam1, lam2){
  y = p*ppois(x,lam1) + (1-p)*ppois(x,lam2)
  return(y)
}

#' Inverse CDF/Quantile function for mixed-Poisson distribution
#'
#' @param y evaluation point (or vector)
#' @param p mixing probability
#' @param lam1 mean of first component
#' @param lam2 mean of second component
#'
#' @return quantile function at y with specified parameters
#' @export
#'

qmixpois = function(y, p, lam1, lam2){
  yl = length(y)
  x = rep(0,yl)
  for (n in 1:yl){
    while(pmixpois(x[n], p, lam1, lam2) <= y[n]){ # R qpois would use <y; this choice makes the function right-continuous; this does not really matter for our model
      x[n] = x[n]+1
    }
  }
  return(x)
}


# .........................................................
# Generate mixed Poisson-AR1 model
# Inputs:
#  T: length of the series
#  p: mixture porbability
#  lam1, lam2: Poisson lambda parameters
#  phi: AR(1) phi parameter
# .........................................................

#' Generate simulated mixed-Poisson-AR1 model
#'
#' @param Tt length of the series
#' @param phi AR(1) parameter
#' @param p mixture probability
#' @param lam1 mean of first component
#' @param lam2 mean of second component
#'
#' @return simulated series of length Tt
#' @export
#'

sim_mixpois_ar1 <- function(Tt,phi,p,lam1,lam2){
  zt <- rep(0,Tt)
  zt[1] <- rnorm(1,0,1)
  for (t in 2:Tt){
    zt[t] = phi*zt[t-1] + sqrt(1-phi^2)*rnorm(1,0,1)
  }
  xt <- qmixpois(pnorm(zt,0,1),p,lam1,lam2)
  return(xt)
}



# kth hermitte coefficient
HermCoefMixedPois_k <- function(lam1, lam2, prob, k){
  #######################################################################
  # PURPOSE    Compute kth Hermite Coefficient. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   lam1      Marginal parameter - first mean
  #   lam2      Marginal parameter - second mean
  #   prob      Marginal parameter - mixing probability
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
  N <- which(round(pmixpois(0:1000, prob, lam1, lam2), 7) == 1)[1]

  # compute terms in the sum of relation (21) in
  terms <- exp((-qnorm(pmixpois(0:N, prob, lam1, lam2))^2)/2) *
    her(qnorm(pmixpois(0:N, prob, lam1, lam2)))

  # take the sum of all terms
  HC_k <- sum(terms) / (sqrt(2*pi) *  factorial(k))
  return(HC_k)
}

# all hermitte coefficients
HermCoefMixedPois <- function(lam1, lam2, prob, maxCoef = 20){
  #######################################################################
  # PURPOSE    Compute all Hermite Coefficients. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   lam1     Marginal parameter - first mean
  #   lam2     Marginal parameter - second mean
  #   prob     Marginal parameter - mixing probability
  #   maxCoef  number of coefficients to return. Default = 20
  #
  # Output
  #   HC       All Hermite coeficients
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #######################################################################

  h = 1:maxCoef #check me
  HC = rep(NA, length(h)) # storage
  for(i in h) {
    HC[i] <- HermCoefMixedPois_k(lam1, lam2, prob, k = i)
  }

  return(HC)

}


# fit independent mixtures of Poissons with MLE
# Adaptation of code from from Zucchini & MacDonald Hidden Markov Models
# for Timeseries exercise 1.3
mixedPois_MLE <- function(x, inital.value){

  xf <- x

    #working to natural
  w2np <- function(lp) {
    rv <- exp(lp)/(1 + sum(exp(lp)))
    c(1 - sum(rv), rv)
  }

  #optimisation function
  of <- function(pv, m, x) {
    if (m == 1)
      return(-sum(dpois(x, exp(pv), log=TRUE)))
    #convert working parameters to natural paramters
    pr <- exp(pv[1:m])
    probs <- w2np(pv[(m+1):(2*m - 1)])
    #calculate -ve log likelihood
    -sum(log(outer(x, pr, dpois) %*% probs))
  }

  #number of distributions to fit
  m <- 2

  #fit using nlm
  fv <- suppressWarnings( nlm(of, inital.value, m, xf, print.level=0) )
  rv <- fv$est

  #lambda estimates
  lam.est <- exp(rv[1:m])
  #mixing distribution
  prob.est <- w2np(rv[(m+1):(2*m-1)])

  # return both mean parameters and smaller probability
  return(c(sort(lam.est), sort(prob.est)[1]))
}







































#---------Gaussian Likelihood function---------#
GaussLogLik_mixedPois = function(theta, data, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE      Compute Gaussian log-likelihood for mixed-Poisson series
  #
  # INPUT
  #   theta      parameter vector containing the marginal and ARMA parameters
  #   data       count series
  #   ARMAorder  order of ARMA model
  #   MaxCdf     cd fwill be computed up to this number (for light tails cdf=1 fast)
  #   nHC        number of HC to be computed
  #
  #
  # Output
  #   loglik     Gaussian log-likelihood
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         April 2020
  # Version      3.6.3
  #====================================================================================#

  # retrieve parameters and sample size
  r = theta[1]
  p = theta[2]
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
  n = length(data)

  # compute truncation of relation (21)
  N <- which(round(pnbinom(1:MaxCdf, r,1-p), 7) == 1)[1]
  # if(length(N)==0 |is.na(N) ){
  #   cat(sprintf("The max cdf value is %f and N=%f", max(round(pnbinom(1:MaxCdf, r,1-p), 7)),N))
  #   stop("Haven't reached upper limit for cdf")
  # }
  if(length(N)==0 |is.na(N) ){
    N =MaxCdf
  }

  #Select the mean value used to demean--sample or true?
  MeanValue = r*p/(1-p)


  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = CovarNegBin(n, r, p, AR, MA, N, nHC)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}


#---------wrapper to fit Gaussian Likelihood function---------#
FitGaussianLik_mixedPois = function(x0, X, LB, UB, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE       Fit the Gaussian log-likelihood for NegBin series
  #
  # INPUT
  #   x0          initial parameters
  #   X           count series
  #   LB          parameter lower bounds
  #   UB          parameter upper bounds
  #   ARMAorder   order of the udnerlying ARMA model
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
  optim.output <- optim(par       = x0,
                        fn        = GaussLogLikNB,
                        data      = X,
                        ARMAorder = ARMAorder,
                        MaxCdf    = MaxCdf,
                        nHC       = nHC,
                        method    = "L-BFGS-B",
                        hessian   = TRUE,
                        lower     = LB,
                        upper     = UB
  )

  nparms  = length(x0)
  ParmEst = matrix(0,nrow=1,ncol=nparms)
  se      = matrix(NA,nrow=1,ncol=nparms)
  loglik  = rep(0,1)

  # save estimates, loglik and standard errors
  ParmEst[,1:nparms]   = optim.output$par
  loglik               = optim.output$value
  #se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))

  All      = cbind(ParmEst, se, loglik)
  return(All)

}

