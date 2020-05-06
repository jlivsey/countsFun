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


###
# fit independent mixtures of Poissons to quakes dataset
# from Zucchini & MacDonald Hidden Markov Models for Timeseries exercise 1.3
# see http://petewerner.blogspot.com/2014/12/fitting-mixture-of-independent-poisson.html
##

# #quakes data
# xf <- c(13, 14, 8, 10, 16, 26, 32, 27, 18, 32, 36, 24, 22, 23, 22, 18, 25, 21, 21, 14, 8, 11, 14, 23, 18, 17, 19, 20, 22, 19, 13, 26, 13, 14, 22, 24, 21, 22, 26, 21, 23, 24, 27, 41, 31, 27, 35, 26, 28, 36, 39, 21, 17, 22, 17, 19, 15, 34, 10, 15, 22, 18, 15, 20, 15, 22, 19, 16, 30, 27, 29, 23, 20, 16, 21, 21, 25, 16, 18, 15, 18, 14, 10, 15, 8, 15, 6, 11, 8, 7, 18, 16, 13, 12, 13, 20, 15, 16, 12, 18, 15, 16, 13, 15, 16, 11, 11)



mixedPois_MLE <- function(x, inital.value){

  xf <- x

  udist <- function(n) rep(1/n, n)

  #natural to working
  n2wp <- function(p) {
    m <- length(p)
    log(p[2:m]/(1 - sum(p[2:m])))
  }

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

  #initial estimates and probabilities for 2, 3 and 4 distributions
  #the lambda values I just guess, and use an uniform distribution
  #for the initial mixing distribution.
  # pv <- c(log(c(10, 20)), n2wp(udist(2)))
  # pv <- c(log(c(10, 15, 20)), n2wp(udist(3)))
  # pv <- c(log(c(5, 15, 20, 30)), n2wp(udist(4)))

  pv <- inital.value


  #number of distributions to fit
  # m <- (length(pv) + 1)/2

  m <- 2

  #fit using nlm
  fv <- suppressWarnings( nlm(of, pv, m, xf, print.level=0) )
  rv <- fv$est

  #lambda estimates
  lam.est <- exp(rv[1:m])
  #mixing distribution
  prob.est <- w2np(rv[(m+1):(2*m-1)])

  return(c(sort(lam.est), sort(prob.est)[1]))
}
