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
    her(qnorm(pmixpois(0:1000, prob, lam1, lam2)))

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
