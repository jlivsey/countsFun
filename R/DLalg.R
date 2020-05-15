#' Univariate Durbin-Levinson algoithm.
#' This implimentation takes a mean value and subtracts it during calculation
#' of one-step-ahead predictors and residuals
#'
#' @param x univariate time series
#' @param g ACVF function where g[1] = gamma(0)
#' @param lam mean
#'
#' @return list with lower triangular matrix of predicition coefficients, mses,
#' one-step-ahead predictions and residuals
#' @export
#'
#' @examples
DLalg <- function(x, g, lam){
  # parameters known from inputs
  n <- length(x)

  # Iniitalize storage
  phi <- matrix(NA, n, n)
  v <- rep(NA, n)

  # inital values
  phi[1, 1] <- g[2]/g[1]
  v[1] <- g[1]

  # main loop
  for(m in 2:n){
    phi[m, m] <- (g[m+1] - sum(phi[m-1, 1:(m-1)] * g[(m-1+1):2])) / v[m-1]
    phi[m, 1:(m-1)] <- phi[m-1, 1:(m-1)] - phi[m, m] * phi[m-1, (m-1):1]
    v[m] <- v[m-1] * (1 - phi[m-1, m-1]^2)
  }

  # Calculate xhat's
  xhat <- rep(NA, n+1) #storage
  xhat[1] <- 0
  for(m in 1:n){
    xhat[m+1] <- sum(phi[m, 1:m] * (x[m:1] - lam))
  }

  # residuals
  e <- (x - lam) - xhat[1:n]

  return(list(
    phi  = phi,
    v    = v,
    xhat = xhat,
    e    = e
  ))

}
