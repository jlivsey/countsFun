# Univariate Durbin-Levinson algoithm.

# notes:
#         1. I'm assuming variance gamma(0) = 1
#         2. gamma vec being passed does not include variance so g[1] = gamma(1)

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
    v[m] <- v[m-1] * (1 - phi[m, m]^2)
  }

  # Calculate xhat's
  xhat <- rep(NA, n+1) #storage
  xhat[1] <- 0
  for(m in 1:n){
    xhat[m+1] <- sum(phi[m, 1:m] * x[m:1])
  }

  # residuals
  e <- x - xhat[1:n]

  return(list(
    phi  = phi,
    v    = v,
    xhat = xhat,
    e    = e
  ))

}
