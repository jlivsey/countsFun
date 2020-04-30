# Purpose: Calculate coverage of CI's calculated with Sandwich estimators vs
#           not using them.
#           Model is Poisson - AR(1)
library(countsFun)
library(numDeriv)

sim_sandy <- function(lam, phi, n = 100, nsim = 200){

  # simulate data
  x <- sim_pois_ar(n = n, phi = phi, lam = lam)

  # set initial values for optim
  inital.values <- c(jitter(lam), jitter(phi))

  # run optim
  out <- optim(par     = inital.values,
               fn      = GaussLogLik,
               data    = x,
               hessian = TRUE)

  # Calulate numerical Hessian
  h <- hessian(func = GaussLogLik, x = out$par, data = x)

  # Calculate numerical gradient
  g <- grad(func = GaussLogLik, x = out$par, data = x)

  # standard errors usual way using hessian
  SE.hess <- sqrt(diag(solve(h)))

  # standard erros with sandwhich method
  M <- (g) %*% t(g) # meat
  B <- h          # bread
  SE.sand <- solve(-M) %*% B %*% solve(-M)



}



