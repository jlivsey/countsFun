# Purpose: Calculate coverage of CI's calculated with Sandwich estimators vs
#           not using them.
#           Model is Poisson - AR(1)
library(countsFun)
library(numDeriv)
library(sandwich)


sim_sandy <- function(lam, phi, n = 100, nsim = 200){

  P <- matrix(NA, nsim, 2) # parameter estimates
  H <- matrix(NA, nsim, 2) # Std error with hessian
  S <- matrix(NA, nsim, 2) # std error with Sandwich estimator

  for(i in 1:nsim){

  if(i<10 || i%%10==0) print(i)

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
  # g <- grad(func = GaussLogLik, x = out$par, data = x)

  # standard erros with sandwhich method
  logfi <- function(theta, idx){
    LAM <- theta[1]
    PHI <- theta[2]
    HC <- HermCoef(LAM)
    ar.acf <- ARMAacf(ar = PHI, lag.max = n+1)[1:(n+1)]
    g <- CountACVF(h = 0:(n-1), myacf = ar.acf, g = HC)
    # g <- ARMAacf(ar = PHI, lag.max = n)[2:(n+1)]
    DLout <- DLalg(x, g, LAM)
    ei <- DLout$e[idx]
    vi <- DLout$v[idx]
    return(-(log(vi) + ei^2/vi)/2)
  }
  gi <- matrix(-99, n, 2)
  for(k in 1:99) gi[k, ] <- grad(func = logfi, x = out$par, idx = k)
  gi <- gi[-100, ]

  A <- h
  B <- t(gi) %*% gi
  V.sand  <- solve(-A) %*% B %*% solve(-A)
  SE.sand <- sqrt(diag(V.sand))

  # standard errors usual way using hessian
  SE.hess <- sqrt(diag(solve(h)))

  P[i, ] <- out$par
  H[i, ] <- SE.hess
  S[i, ] <- SE.sand

  }

  return(list(P = P,
              H = H,
              S = S))

}



