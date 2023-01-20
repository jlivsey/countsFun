# Purpose: Calculate coverage of CI's calculated with Sandwich estimators vs
#           not using them.
#           Model is Poisson - AR(1)
library(countsFun)
library(numDeriv)
library(sandwich)

sim_sandy_fixedPhi <- function(lam, phi, n = 100, nsim = 200){

  P <- matrix(NA, nsim, 1) # parameter estimates
  H <- matrix(NA, nsim, 1) # Std error with hessian
  S <- matrix(NA, nsim, 1) # std error with Sandwich estimator

  # Main loop
  for(i in 1:nsim){

    # Print if i less than 10 or multiple of 10 to track progress
    if(i<10 || i%%10==0) print(i)

  # simulate data
  if(phi == 0){
    x <- rpois(n, lam)
  }else{
    x <- sim_pois_ar(n = n, phi = phi, lam = lam)
  }

  # set initial values for optim
  inital.values <- c(jitter(lam))

  # run optim
  out <- optim(par     = inital.values,
               fn      = GaussLogLik_fixedPhi,
               phi_fixed = phi,
               data    = x,
               hessian = TRUE)

  # Calulate numerical Hessian
  h <- hessian(func = GaussLogLik_fixedPhi,
               x = out$par,
               data = x,
               phi_fixed = phi)

  # Calculate numerical gradient
  g <- grad(func = GaussLogLik_fixedPhi,
            x = out$par,
            data = x,
            phi_fixed = phi)

  # standard erros with sandwhich method
  logfi <- function(theta, idx){
    LAM <- theta[1]
    PHI <- phi
    HC <- HermCoef(LAM)
    ar.acf <- ARMAacf(ar = PHI, lag.max = n+1)[1:(n+1)]
    g <- CountACVF(h = 0:(n-1), myacf = ar.acf, g = HC)
    # g <- ARMAacf(ar = PHI, lag.max = n)[2:(n+1)]
    DLout <- DLalg(x, g, lam = LAM)
    ei <- DLout$e[idx]
    vi <- DLout$v[idx]
    return(-(log(vi) + ei^2/vi)/2)
  }
  gi <- matrix(-99, n, 1)
  for(k in 1:99) gi[k, ] <- grad(func = logfi, x = out$par, idx = k)

  # AD-HOC REMOVE
  gi <- gi[-100]

  A <- h
  B <- t(gi) %*% gi
  V.sand  <- solve(-A) %*% B %*% solve(-A)
  SE.sand <- sqrt(diag(V.sand))

  # standard errors usual way using hessian
  SE.hess <- sqrt(diag(solve(h)))

  # Outputs
  P[i, ] <- out$par
  H[i, ] <- SE.hess
  S[i, ] <- SE.sand

  }

  return(list(P = P,
              H = H,
              S = S))

}



