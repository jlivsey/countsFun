# Purpose: Calculate coverage of CI's calculated with Sandwich estimators vs
#           not using them.
#           Model is Poisson - AR(1)
library(countsFun)
library(numDeriv)
library(sandwich)

sim_sandy <- function(lam, phi, n = 100, nsim = 200){

  # simulate data
  x <- sim_pois_ar(n = n, phi = phi, lam = lam)

  x <- rpois(n, lam)

  # set initial values for optim
  inital.values <- c(jitter(lam), jitter(phi))
  inital.values <- c(jitter(lam))


  # run optim
  out <- optim(par     = inital.values,
               fn      = GaussLogLik_fixedPhi,
               phi_fixed = 0,
               data    = x,
               hessian = TRUE)

  # Calulate numerical Hessian
  h <- hessian(func = GaussLogLik_fixedPhi,
               x = out$par,
               data = x,
               phi_fixed = 0)

  # Calculate numerical gradient
  g <- grad(func = GaussLogLik_fixedPhi,
            x = out$par,
            data = x,
            phi_fixed = 0 )

  # standard errors usual way using hessian
  SE.hess <- sqrt(diag(solve(h)))

  # standard erros with sandwhich method
  logfi <- function(theta, idx){
    LAM <- theta[1]
    PHI <- 0
    HC <- HermCoef(LAM)
    ar.acf <- ARMAacf(ar = PHI, lag.max = n+1)[1:(n+1)]
    g <- CountACVF(h = 0:(n-1), myacf = ar.acf, g = HC)
    # g <- ARMAacf(ar = PHI, lag.max = n)[2:(n+1)]
    DLout <- DLalg(x, g)
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
  # Does true parameter fall in CI


}



