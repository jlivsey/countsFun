#======================================================================================================#
#Purpose:   Function to evalutate Sandwich Std Errors for Count Models for NegBinom-AR(1) model
#
# Author:   Stefanos Kechagias and James Livsey
# Team:     Vladas Pipiras, James Livsey, Stefanos Kechagias, Robert Lund, Yisu Jia
# Date:     July 2020
#=====================================================================================================#

library(countsFun)
library(numDeriv)
library(sandwich)


sand <- function(ARMAorder, thetaEst, data){

  # vvvv Temporary inputs to test function vvvv
  ARMAorder <- c(1, 0)
  theta <- c(5, 1/4, 1/2)

  r <- theta[1]
  p <- theta[2]
  phi <- theta[3]

  data <- sim_negbin(n = 200, ARMAmodel = list(phi, 0), r = r, p = p)

  n <- length(data)

  theta.est <- c(jitter(theta[1]),
                 jitter(theta[2]),
                 jitter(theta[3]))

  # ^^^^ END of temporary inputs ^^^^

  # Calulate numerical Hessian
  h <- hessian(func = GaussLogLikNB, x = theta.est,
               ARMAorder = ARMAorder, # additional arg for GaussLogLikNB
               MaxCdf = 1000,         # additional arg for GaussLogLikNB
               data = data,           # additional arg for GaussLogLikNB
               nHC = 20)              # additional arg for GaussLogLikNB

  # standard erros with sandwhich method following notation of Freedman (2006)
  logfi <- function(theta, idx){
    p <- theta[1]
    r <- theta[2]
    phi <- theta[3]
    HC <- HermCoefNegBin(r = r, p = p, N = 1000, nHC = 20)
    ar.acf <- ARMAacf(ar = phi, lag.max = n+1)[1:(n+1)]
    g <- CountACVF(h = 0:(n-1), myacf = ar.acf, g = HC)
    DLout <- DLalg(data, g, p*r / (1-p))
    ei <- DLout$e[idx]
    vi <- DLout$v[idx]
    return(-(log(vi) + ei^2/vi)/2)
  }
  gi <- matrix(-99, n, length(theta))
  for(k in 1:99){
    print(k)
    gi[k, ] <- grad(func = logfi, x = theta.est, idx = k)
  }
  gi <- gi[-n, ]

  # Calculate Cov matrix
  A <- h
  B <- t(gi) %*% gi
  V.sand  <- solve(-A) %*% B %*% solve(-A)
  SE.sand <- sqrt(diag(V.sand))

  # standard errors usual way using hessian
  SE.hess <- sqrt(diag(solve(h)))

}
