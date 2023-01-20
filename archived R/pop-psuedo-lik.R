#
# myV <- function(lam.tru, phi.tru, lam.free, phi.free){
#
#   # For testing purposes
#   # lam.tru  <- 3
#   # phi.tru  <- 1/2
#   # lam.free <- 3
#   # phi.free <- 1/2
#
#   # Calculate gam.x of the true process
#   hc <- HermCoef(lam.tru)
#   gam.z <- ARMAacf(ar = phi.tru, lag.max = n)
#   gam.x <- CountACVF(h = 1:n, myacf = gam.z, g = hc)
#   gam.x[1] <- 1
#   gam.x.tru <- gam.x # Defined here
#
#   # Calculate gam.x of the free process
#   hc <- HermCoef(lam.free)
#   gam.z <- ARMAacf(ar = phi.free, lag.max = n)
#   gam.x <- CountACVF(h = 1:n, myacf = gam.z, g = hc)
#   gam.x[1] <- 1
#   gam.x.free <- gam.x # defined here
#
#   # YW vector and matrix of autocovariances for true and free process
#   gam0 <- gam.x.tru[-1]
#   gam1 <- gam.x.free[-1]
#   G0 <- toeplitz(gam.x.tru[-(n-1)])
#   G1 <- toeplitz(gam.x.free[-(n-1)])
#
#   # Terms in error expression
#   term1 <- lam.tru
#   term2 <- -2 * lam.tru * t(gam0) %*% solve(G1) %*% gam1
#   term3 <- lam.tru * t(gam1) %*% solve(G1) %*% G0 %*% solve(G1) %*% gam1
#
#   # output
#   V <- term1 + term2 + term3
#   return(V)
# }
#
#
#
#
