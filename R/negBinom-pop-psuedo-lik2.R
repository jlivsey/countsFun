myV_negBinom_2 <- function(size.tru,
                         prob.tru,
                         phi.tru,
                         size.free,
                         prob.free,
                         phi.free){

    # For testing only
    # size.tru <- 3
    # prob.tru <- 1/4
    # phi.tru <- 1/2
    # size.free <- 3
    # prob.free <- 1/4
    # phi.free<- 1/2


    # Calculate gam.x of the true process
    hc <- HermCoefNegBin(size.tru, prob.tru)
    gam.z <- ARMAacf(ar = phi.tru, lag.max = n)
    gam.x <- CountACVF(h = 1:n, myacf = gam.z, g = hc)
    gam.x.tru <- gam.x # define here

    # repeate to get free parameter ACF
    hc <- HermCoefNegBin(size.free, prob.free)
    gam.z <- ARMAacf(ar = phi.free, lag.max = n)
    gam.x <- CountACVF(h = 1:n, myacf = gam.z, g = hc)
    gam.x.free <- gam.x # defined here

    # Yule-Walker vector and matrix
    gam0 <- gam.x.tru[-1]
    gam1 <- gam.x.free[-1]
    G0 <- toeplitz(gam.x.tru[-(n-1)])
    G1 <- toeplitz(gam.x.free[-(n-1)])

    # Variance of each negative binomial process
    nb.var.tru <- size.tru * (1 - prob.tru) / prob.tru^2
    nb.var.free <- size.free * (1 - prob.free) / prob.free^2

    # terms for the output
    term1 <- nb.var.tru
    term2 <- -2 * nb.var.tru * t(gam0) %*% solve(G1) %*% gam1
    term3 <- nb.var.tru * t(gam1) %*% solve(G1) %*% G0 %*% solve(G1) %*% gam1

    # output
    V <- (term1 + term2 + term3) / det(G1)^(1/n) + log(det(G1))/n
    return(V)
}
