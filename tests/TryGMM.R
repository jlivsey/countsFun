# example on GMM for normal data
library(gmm)

# moment conditions
g1 <- function(theta,x) {
     m1 <- (theta[1]-x)
     m2 <- (theta[2]^2 - (x - theta[1])^2)
     m3 <- x^3-theta[1]*(theta[1]^2+3*theta[2]^2)
     f <- cbind(m1,m2,m3)
     return(f)
     }
# derivative
Dg <- function(theta, x) {
  G <- matrix(c(1, 2 * (-theta[1] + mean(x)), -3 * theta[1]^2 - 3 * theta[2]^2,
                  + 0, 2 * theta[2], -6 * theta[1] * theta[2]), nrow = 3, ncol = 2)
  return(G)
}


set.seed(123)
n <- 200
x1 <- rnorm(n, mean = 4, sd = 2)
print(res <- gmm(g1, x1, c(mu = 0, sig = 0)))
#=====================================================================================#

library(gmm)
g2 = function(theta,x){
  l1 = theta[1]
  l2 = theta[2]
  p = theta[3]

  m1 = p*l1 + (1-p)*l2 - x
  m2 = p*l1^2 + (1-p)*l2^2 + x - x^2
  m3 = p*l1^3 + (1-p)*l2^3 + 3*x^2 - 2*x - x^3

  f = cbind(m1,m2,m3)
  return(f)
  }

set.seed(1)
n <- 200
x2 <- rmixpois(n,p=0.25,lam1=2, lam2=5)
p0 = c(3,3,0.2)
print(res <- gmm(g2, x2, p0)
print(res <- gel(g2, x2, p0, type = "ETEL")







