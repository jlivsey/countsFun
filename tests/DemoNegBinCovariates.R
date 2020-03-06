# Compare Runtimes of the NegBinomial likelihood evaluation with and without covariates
# I need to vectorize the covariate version. It is too slow.

# libraries
library(countsFun)
library(tictoc)

source('C:/Users/stefa/Desktop/countsFun/R/negBinom-functions.R')

# sample size
n = 10

# generate a regressor dummy variable
X  = cbind(rep(1,n),rbinom(n, 1 ,0.4))

# specify parameters
r = 4
b0 = 1
b1 = 0
b = rbind(b0,b1)
m = exp(X%*%b)
phi = 0.4

# parameters to be passed in likelihood
theta = c(b,r,phi)

# classic parametrization
p = m[1]/(m[1]+r)
theta1 = c(r,p,phi)


# generate negative binomial AR series
y  =sim_negbin_ar_2(n, phi, r, m)

# test likelihood
tic()
GaussLogLikNB_2(theta, y, X)
toc()


tic()
GaussLogLikNB(theta1, y)
toc()
