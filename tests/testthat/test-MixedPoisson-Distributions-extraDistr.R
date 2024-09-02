#===========================================================================================#
# Purpose: Compare our mixed Poisson distribution functions with the ones form extraDistr
#          for inverse cdf extraDistr doesn't have something
#          spatstat package also has some implementation of mix Poisson but it is a different
#          parametrization
#
# Author: Stefanos Kechagias
# Date:   September 2024
#===========================================================================================#

test_that("Implementations of Mixed Poisson Distributions", {
# Constant parameters
lambda1 = 2
lambda2 = 3
p = 0.4
x = 1

a1 = dmixpois(x, c(lambda1,lambda2), c(p,1-p))
a2 = dmixpois1(x,lambda1,lambda2,p)
expect_equal(a1, a2)

# vector input and constant parameters
p = 0.4
x = c(1,1)
lambda1 = 2
lambda2 = 3

a1 = dmixpois(x, c(lambda1,lambda2), c(p,1-p))
a2 = dmixpois1(x,lambda1,lambda2,p)

expect_equal(a1, a2)
# vector input and vector constant parameters
p = 0.4
x = c(1,1)
lambda1 = c(2,1)
lambda2 = c(3,4)

a1 = dmixpois(x, matrix(c(lambda1,lambda2),ncol=2), c(p,1-p))
a2 = dmixpois1(x,lambda1,lambda2,p)
expect_equal(a1, a2)

#########################################################


# set linear predictor parameters
b0_1 = 2
b1_1 = 0.5

b0_2 = 3
b1_2 =  1

# set mixing probability
prob = 0.4

# gather Marginal Parameters
MargParm = c(b0_1, b1_1, b0_2, b1_2,prob)

# set other model details
n = 50
CountDist = "Mixed Poisson"
ARParm = NULL
MAParm = NULL
#ARMAModel      = c(length(ARParm),length(MAParm))

# generate regressor
x = runif(n)
Regressor = cbind(1,x)

# coimpute the Poisson means
lambda1 = exp(Regressor%*%c(b0_1,b1_1))
lambda2 = exp(Regressor%*%c(b0_2,b1_2))

# collect  Marginal Parameters
DynamMargParm = matrix(c(lambda1,lambda2),ncol=2)
ConstMargParm = prob

# generate some counts
set.seed(2)
x = round(10*pnorm(rnorm(n)))+1

# compute PMF using both functions
a1 = dmixpois(x, DynamMargParm, c(prob,1-prob))
a2 = dmixpois1(x,lambda1,lambda2,prob)

expect_equal(a1, a2)
###########################################################

# compare cdfs with linear predictors
a1 = pmixpois(x, DynamMargParm, c(prob,1-prob))
a2 = pmixpois1(x, DynamMargParm[,1], DynamMargParm[,2], prob)
expect_equal(a1, a2)

#########################################################################

# set linear predictor parameters
n = 10
b0_1 = 2
b1_1 = 0.5

b0_2 = 3
b1_2 =  1

prob = 0.4
# generate regressor
set.seed(1)
x = round(10*pnorm(rnorm(n)))+1
Regressor = cbind(1,x)

# coimpute the Poisson means
lambda1 = exp(Regressor%*%c(b0_1,b1_1))
lambda2 = exp(Regressor%*%c(b0_2,b1_2))

p = pmixpois(x,lambda1,lambda2,prob)

# Both of the functions below are our implementations
a1 = qmixpois(p,lambda1,lambda2,prob)

a2 = qmixpois1(p,lambda1,lambda2,prob)
expect_equal(a1, a2)

})

























