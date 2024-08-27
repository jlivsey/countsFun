#============================================================================================#
# PURPOSE: Check that the VGAM distribution functions of GenPois and the ones I have been using
#          are the same.-- I would like to switch into using the VGAM ones after adequate testing.
#
# Author: Stefanos Kechagias
# Date:   August 2024
#============================================================================================#


test_that("VGAM versus our implementation of GenPois", {

# Specify model and methods
n              = 100
CountDist      = "Generalized Poisson"
alpha          = 1
mu             = 3
MargParm       = c(alpha,mu)
ARParm         = 0.75
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))

# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm)

# Run the wrapper
mod = ModelScheme(DependentVar   = DependentVar,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel)

# get Initial Estimate
theta1 = InitialEstimates(mod)

# compute loglik
set.seed(1)
a1 = ParticleFilter_Res_ARMA(theta1, mod)

# VGAM distributions
CountDist      = "Generalized Poisson 2"

# Run the wrapper
mod = ModelScheme(DependentVar   = DependentVar,
                  CountDist      = CountDist,
                  ARMAModel      = ARMAModel)

#
theta2 = InitialEstimates(mod)

set.seed(1)
a2 = ParticleFilter_Res_ARMA(theta2, mod)

set.seed(1)
# Run the wrapper
mylgc = lgc(DependentVar   = DependentVar,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel)


# check the likelihoods
expect_equal(a1,a2,tolerance = 10^(-10))

# check the likelihoods through the lgc wrapper
expect_equal(a1,as.numeric(mylgc$FitStatistics[1]),tolerance = 10^(-10))

# check the initial estimates
expect_equal(theta1[1],theta2[1],tolerance = 10^(-10))
expect_equal(theta1[2],theta2[2],tolerance = 10^(-10))
expect_equal(theta1[3],theta2[3],tolerance = 10^(-10))



})

test_that("VGAM handles NA bettter", {

  # select a setting that our implementation of GenPois can not compute
  a = 0.5
  p = 0.99
  mu = 15
  x1 = qgenpois2(p,mu,a)
  x2 = qGpois(p,a,mu)

  # qgenpois2 works
  expect_identical(x1,162)

  # qGpois produces NA
  expect_identical(is.na(x2),TRUE)

})
