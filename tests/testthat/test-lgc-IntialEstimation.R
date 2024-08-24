#==========================================================================================#
# Purpose: Test the bias of initial estimation. I ll do 20 realizations with sample size =100
#           from each distribution and check that relative bias is not too bad.
#
# Author: Stefanos Kechagias
#==========================================================================================#

test_that("Initial Estimation for Poisson-AR(2)", {

# regressor variable with intercept
Regressor      = NULL
n              = 100
ARMAModel      = c(2,0)
ARParm         = c(0.5, 0.2)
MAParm         = NULL
nsim           = 20
CountDist      = "Poisson"
MargParm       = 3
TrueParam      = c(MargParm,ARParm, MAParm)
theta = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
for (i in 1:nsim){
  # simulate data
  set.seed(i)
  DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

  # call the wrapper function with less arguments
  mod = ModelScheme(DependentVar   = DependentVar,
            Regressor      = Regressor,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel)

  theta[i,] = InitialEstimates(mod)
}

BIAS = colMeans(theta) - TrueParam
RELATIVEBIAS = BIAS/TrueParam
# check that relative bias is small for all
expect_identical(prod(abs(RELATIVEBIAS)<0.1), 1)

# BIAS
# 0.02050000 -0.03782252 -0.01822419
# RELATIVEBIAS
# 0.006833333 -0.075645049 -0.091120948
})

test_that("Initial Estimation for Negative Binomial-MA(2)", {

  # regressor variable with intercept
  Regressor      = NULL
  n              = 100
  ARMAModel      = c(0,2)
  ARParm         = NULL
  MAParm         = c(0.5, 0.2)
  nsim           = 20
  CountDist      = "Negative Binomial"
  MargParm       = c(5,0.3)
  TrueParam      = c(MargParm,ARParm, MAParm)
  theta          = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
  for (i in 1:nsim){
    # simulate data
    set.seed(i)
    DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

    # call the wrapper function with less arguments
    mod = ModelScheme(DependentVar   = DependentVar,
                      Regressor      = Regressor,
                      CountDist      = CountDist,
                      ARMAModel      = ARMAModel)

    theta[i,] = InitialEstimates(mod)
  }

  BIAS = colMeans(theta) - TrueParam
  RELATIVEBIAS = BIAS/TrueParam
  expect_equal(RELATIVEBIAS[1], 0.29426939, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[2], -0.05032221, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[3], -0.14160890, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[4], 0.02372624, tolerance = 10^(-5))

  # BIAS
  # 0.02050000 -0.03782252 -0.01822419
  # RELATIVEBIAS
  # 0.006833333 -0.075645049 -0.091120948
})

test_that("Initial Estimation for Mixed Poisson-ARMA(1,1)", {

  # regressor variable with intercept
  Regressor      = NULL
  n              = 100
  ARMAModel      = c(1,1)
  ARParm         = 0.7
  MAParm         = 0.2
  nsim           = 20
  CountDist      = "Mixed Poisson"
  MargParm       = c(2,10,0.5)
  theta          = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
  TrueParam      = c(MargParm,ARParm, MAParm)
  for (i in 1:nsim){
    # simulate data
    set.seed(i)
    DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

    # call the wrapper function with less arguments
    mod = ModelScheme(DependentVar   = DependentVar,
                      Regressor      = Regressor,
                      CountDist      = CountDist,
                      ARMAModel      = ARMAModel)

    theta[i,] = InitialEstimates(mod)
  }

  BIAS = colMeans(theta) - TrueParam
  RELATIVEBIAS = BIAS/TrueParam
  expect_equal(RELATIVEBIAS[1], -0.04480000, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[2], -0.00850000, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[3], -0.05260000, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[4], -0.08644158, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[5], -0.12854353, tolerance = 10^(-5))

  # BIAS
  # 0.02050000 -0.03782252 -0.01822419
  # RELATIVEBIAS
  # 0.006833333 -0.075645049 -0.091120948
})

test_that("Initial Estimation for Generalized Poisson-AR(2)", {

  # regressor variable with intercept
  Regressor      = NULL
  n              = 100
  ARMAModel      = c(2,0)
  ARParm         = c(0.5, 0.2)
  MAParm         = NULL
  nsim           = 20
  CountDist      = "Generalized Poisson"
  MargParm       = c(0.3,5)
  theta          = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
  TrueParam      = c(MargParm,ARParm, MAParm)

  for (i in 1:nsim){
    # simulate data
    set.seed(i)
    DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

    # call the wrapper function with less arguments
    mod = ModelScheme(DependentVar   = DependentVar,
                      Regressor      = Regressor,
                      CountDist      = CountDist,
                      ARMAModel      = ARMAModel)

    theta[i,] = InitialEstimates(mod)
  }

  BIAS = colMeans(theta) - TrueParam
  RELATIVEBIAS = BIAS/TrueParam

  expect_equal(RELATIVEBIAS[1], -0.02076303, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[2], 0.01380000, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[3], -0.18200882, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[4], -0.12706686, tolerance = 10^(-5))


})

test_that("Initial Estimation for Generalized Poisson 2-AR(2)", {

  # regressor variable with intercept
  Regressor      = NULL
  n              = 100
  ARMAModel      = c(2,0)
  ARParm         = c(0.5, 0.2)
  MAParm         = NULL
  nsim           = 20
  MargParm       = c(0.3,5)
  theta          = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
  TrueParam      = c(MargParm,ARParm, MAParm)

  for (i in 1:nsim){
    # simulate data
    set.seed(i)
    # at this point in time I haven added the Gen Pois 2 in the sim_lgc
    CountDist      = "Generalized Poisson"
    DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)
    CountDist      = "Generalized Poisson 2"

    # call the wrapper function with less arguments
    mod = ModelScheme(DependentVar   = DependentVar,
                      Regressor      = Regressor,
                      CountDist      = CountDist,
                      ARMAModel      = ARMAModel)

    theta[i,] = InitialEstimates(mod)
  }

  BIAS = colMeans(theta) - TrueParam
  RELATIVEBIAS = BIAS/TrueParam

  expect_equal(RELATIVEBIAS[1], -0.02076303, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[2], 0.01380000, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[3], -0.18200882, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[4], -0.12706686, tolerance = 10^(-5))


})

# test_that("Initial Estimation for Generalized Poisson-AR(2) with Regressor", {
#
#   # regressor variable with intercept
#   n              = 100
#   b0             = 0.5
#   b1             = 2
#   alpha          = 0.5
#   MargParm       = c(alpha, b0,b1)
#   e              = rbinom(n,1,0.2)
#   Regressor      = e
#   ARMAModel      = c(2,0)
#   ARParm         = c(0.5, 0.2)
#   MAParm         = NULL
#   nsim           = 20
#   CountDist      = "Generalized Poisson"
#   theta          = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
#   TrueParam      = c(MargParm,ARParm, MAParm)
#
#   for (i in 1:nsim){
#     # simulate data
#     set.seed(i)
#     DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)
#
#     # call the wrapper function with less arguments
#     mod = ModelScheme(DependentVar   = DependentVar,
#                       Regressor      = Regressor,
#                       CountDist      = CountDist,
#                       ARMAModel      = ARMAModel)
#
#     theta[i,] = InitialEstimates(mod)
#   }
#
#   BIAS = colMeans(theta) - TrueParam
#   RELATIVEBIAS = BIAS/TrueParam
#
#   expect_equal(RELATIVEBIAS[1], -0.02076303, tolerance = 10^(-5))
#   expect_equal(RELATIVEBIAS[2], 0.01380000, tolerance = 10^(-5))
#   expect_equal(RELATIVEBIAS[3], -0.18200882, tolerance = 10^(-5))
#   expect_equal(RELATIVEBIAS[4], -0.12706686, tolerance = 10^(-5))
#
#
# })
