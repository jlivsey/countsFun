#====================================================================================================#
# Purpose: Test the performance of initial estimation. Not extensive.
#
# Author: Stefanos Kechagias
#==========================================================================================#

#--------------------------- without Regressor ------------------------------------------#
test_that("Initial Estimation for Poisson-AR(2)", {

# regressor variable with intercept
Regressor      = NULL
n              = 100
ARParm         = c(0.5, 0.2)
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))
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
  ARParm         = NULL
  MAParm         = c(0.5, 0.2)
  ARMAModel      = c(length(ARParm),length(MAParm))
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
  expect_equal(RELATIVEBIAS[3], -0.10107972, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[4], 0.10518289, tolerance = 10^(-5))

  # BIAS
  # 0.02050000 -0.03782252 -0.01822419
  # RELATIVEBIAS
  # 0.006833333 -0.075645049 -0.091120948
})

test_that("Initial Estimation for Mixed Poisson-ARMA(1,1)", {

  # regressor variable with intercept
  Regressor      = NULL
  n              = 100

  ARParm         = 0.7
  MAParm         = 0.2
  ARMAModel      = c(length(ARParm),length(MAParm))

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
  expect_equal(RELATIVEBIAS[4], -0.07637862, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[5], -0.06920126, tolerance = 10^(-5))

  # BIAS
  # 0.02050000 -0.03782252 -0.01822419
  # RELATIVEBIAS
  # 0.006833333 -0.075645049 -0.091120948
})

test_that("Initial Estimation for Generalized Poisson-AR(2)", {

  # regressor variable with intercept
  Regressor      = NULL
  n              = 100
  ARParm         = c(0.5, 0.2)
  MAParm         = NULL
  ARMAModel      = c(length(ARParm),length(MAParm))
  CountDist      = "Generalized Poisson 2"
  alpha          = 0.3
  mu             = 5
  MargParm       = c(alpha,mu)
  TrueParam      = c(MargParm,ARParm, MAParm)
  nsim           = 20
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

  expect_equal(RELATIVEBIAS[1], -0.02076303, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[2], 0.01380000, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[3], -0.08716337, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[4], -0.07661555, tolerance = 10^(-5))


})

test_that("Initial Estimation for Binomial-AR(2)", {
n              = 100
Regressor      = NULL
ARParm         = c(0.8, -0.5)
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))
nsim           = 20
CountDist      = "Binomial"
ntrials        = 5
p              = 0.3
MargParm       = p
TrueParam      = c(MargParm,ARParm, MAParm)
theta = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))

for (i in 1:nsim){
  # simulate data
  set.seed(i)
  DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor,ntrials)

  # call the wrapper function with less arguments
  mod = ModelScheme(DependentVar   = DependentVar,
                    Regressor      = Regressor,
                    CountDist      = CountDist,
                    ARMAModel      = ARMAModel,
                      ntrials      = ntrials)

  theta[i,] = InitialEstimates(mod)
}

BIAS = colMeans(theta) - TrueParam
RELATIVEBIAS = BIAS/TrueParam


expect_equal(RELATIVEBIAS, c(-0.0003333333, -0.1532225884, -0.2103357704))
})

test_that("Initial Estimation for ZIP-ARMA(1,1)", {

  # regressor variable with intercept
  n              = 50
  MAParm         = 0.2
  nsim           = 20
  ARParm         = 0.5
  ARMAModel      = c(length(ARParm),length(MAParm))
  CountDist      = "ZIP"
  lambda         = 3
  p              = 0.2
  MargParm       = c(lambda,p)
  theta          = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
  TrueParam      = c(MargParm,ARParm, MAParm)


  for (i in 1:nsim){
    # simulate data
    set.seed(i)
    DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm)

    # call the wrapper function with less arguments
    mod = ModelScheme(DependentVar   = DependentVar,
                      CountDist      = CountDist,
                      ARMAModel      = ARMAModel)

    theta[i,] = InitialEstimates(mod)
  }

  BIAS = colMeans(theta) - TrueParam
  RELATIVEBIAS = BIAS/TrueParam
  expect_equal(RELATIVEBIAS[1], 0.04668120,tolerance=10^(-5))
  expect_equal(RELATIVEBIAS[2], -0.04539234,tolerance=10^(-5))
  expect_equal(RELATIVEBIAS[3], -0.22390362,tolerance=10^(-5))
  expect_equal(RELATIVEBIAS[4], 0.24277838,tolerance=10^(-5))

})


#--------------------------- with Regressor ------------------------------------------#

test_that("Initial Estimation for Generalized Poisson-AR(1) with Regressor", {
# set parameters
CountDist      = "Generalized Poisson 2"
alpha          = 1.5
b0             = 1
b1             = 4
MargParm       = c(b0,b1,alpha)
ARParm         = 0.6
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))
TrueParam      = c(MargParm,ARParm, MAParm)
n              = 100
#e              = rbinom(n,1,0.1)
#Regressor      = cbind(rep(1,n),e)
set.seed(3)
e              = runif(n)
Regressor      = cbind(rep(1,n),e)
nsim           = 20
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

# check me: this doesnt seem to work as well. what about IYW?
expect_equal(RELATIVEBIAS[1], -0.11144035 , tolerance = 10^(-5))
expect_equal(RELATIVEBIAS[2], 0.08922739 , tolerance = 10^(-5))
expect_equal(RELATIVEBIAS[3], -0.03059165, tolerance = 10^(-5))
expect_equal(RELATIVEBIAS[4], -0.22826401, tolerance = 10^(-5))

})

test_that("Initial Estimation for Poisson-AR(2) with Regressor", {

  # regressor variable with intercept

  n              = 100
  set.seed(3)
  e              = runif(n)
  Regressor      = cbind(rep(1,n),e)
  ARMAModel      = c(2,0)
  ARParm         = c(0.5, 0.2)
  MAParm         = NULL
  nsim           = 20
  CountDist      = "Poisson"
  b0             = 1
  b1             = 4
  MargParm       = c(b0,b1)
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

  expect_equal(RELATIVEBIAS[1], -0.0012539700, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[2], 0.0004266268, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[3],  -0.067607045, tolerance = 10^(-5))
  expect_equal(RELATIVEBIAS[4], -0.0982610050, tolerance = 10^(-5))

})

test_that("Initial Estimation for Binomial-MA(3) with Regressor", {

  # set sample size and model Parameters
  n              = 100
  ARParm         = NULL
  #MAParm         = c(0.8, -0.5, 0.4)
  #MAParm         = c(0.8)
  MAParm          = 0.4
  #z  =arima.sim(model = list( ar = ARParm, ma=MAParm), n = n)

  ARMAModel      = c(length(ARParm),length(MAParm))
  nsim           = 20
  CountDist      = "Binomial"

  # select parameters for linear predictor
  b0    = 0.2
  b1    = 2
  beta  = c(b0,b1)

  # set the constant binomial parameter
  ntrials        = 1

  # collect Marginal and True parameters
  MargParm       = c(b0,b1)
  TrueParam      = c(MargParm,ARParm, MAParm)

  # allocate memory for estimates
  theta = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
  set.seed(3)
  Regressor = cbind(1,rnorm(n,0,1))

    for (i in 1:nsim){
    # simulate data
    set.seed(i)
    # Generate a regressor that will be used as a linear predictor
    DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor,ntrials)

    # call the wrapper function with less arguments
    mod = ModelScheme(DependentVar   = DependentVar,
                      Regressor      = Regressor,
                      CountDist      = CountDist,
                      ARMAModel      = ARMAModel,
                      ntrials      = ntrials)

    theta[i,] = InitialEstimates(mod)
  }

  BIAS = colMeans(theta) - TrueParam
  RELATIVEBIAS = BIAS/TrueParam


  expect_equal(RELATIVEBIAS, c( -0.20985523,  0.0252643, -0.67265663))
})

test_that("Initial Estimation for ZIP-AR(1) with Rerggresor", {

  n              = 50
  nsim           = 20
  ARParm         = 0.5
  MAParm         = NULL
  ARMAModel      = c(length(ARParm),length(MAParm))
  CountDist      = "ZIP"

  # Set parameter for the Poisson link
  b0             = 2
  b1             = 4
  beta           = c(b0,b1)
  p              = 0.4
  MargParm       = c(b0,b1,p)
  theta          = matrix(NA,nrow=nsim,ncol=length(MargParm)+sum(ARMAModel))
  TrueParam      = c(MargParm,ARParm, MAParm)

  # Generate a regressor
  set.seed(3)
  Regressor  = cbind(1,runif(n))


  for (i in 1:nsim){
    # simulate data
    set.seed(i)
    DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

    # call the wrapper function with less arguments
    mod = ModelScheme(DependentVar   = DependentVar,
                      CountDist      = CountDist,
                      ARMAModel      = ARMAModel,
                      Regressor      = Regressor)

    theta[i,] = InitialEstimates(mod)
  }

  BIAS = colMeans(theta) - TrueParam
  RELATIVEBIAS = BIAS/TrueParam
  expect_equal(RELATIVEBIAS[1], 0.016037874,tolerance=10^(-5))
  expect_equal(RELATIVEBIAS[2], -0.009261874,tolerance=10^(-5))
  expect_equal(RELATIVEBIAS[3], -0.052515431,tolerance=10^(-5))
  expect_equal(RELATIVEBIAS[4], -0.208396106,tolerance=10^(-5))

})
