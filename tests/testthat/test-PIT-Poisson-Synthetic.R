#======================================================================================================#
#Purpose:   Evaluate the Dominick data while providing fewer arguments.
#
# Author:   Stefanos Kechagias
# Date:     April 2026
#=====================================================================================================#

test_that("PIT histogram is approximately uniform", {

  # regressor variable with intercept
  CountDist      = "Poisson"
  Regressor      = NULL
  Intercept      = NULL
  n              = 200
  ARMAModel      = c(2,0)
  MargParm       = 3
  ARParm         = c(0.5, 0.2)
  MAParm         = NULL
  OptMethod      = "L-BFGS-B"
  Task           = "Optimization"

  # specify the regression formula (no regressors here)
  RegModel       = DependentVar ~ 1

  # simulate data
  set.seed(2)
  DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, RegModel)
  df = data.frame(DependentVar)


  # fit the model
  fit = lgc(RegModel = RegModel,
            df = df,
            CountDist = CountDist,
            ARMAModel = ARMAModel,
            Task = Task)

  # Retrieve the model Specification
  ModelSpec = model(fit)

  # set a parameter point
  theta = c( MargParm, ARParm, MAParm)

  # compute predictive distribution
  predDist = pred_dist(theta, ModelSpec)

  H = 10
  # compute PIT values
  PIT = pit(H, predDist)

  # approx chi sq (I dont ahve independence)
  chisq_stat <- sum((n * PIT - n / H)^2 / (n / H))
  expect_lt(chisq_stat, qchisq(0.95, df = H - 1))
})
