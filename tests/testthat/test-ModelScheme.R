test_that("parsing works ok", {

  # Specify model and methods
  n              = 200
  # Regressor      = cbind(rep(1,n),rbinom(n,1,0.25))
  Regressor      = NULL
  CountDist      = "Poisson"
  MargParm       = 3
  ARParm         = c(0.8, -0.25)
  MAParm         = c(0.2, 0.5,0.2, 0.1)
  #ARParm         = NULL
  #MAParm         = NULL
  ARMAModel      = c(length(ARParm),length(MAParm))
  ParticleNumber = 1
  epsilon        = 0.5
  EstMethod      = "PFR"
  TrueParam      = c(MargParm,ARParm,MAParm)
  initialParam   = TrueParam
  Task           = 'Optimization'
  SampleSize     = NULL
  OptMethod      = "bobyqa"
  OutputType     = "data.frame"
  ParamScheme    = NULL
  maxdiff        = 10^(-6)
  # simulate data
  set.seed(2)
  DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

  mod = ModelScheme(DependentVar, Regressor, EstMethod, ARMAModel, CountDist,ParticleNumber, epsilon,
                    initialParam, TrueParam, Task,SampleSize, OptMethod, OutputType, ParamScheme, maxdiff)

  # the following will check only the variables that enter the ModelScheme function as inputs '
  # and come out directly as outputs. I am not checking here if other variables that are computed
  # inside the ModelScheme function are correct.
  expect_equal(mod$Regressor, Regressor)
  expect_equal(mod$CountDist, CountDist)
  expect_equal(mod$maxdiff, maxdiff)
  expect_equal(mod$ParticleNumber, ParticleNumber)
  expect_equal(mod$OptMethod, OptMethod)
  expect_equal(mod$ParamScheme, ParamScheme)
  expect_equal(mod$TrueParam, TrueParam)
  expect_equal(mod$epsilon, epsilon)
  expect_equal(mod$SampleSize, SampleSize)
  expect_equal(mod$Task, Task)
  expect_equal(as.numeric(mod$initialParam), initialParam)
  expect_equal(mod$EstMethod, EstMethod)
  expect_equal(mod$ARMAModel, ARMAModel)
  expect_equal(mod$OutputType, OutputType)
})


