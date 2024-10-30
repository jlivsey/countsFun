test_that("LGC wrapper works for evaluation MixPois AR(1)", {

  # Specify model and methods
  n              = 200
  Regressor      = NULL
  CountDist      = "Mixed Poisson"
  MargParm       = c(2, 5, 0.7)
  ARParm         = 0.75
  MAParm         = NULL
  ARMAModel      = c(length(ARParm),length(MAParm))
  ParticleNumber = 5
  epsilon        = 0.5
  EstMethod      = "PFR"
  TrueParam      = NULL
  initialParam   = NULL
  Task           = 'Evaluation'
  SampleSize     = NULL
  nsim           = NULL
  no_cores       = NULL
  OptMethod      = "bobyqa"
  OutputType     = "list"
  ParamScheme    = NULL
  maxdiff        = 10^(-8)

  # simulate data
  set.seed(2)
  DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm)

  # save the data in a data frame
  df = data.frame(DependentVar)

  # specify the regression model
  formula = DependentVar~0

  # Run the wrapper
  mylgc = lgc(formula        = formula,
              data           = df,
              EstMethod      = EstMethod,
              CountDist      = CountDist,
              ARMAModel      = ARMAModel,
              ParticleNumber = ParticleNumber,
              epsilon        = epsilon,
              initialParam   = initialParam,
              TrueParam      = TrueParam,
              Task           = Task,
              SampleSize     = SampleSize,
              nsim           = nsim,
              no_cores       = no_cores,
              OptMethod      = OptMethod,
              OutputType     = OutputType,
              ParamScheme    = ParamScheme,
              maxdiff        = maxdiff)
  est = c(1.602, 4.163,0.507,0.67716673)

  expect_equal(mylgc$FitStatistics[[1]], 369.35048)
  expect_equal(as.numeric(mylgc$ParamEstimates), est)


})
