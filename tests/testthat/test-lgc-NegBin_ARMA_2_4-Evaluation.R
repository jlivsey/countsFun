test_that("LGC wrapper works for evaluation NegBin ARMA(2,4)", {

  # Specify model and methods
  n              = 200
  Regressor      = NULL
  CountDist      = "Negative Binomial"
  MargParm       = c(3,0.8)
  ARParm         = c(0.4, 0.4)
  MAParm         = c(0.3, 0.5, -0.2, 0.1)
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

  est = c(3.02527116, 0.80644794 , 0.17345437 , 0.48201984 , 0.47124456, 0.51209982 , -0.10668423, 0.008958203)

  expect_equal(mylgc$FitStatistics[[1]], 560.81726 )
  expect_equal(as.numeric(mylgc$ParamEstimates), est)
})
