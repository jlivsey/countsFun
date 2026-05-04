test_that("Basic_S3_methods_for_Evaluation_Task", {

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

  # specify the regression formula (no regressors here)
  RegModel       = DependentVar ~ 0

  # simulate data
  set.seed(2)
  DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, RegModel)
  df = data.frame(DependentVar)

  # Run the wrapper
  mylgc = lgc(RegModel       = RegModel,
              df           = df,
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

  expect_equal(as.numeric( AIC(mylgc)), 746.710518)
  expect_equal(as.numeric( BIC(mylgc)), 759.9037875)
  expect_equal(as.numeric( logLik(mylgc)), 369.355259)
  expect_equal(model(mylgc), mylgc$mod)
  expect_equal(colnames(coefficients(mylgc))[1], "lambda_1")
  expect_equal(colnames(coefficients(mylgc))[2], "lambda_2")
  expect_equal(colnames(coefficients(mylgc))[3], "p")
  expect_equal(colnames(coefficients(mylgc))[4], "AR_1")
  expect_equal(as.numeric(coefficients(mylgc))[1], 1.599451387)
  expect_equal(as.numeric(coefficients(mylgc))[2], 4.159660347)
  expect_equal(as.numeric(coefficients(mylgc))[3], 0.5056847959)
  expect_equal(as.numeric(coefficients(mylgc))[4], 0.6771269056)


})
