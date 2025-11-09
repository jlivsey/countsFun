test_that("Likelihood Poisson-AR(1)", {


  # Specify model and methods
  n              = 200
  CountDist      = "Poisson"
  MargParm       = 3
  ARParm         = c(0.8, -0.25)
  MAParm         = c(0.2, 0.5,0.2, 0.1)
  #ARParm         = NULL
  #MAParm         = NULL
  ARMAModel      = c(length(ARParm),length(MAParm))
  ParticleNumber = 10
  epsilon        = 0.5
  EstMethod      = "PFR"
  TrueParam      = c(MargParm,ARParm,MAParm)
  initialParam   = TrueParam
  Task           = 'Optimization'
  SampleSize     = NULL
  nsim           = NULL
  no_cores       = NULL
  OptMethod      = "bobyqa"
  OutputType     = "data.frame"
  ParamScheme    = NULL
  maxdiff        = 10^(-6)


  # specify the regression formula (no regressors here)
  RegModel       = DependentVar ~ 1

  # simulate data
  set.seed(2)
  DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, RegModel)
  df = data.frame(DependentVar)

    # populate a list with the model characteristics
  mod = ModelScheme(RegModel = RegModel,
                    df = df,
                    EstMethod = EstMethod,
                    CountDist = CountDist,
                    ARMAModel = ARMAModel,
                    OptMethod = OptMethod,
                    ParticleNumber = ParticleNumber,
                    epsilon = epsilon,
                    initialParam = initialParam,
                    TrueParam = TrueParam,
                    Task = Task,
                    SampleSize = SampleSize,
                    OutputType = OutputType,
                    ParamScheme = ParamScheme,
                    maxdiff = maxdiff)


  # evaluate the likleihood
  a1 = ParticleFilter(initialParam,mod)

  # for the set.seed(2) I will get the following likelihood: 302.57734
  expect_equal(a1,289.4542947,tolerance=10^(-5))


})
