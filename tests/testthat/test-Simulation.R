#------------------------------------------------------------------------------------------#
# PURPOSE:        Testing lgc for Simulation Task
#
# OUTPUT
#   df            parameter estimates, true values, standard errors, likelihood value and
#                 information on optimization results
#
# DATE:           June 2020
#------------------------------------------------------------------------------------------#

test_that("Simulation for Poisson-AR(1)", {


  #-------------------------- Specify model and other parameters --------------------------#
  SampleSize     = 50
  CountDist      = "Poisson"
  MargParm       = 3
  ARParm         = 0.5
  MAParm         = 0.2
  ARMAModel      = c(length(ARParm),length(MAParm))
  ParticleNumber = 5
  epsilon        = 0.5
  EstMethod      = "PFR"
  initialParam   = NULL
  TrueParam      = c(MargParm,ARParm,MAParm)
  Task           = 'Simulation'
  nsim           = 2
  no_cores       = 1
  OptMethod      = "L-BFGS-B"
  data           = NULL
  RegModel       = DependentVar~0
  set.seed(203)
  #----------------------------------------------------------------------------------------#


  #--------------------------------   Run the main wrapper   ------------------------------#
  df = lgc(RegModel = RegModel,
                 df = df,
          EstMethod = EstMethod,
          CountDist = CountDist,
          ARMAModel = ARMAModel,
     ParticleNumber = ParticleNumber,
            epsilon = epsilon,
       initialParam = initialParam,
          TrueParam = TrueParam,
               Task = Task,
         SampleSize = SampleSize,
               nsim = nsim,
           no_cores = no_cores,
          OptMethod = OptMethod)
  #----------------------------------------------------------------------------------------#


  #------------------------------------   TEST VALUES -------------------------------------#
  ## ---- Row 1 ----
  expect_equal(df[1, "True_lambda"],       3)
  expect_equal(df[1, "lambda"],   3.35522,  tolerance = 1e-5)
  expect_equal(df[1, "se(lambda)"],    0.458925, tolerance = 1e-5)

  expect_equal(df[1, "True_AR_1"],         0.5)
  expect_equal(df[1, "AR_1"],     0.159887, tolerance = 1e-5)
  expect_equal(df[1, "se(AR_1)"],      0.279380, tolerance = 1e-5)

  expect_equal(df[1, "True_MA_1"],         0.2)
  expect_equal(df[1, "MA_1"],     0.5513724, tolerance = 1e-5)
  expect_equal(df[1, "se(MA_1)"],      0.237455,  tolerance = 1e-5)

  expect_equal(df[1, "SampleSize"],            50)
  expect_equal(df[1, "loglik"],       91.3924,  tolerance = 1e-5)
  expect_equal(df[1, "ConvergeStatus"], 0)
  expect_equal(df[1, "kkt1"],         1)
  expect_equal(df[1, "kkt2"],         1)


  ## ---- Row 2 ----
  expect_equal(df[2, "True_lambda"],       3)
  expect_equal(df[2, "lambda"],   3.05718,  tolerance = 1e-5)
  expect_equal(df[2, "se(lambda)"],    0.447802, tolerance = 1e-5)

  expect_equal(df[2, "True_AR_1"],         0.5)
  expect_equal(df[2, "AR_1"],     0.440082, tolerance = 1e-5)
  expect_equal(df[2, "se(AR_1)"],      0.241515, tolerance = 1e-5)

  expect_equal(df[2, "True_MA_1"],         0.2)
  expect_equal(df[2, "MA_1"],     0.0707647, tolerance = 1e-5)
  expect_equal(df[2, "se(MA_1)"],      0.255301,  tolerance = 1e-5)

  expect_equal(df[2, "SampleSize"],            50)
  expect_equal(df[2, "loglik"],       91.0961,  tolerance = 1e-5)
  expect_equal(df[2, "ConvergeStatus"], 0)
  expect_equal(df[2, "kkt1"],         1)
  expect_equal(df[2, "kkt2"],         1)
  #----------------------------------------------------------------------------------------#


  })
