#============================================================================================#
# PURPOSE: Check that the new sim_lgc function yields the same results as sim_lgc_old
#
# Author: Stefanos Kechagias
# Date:   August 2024
#============================================================================================#

test_that("New synthesis-GenPois_AR1_NoReg", {

  # set parameters
  CountDist      = "Generalized Poisson"
  alpha          = 1
  mu             = 3
  MargParm       = c(alpha,mu)
  ARParm         = 0.75
  MAParm         = NULL
  ARMAModel      = c(length(ARParm),length(MAParm))
  Task           = "Synthesis"
  TrueParam      = c(MargParm,ARParm, MAParm)
  SampleSize     = 50
  Intercept      = FALSE

  # simulate data with the old lgc function
  set.seed(1)
  x1   = sim_lgc_old(SampleSize, CountDist, MargParm, ARParm, MAParm, Regressor = NULL, Intercept)

  # synthesize the data using the new sim_lgc function
  set.seed(1)
  x2  = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm, Regressor = NULL)

  expect_identical(x1, x2)


  })

test_that("New synthesis-GenPois_AR1_Reg", {

  # set parameters
  CountDist      = "Generalized Poisson"
  alpha          = 1
  b0             = 0.5
  b1             = 2
  MargParm       = c(b0,b1,alpha)
  ARParm         = 0.75
  MAParm         = NULL
  ARMAModel      = c(length(ARParm),length(MAParm))
  Task           = "Synthesis"
  TrueParam      = c(MargParm,ARParm, MAParm)
  SampleSize     = 50
  e              = rbinom(SampleSize,1,0.1)
  #Regressor      = cbind(rep(1,SampleSize),e)
  Regressor      = e
  Intercept      = TRUE

  # simulate data with the old lgc function
  set.seed(1)
  x1   = sim_lgc_old(SampleSize, CountDist, MargParm, ARParm, MAParm, Regressor,Intercept)


  # synthesize the data using the new sim_lgc function
  set.seed(1)
  x2  = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm, Regressor,Intercept)

  expect_identical(x1, x2)



})

