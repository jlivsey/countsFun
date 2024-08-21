test_that("Compute the Pred distribution for the Dominick Data", {


# regressor variable with intercept
CountDist      = "Poisson"
Regressor      = NULL
n              = 100
ARMAModel      = c(2,0)
MargParm       = 3
ARParm         = c(0.5, 0.2)
MAParm         = NULL
OptMethod      = "L-BFGS-B"
Task           = "Optimization"

# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

# # call the wrapper function with less arguments
# mylgc = lgc(DependentVar   = DependentVar,
#             Regressor      = Regressor,
#             CountDist      = CountDist,
#             ARMAModel      = ARMAModel,
#             OptMethod      = OptMethod,
#                  Task      = Task)


# call the wrapper function with less arguments
mod = ModelScheme(DependentVar   = DependentVar,
                  Regressor      = Regressor,
                  CountDist      = CountDist,
                  ARMAModel      = ARMAModel,
                  OptMethod      = OptMethod,
                       Task      = Task)

# theta = mylgc$ParamEstimates
theta = c( MargParm, ARParm, MAParm)

# compute predictive distribution
predDist = PDvalues(theta, mod)

H = 10
# compute PIT values
PIT = PITvalues(H, predDist)

# saved expected values
ExpValues = c(0.11315799, 0.11552394, 0.10897874, 0.11410315, 0.09598454, 0.07847884, 0.06713758,
0.07751391, 0.09185189, 0.13726942)

for (i in 1:10){
  expect_equal(PIT[i], ExpValues[i], tolerance = 10^(-5))
}


})
