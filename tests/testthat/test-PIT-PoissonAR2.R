test_that("Compute the Pred distribution for the Dominick Data", {


# regressor variable with intercept
CountDist      = "Poisson"
Regressor      = NULL
Intercept      = NULL
n              = 100
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

# saved expected values
ExpValues = c(0.11387452, 0.12087723, 0.11022829, 0.10980887, 0.09137042, 0.07570144, 0.06788390, 0.07421936, 0.09515664, 0.14087933)

for (i in 1:10){
  expect_equal(PIT[i], ExpValues[i], tolerance = 10^(-5))
}


})
