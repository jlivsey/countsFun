
test_that("The ParticleFilter_Res_ARMA is the same as ParticleFilter", {

# Specify model and methods
n              = 200
# Regressor      = cbind(rep(1,n),rbinom(n,1,0.25))
Regressor      = rnorm(n)
Intercept      = FALSE
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
Task           = 'Evaluation'
OptMethod      = "bobyqa"
OutputType     = "data.frame"
maxdiff        = 10^(-6)

# specify the regression formula (no regressors here)
RegModel     = DependentVar ~ 0

# simulate data
set.seed(2)
DependentVar = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, RegModel)
df         = data.frame(DependentVar)

# Parse the model
mod = ModelScheme(RegModel = RegModel,
                        df = df,
                 CountDist = CountDist,
                 EstMethod = EstMethod,
                 ARMAModel = ARMAModel,
            ParticleNumber = ParticleNumber,
                   epsilon = epsilon,
              initialParam = initialParam,
                 TrueParam = TrueParam,
                      Task = Task,
                SampleSize = SampleSize,
                 OptMethod = OptMethod,
                OutputType = OutputType,
                   maxdiff = maxdiff)

# evaluate the likelihood
a1 = ParticleFilter_Res_ARMA(initialParam,mod)
a2 = ParticleFilter(initialParam,mod)
expect_equal(a1,a2, tolerance = 10^(-4))

})


