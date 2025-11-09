
test_that("The likelihood function ParticleFilter_Res_ARMA works for a 'misspecified'
          model on the Dominick data.", {

# load the data
data(drinksales)

# deal with a smaller sample size that we considered in the JASA paper
n = 104
Smallsales  = drinksales[1:n,]

# regressor variable with intercept
DependentVar   = Smallsales$MOVE
df             = data.frame(DependentVar)
CountDist      = "Poisson"
ARMAModel      = c(2,0)
#initialParam   = c(2.1756853 , 1.2048704, -0.3875602, 0.0603419 )
initialParam   = c(10, -0.3875602, 0.0603419 )
verbose        = FALSE

# specify the regression formula (no regressors here - although EDA susggest I should use some)
RegModel     = DependentVar ~ 0

# Run the wrapper
mod = ModelScheme(RegModel = RegModel,
                        df = df,
                 CountDist = CountDist,
                 ARMAModel = ARMAModel,
              initialParam = initialParam,
                   verbose = verbose)

# get initial estimate
if (is.null(mod$initialParam)){
  mod$initialParam = InitialEstimates(mod)
}
theta = mod$initialParam

set.seed(1)
expect_equal(ParticleFilter_Res_ARMA(theta,mod), 100000000, tolerance = 10^(-5))

# after changing the function to compute limits ion log space when C_1 = 1 the value changes a
# bit
#expect_equal(ParticleFilter(theta,mod), 1442.823, tolerance = 10^(-5))
expect_equal(ParticleFilter(theta,mod), 1442.785, tolerance = 10^(-5))

})
