
test_that("The likelihood function ParticleFilter_Res_ARMA works for a 'misspecified'
          model on the Dominick data.", {

# load libraries
#library(optimx)
#library(ltsa)
#require(countsFun)
#library(itsmr)
#library(tictoc)

# load the data
#mysales = read.csv("https://raw.githubusercontent.com/jlivsey/countsFun/master/data/MySelectedSeries.csv")
data(MySelectedSeries)

# deal with a smaller sample size that we considered in the JASA paper
n = 104
Smallsales  = mysales[1:n,]

# regressor variable with intercept
DependentVar   = Smallsales$MOVE
Regressor      = cbind(rep(1,length(Smallsales$Buy)),Smallsales$Buy)
CountDist      = "Poisson"
ARMAModel      = c(2,0)
OptMethod      = "L-BFGS-B"
#initialParam   = c(2.1756853 , 1.2048704, -0.3875602, 0.0603419 )
initialParam   = c(10, -0.3875602, 0.0603419 )

# call the wrapper function with less arguments
mod = ModelScheme(DependentVar   = DependentVar,
                  Regressor      = NULL,
                  CountDist      = CountDist,
                  ARMAModel      = ARMAModel,
                  OptMethod      = OptMethod,
                  initialParam = initialParam)

# get initial estimate
if (is.null(mod$initialParam)){
  mod$initialParam = InitialEstimates(mod)
}
theta = mod$initialParam

set.seed(1)
expect_equal(ParticleFilter_Res_ARMA(theta,mod), 100000000, tolerance = 10^(-5))
# expect_equal(ParticleFilter_Res_ARMA_MisSpec(theta,mod), 1680.537826, tolerance = 10^(-5))


})
