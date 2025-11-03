test_that("Misspecification truncnorm used", {

  # purpose: this is an example where the truncnorm library will be used when sampling particles.
  # see #27 and commit on Sep 27, 2025.

# Generate Mixed Poisson model
set.seed(2)
n              = 50
CountDist      = "Mixed Poisson"
MargParm       = c(2, 5, 0.7)
ARParm         = 0.75
MAParm         = NULL
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm)

# generate a zero inflated (p =0.9) poisson 50 regressor
Regressor      = ifelse(rbinom(n, size = 1, prob = 0.9) > 0, 0, rpois(n, lambda = 50))
DependentVar    = DependentVar + Regressor
#plot(1:n,DependentVar, type="l")

# Set parameters for model to fit
CountDist      = "Poisson"
MargParm       = 55
ARParm         = 0.75
MAParm         = -0.3
ARMAModel      = c(length(ARParm),length(MAParm))
Task           = 'Optimization'
OptMethod      = 'L-BFGS-B'
initialParam   = c(MargParm,ARParm,MAParm)
theta          = initialParam
verbose        = TRUE
# call the wrapper function with less arguments
mod = ModelScheme(DependentVar   = DependentVar,
                  Regressor      = NULL,
                  CountDist      = CountDist,
                  ARMAModel      = ARMAModel,
                  OptMethod      = OptMethod,
                  initialParam   = initialParam,
                  verbose        = verbose)

theta = c(54.99986, 0.8015283, -0.2920431)
a1 = ParticleFilter_Res_ARMA(theta,mod)
a2 = ParticleFilter(theta,mod)

expect_equal(a1,631.808,tolerance=10^(-3))
expect_equal(a2,631.6743,tolerance=10^(-3))


})
