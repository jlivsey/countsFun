
test_that("VGAM versus our implementation of GenPois", {

# Specify model and methods
n              = 100
CountDist      = "Generalized Poisson"
alpha          = 1
mu             = 3
MargParm       = c(alpha,mu)
ARParm         = 0.75
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))

# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm)

# Run the wrapper
mod = ModelScheme(DependentVar   = DependentVar,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel)

theta1 = InitialEstimates(mod)
theta1
set.seed(1)
a1 = ParticleFilter_Res_ARMA(theta1, mod)

CountDist      = "Generalized Poisson 2"

# Run the wrapper
mod = ModelScheme(DependentVar   = DependentVar,
                  CountDist      = CountDist,
                  ARMAModel      = ARMAModel)


theta2 = InitialEstimates(mod)
theta2
set.seed(1)
a2 = ParticleFilter_Res_ARMA(theta2, mod)

set.seed(1)
# Run the wrapper
mylgc = lgc(DependentVar   = DependentVar,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel)


# check the likelihoods
expect_equal(a1,a2,tolerance = 10^(-10))

# check the likelihoods (directly) and throughthe lgc wrapper
expect_equal(a1,as.numeric(mylgc$FitStatistics[1]),tolerance = 10^(-10))

# check the initial estimates
expect_equal(theta1[1],theta2[1],tolerance = 10^(-10))
expect_equal(theta1[2],theta2[2],tolerance = 10^(-10))
expect_equal(theta1[3],theta2[3],tolerance = 10^(-10))



})
