#---------------------------------------------------------------------------------#
# Optim works but my lgc does not
#---------------------------------------------------------------------------------#


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
n              = 50
CountDist      = "Poisson"
MargParm       = 55
ARParm         = 0.75
MAParm         = -0.3
ARMAModel      = c(length(ARParm),length(MAParm))
Task           = 'Optimization'
OptMethod      = 'L-BFGS-B'
initialParam   = c(MargParm,ARParm,MAParm)
theta          = initialParam

# call the wrapper function with less arguments
mod = ModelScheme(DependentVar   = DependentVar,
                  Regressor      = NULL,
                  CountDist      = CountDist,
                  ARMAModel      = ARMAModel,
                  OptMethod      = OptMethod,
                  initialParam = initialParam)


# optim.output <- optimx(
#   par     = theta,
#   fn      = ParticleFilter_Res_ARMA,
#   lower   = mod$LB,
#   upper   = mod$UB,
#   hessian = TRUE,
#   method  = mod$OptMethod,
#   mod     = mod)


mylgc = lgc(DependentVar   = DependentVar,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            initialParam   = initialParam,
            Task           = Task)
