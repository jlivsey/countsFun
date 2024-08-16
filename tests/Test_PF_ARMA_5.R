# Show that mixed poisson doesnt have an issue with L-BFGS-B - things
# seem to break with bobyqa
n              = 200
Regressor      = NULL
CountDist      = "Mixed Poisson"
MargParm       = c(2, 5, 0.7)
ARParm         = 0.5
MAParm         = NULL
ARMAModel      = c(length(ARParm),length(MAParm))
ParticleNumber = 5
epsilon        = 0.5
EstMethod      = "PFR"
TrueParam      = NULL
initialParam   = NULL
Task           = 'Optimization'
SampleSize     = NULL
nsim           = NULL
no_cores       = NULL
OptMethod      = "L-BFGS-B"
OutputType     = "list"
ParamScheme    = NULL
maxdiff        = 10^(-8)

# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm)


mod = ModelScheme(DependentVar, Regressor, EstMethod, ARMAModel, CountDist,ParticleNumber, epsilon,
                  initialParam, TrueParam, Task,SampleSize, OptMethod, OutputType, ParamScheme, maxdiff)


if (is.null(initialParam)){
  theta  = InitEst = InitialEstimates(mod)
  mod$initialParam = InitEst
}else{
  theta  = InitEst = mod$initialParam
}


# Run the wrapper
mylgc = lgc(DependentVar   = DependentVar,
            Regressor      = Regressor,
            EstMethod      = EstMethod,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            ParticleNumber = ParticleNumber,
            epsilon        = epsilon,
            initialParam   = initialParam,
            TrueParam      = TrueParam,
            Task           = Task,
            SampleSize     = SampleSize,
            nsim           = nsim,
            no_cores       = no_cores,
            OptMethod      = OptMethod,
            OutputType     = OutputType,
            ParamScheme    = ParamScheme,
            maxdiff        = maxdiff)
mylgc
