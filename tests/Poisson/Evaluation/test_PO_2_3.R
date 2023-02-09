# load necessary libraries
library(optimx)
library(ltsa)
library(countsFun)
library(itsmr)
library(tictoc)


# Specify model and methods
n              = 200
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 3
ARParm         = c(0.4, 0.4)
MAParm         = c(0.3, 0.5, -0.2, 0.1)
ARMAModel      = c(length(ARParm),length(MAParm))
ParticleNumber = 5
epsilon        = 0.5
EstMethod      = "PFR"
TrueParam      = NULL
initialParam   = NULL
Task           = 'Evaluation'
SampleSize     = NULL
nsim           = NULL
no_cores       = NULL
OptMethod      = "bobyqa"
OutputType     = "list"
ParamScheme    = NULL
maxdiff        = 10^(-8)


# simulate data
set.seed(2)
DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor)

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





