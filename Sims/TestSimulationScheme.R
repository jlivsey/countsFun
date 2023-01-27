#-------------------------------------------------------------------------#
# Purpose: Fit Synthetic data using the lgc wrapper
#
#
# Date: Jan 2023
#-------------------------------------------------------------------------#

# load necessary libraries
library(optimx)
library(ltsa)
library(countsFun)
library(doParallel)
library(foreach)

# Specify model and methods
SampleSize     = 100
Regressor      = NULL
CountDist      = "Poisson"
MargParm       = 10
ARParm         = c(0.75, 0.2)
MAParm         = NULL
EstMethod      = "PFR"
ParticleNumber = 10
epsilon        = 0.5
TrueParam      = c(MargParm, ARParm, MAParm)
initialParam   = NULL
OptMethod      = "bobyqa"
ParamScheme    = 1
nsim           = 2
Task           = 'Simulation'
OutputType     = "list"
no_cores       = NULL
DependentVar   = NULL
ARMAModel      = c(length(ARParm), length(MAParm))

# Run the Simulation
mysim = lgc(DependentVar   = DependentVar,
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
            ParamScheme    = ParamScheme)
