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
ARParm         = c(0.75)
MAParm         = NULL
EstMethod      = "PFR"
ParticleNumber = 10
epsilon        = 0.5
#initialParam   = c(MargParm, ARParm, MAParm)
initialParam   = NULL
OptMethod      = "bobyqa"
ParamScheme    = 1
nsim           = 2
Task           = 'Optimization'
OutputType     = "data.frame"
no_cores       = NULL

# Run the Simulation
mysim = LGCSimulation(nsim           = nsim,
                      SampleSize     = SampleSize,
                      Regressor      = Regressor,
                      EstMethod      = EstMethod,
                      CountDist      = CountDist,
                      MargParm       = MargParm,
                      ARParm         = ARParm,
                      MAParm         = MAParm,
                      ParticleNumber = ParticleNumber,
                      epsilon        = epsilon,
                      initialParam   = initialParam,
                      OptMethod      = OptMethod,
                      Task           = Task,
                      OutputType     = OutputType,
                      ParamScheme    = ParamScheme,
                      no_cores       = no_cores
                      )
