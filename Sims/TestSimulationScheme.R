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
ARMAorder      = c(length(ARParm),length(MAParm))
ParticleNumber = 10
epsilon        = 0.5
initialParam   = c(MargParm, ARParm, MAParm)
# initialParam   = NULL
theta          = initialParam
OptMethod      = "bobyqa"
ParamScheme    = 1
nsim           = 1


no_cores      = detectCores() - 1                  # Select the number of cores
##-------------------------------------------------------------------------------------------------#

# generate all the realizations and save in a list
DependentVar <- list()
for (i in 1:nsim) {
  set.seed(i)
  DependentVar[[i]] = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm)
}

ARMAorder = c(length(ARParm), length(MAParm))

# initiate and register the cluster
cl <- makeCluster(no_cores)

#clusterSetRNGStream(cl, 1001) #make the bootstrapping exactly the same as above to equate computation time
registerDoParallel(cl)

# run foreach
all = foreach(index = 1:nsim,
              .packages=c('countsFun', 'optimx', 'itsmr', 'ltsa'))  %dopar%  {
                lgc( DependentVar[[index]],
                    Regressor,
                    EstMethod,
                    CountDist,
                    ARMAorder,
                    ParticleNumber,
                    epsilon,
                    initialParam,
                    OptMethod
                )
              }

stopCluster(cl)





# mysim = LGCSimulation(nsim           = nsim,
#                       SampleSize     = SampleSize,
#                       Regressor      = Regressor,
#                       EstMethod      = EstMethod,
#                       CountDist      = CountDist,
#                       MargParm       = MargParm,
#                       ARParm         = ARParm,
#                       MAParm         = MAParm,
#                       ParticleNumber = ParticleNumber,
#                       epsilon        = epsilon,
#                       initialParam   = initialParam,
#                       OptMethod      = OptMethod,
#                       ParamScheme    = ParamScheme
#                       )
#
#
#
#
