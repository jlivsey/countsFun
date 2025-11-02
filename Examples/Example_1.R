#=====================================================================================#
# Purpose: Write a simple example of generating and fitting a Poisson-AR(1) model
#
# Date:   October 2025.
#=====================================================================================#

# specify model to synthesize from
SampleSize     = 200
CountDist      = "Poisson"
MargParm       = 3
ARParm         = 0.5
MAParm         = NULL

# synthesize data
set.seed(2)
DependentVar   = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm)
data           = data.frame(DependentVar)

# specify model to be fitted
formula        = DependentVar~0
ARMAModel      = c(length(ARParm),length(MAParm))

# specify Task and and optimization parameters
Task           = 'Optimization'
initialParam   = c(5,0.1)
OptMethod      = "L-BFGS-B"

# run the main wrapper
fit = lgc(
       formula = formula,
          data = data,
     CountDist = CountDist,
     ARMAModel = ARMAModel,
  initialParam = initialParam,
          Task = Task,
     OptMethod = OptMethod
)

# print the result
print(fit)
