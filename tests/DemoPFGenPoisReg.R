# libraries
library(countsFun)
library(tictoc)
library(lavaSearch2)
library(optimx)
library(FitAR)

symmetrize <- lavaSearch2:::symmetrize
# sample size
n = 200

# generate a regressor dummy variable
Regressor  = cbind(rep(1,n),rbinom(n, 1 ,0.4))

# specify parameters
b0 = 1
b1 = 2
b = rbind(b0,b1)
lambda = exp(Regressor%*%b)
omega = 0.5
phi = -0.5
MaxCdf = 100
nHC = 30

# list with true ARMA parameters
ARMAmodel = list(NULL,NULL)
ARMAmodel[[1]] = phi
ARMAorder = c(1,0)
ParticleNumber = 10
CountDist = "Generalized Poisson"
epsilon=1
# parameters to be passed in likelihood
theta = c(b,omega,phi)

# generate Gen Pois series
data = sim_genpois(n, ARMAmodel, omega, lambda)

tic()
a = ParticleFilterRes_Reg(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon)
toc()
