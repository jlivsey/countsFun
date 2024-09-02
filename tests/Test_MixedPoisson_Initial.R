#======================================================================================#
# Purpose: Get initital estimates for mixed poisson

#======================================================================================#

# set linear predictor parameters
b0_1 = 2
b1_1 = 0.5

b0_2 = 3
b1_2 =  1

# set mixing probability
prob = 0.4

# gather Marginal Parameters
MargParm = c(b0_1, b1_1, b0_2, b1_2,prob)

# set other model details
n = 50
CountDist = "Mixed Poisson"
ARParm = NULL
MAParm = NULL
ARMAModel      = c(length(ARParm),length(MAParm))

# generate regressor
x = runif(n)
Regressor = cbind(1,x)

# coimpute the Poisson means
lambda1 = exp(Regressor%*%c(b0_1,b1_1))
lambda2 = exp(Regressor%*%c(b0_2,b1_2))

# generate mixing Bernoulli variable
w = rbinom(n,1,prob)

# generate Poisson mixture
Comp1 = rpois(n,lambda1)
Comp2 = rpois(n,lambda2)
y = w*Comp1 + (1-w)*Comp2

# gather the model details
mod = ModelScheme(DependentVar   = y,
                  CountDist      = CountDist,
                  ARMAModel      = ARMAModel,
                  Regressor      = Regressor)

# compute Initial Estimates
theta = InitialEstimates(mod)
theta

# compute likelihood
ParticleFilter_Res_ARMA(theta, mod)
