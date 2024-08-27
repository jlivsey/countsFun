#library(tidyverse)

# Set Sample size
n = 50

# Generate a regressor that will be used as a linear predictor
Regressor = cbind(1,rnorm(n,0,1))

# select parameters for linear predictor
b0  = 1.5
b1  = 3
beta  = c(b0,b1)

# set the constant binomial parameter
r   = 4

CountDist = "Binomial"
MargParm = c(b0,b1)
ARParm   = 0.5
MAParm   = NULL
ntrials  = r

# simulate the data with time dependene
y = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor, ntrials)


## compute the probabilities so that logit(p_t) = bo + b1*Regressor
# p_t = plogis(Regressor%*%beta)
## simulate the data (no time dependence)
# y = rbinom(n = N, size = r, prob = p_t)

# fit the model
glmfit = glm(cbind(y,ntrials-y) ~ Regressor[,2] , family = 'binomial')

# check the coefficients
coef(glmfit)

