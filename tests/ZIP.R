#===========================================================================================#
# PURPOSE: Demo of ZIP fit via glm packages VGAM and pscl
# Author: Stefanos Kechagias
#===========================================================================================#

# Set Sample Size
n = 200

# Set the inflation probability parameter
p = 0.85

# Set parameter for the Poisson link
b0 = 2
b1 = 4
beta = c(b0,b1)

# Generate a regressor
set.seed(1)
Regressor  = cbind(1,runif(n))

# Compute the Poisson mean using log link
lambda_t = exp(Regressor%*%beta)

# generate some data
x = rzipois(n, lambda_t, p)

# fit using the pscl package
m1 <- zeroinfl(x ~ Regressor[,2] | 1,)
coef(m1,matrix=TRUE)

# fit using VGAM
fit = vglm( x ~ Regressor[,2],zipoisson(zero=1))
coef(fit,matrix=TRUE)

