#======================================================================================#
# Purpose: Test NegBin-White Noise fit: GLM vs LGC
#
# Author: Stefanos Kechagias
#======================================================================================#

#--------------------Generate using base R function--------------
test_that("NegBin fit GLM vs LGC", {

set.seed(123)

# True parameters
n  <- 200
b0 <- 1.0     # intercept
b1 <- 0.5     # slope
theta_true <- 2  # dispersion (size parameter)


# Generate regressor
x <- rnorm(n, mean = 0, sd = 1)

# Generate mean
mu <- exp(b0 + b1 * x)

# Generate response (NB)
y <- rnbinom(n, size = theta_true, mu = mu)

# Fit model using GLM
fit1 <- glm.nb(y ~ x)
glmESts = c(fit1$coefficients, fit1$theta)

# fit model using lgc starting from a point away from the truth
RegModel       = y ~ 1 + x
df = data.frame(y,x)
CountDist = "Negative Binomial"
ARMAModel = c(0,0)
Task = 'Optimization'
initialParam   = c(0,1,1)

# populate a list with the model characteristics to compute initial estimates
mod_test = ModelSpec(RegModel = RegModel,
                df = df,
                CountDist = CountDist,
                ARMAModel = ARMAModel)

InitEsts = InitialEstimates(mod_test)

# fit  model using lgc
fit = lgc(RegModel = RegModel,
          df = df,
          CountDist = CountDist,
          ARMAModel = ARMAModel,
          initialParam = initialParam,
          Task = Task)

PFEsts = coefficients(fit)

# compare LGC initial with GLM - should be identical
expect_equal(as.numeric(glmESts), as.numeric(InitEsts))

# compare LGC with GLM
expect_equal(as.numeric(glmESts), as.numeric(PFEsts), tolerance = 0.00001)

})

#--------------------Generate using sim_lgc function--------------
