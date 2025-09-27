
#================================================================================#
# PURPOSE: Compute Predictive Distribution and PIT for Dominick data
# DATE: July 2020
# NOTES: see relation (36) in JASA paper
# Author: Stefanos Kechagias
#================================================================================#

test_that("Compute the Pred distribution for the Dominick Data", {


set.seed(1)
# load the data
data(drinksales)

# attach the datafrmae
n = 104
Smallsales  = drinksales[1:n,]
MOVE = Smallsales$MOVE
Buy = Smallsales$Buy


# regressor variable with intercept
DependentVar   = MOVE
Regressor      = Smallsales$Buy
Intercept      = TRUE
CountDist      = "Negative Binomial"
ARMAModel      = c(2,0)
OptMethod      = "L-BFGS-B"
initialParam   = c(2.1756853 , 1.2048704,0.5, -0.3875602, 0.0603419 )

# call the wrapper function with less arguments
mod = ModelScheme(DependentVar   = DependentVar,
            Regressor      = Regressor,
            Intercept      = Intercept,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            OptMethod      = OptMethod,
            initialParam = initialParam)


theta = initialParam

# compute predictive distribution
predDist = PDvalues(theta, mod)

# select two times and also the sum
expect_equal(predDist[1,30], 0.5402, tolerance = 10^(-4))
expect_equal(predDist[2,30], 0.5548, tolerance = 10^(-4))
expect_equal(predDist[1,60], 0.3811, tolerance = 10^(-4))
expect_equal(predDist[2,60], 0.5277, tolerance = 10^(-4))
expect_equal(sum(predDist), 82.914, tolerance = 10^(-3))

})

