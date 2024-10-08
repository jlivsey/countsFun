
#================================================================================#
# PURPOSE: Compute Predictive Distribution and PIT for Dominick data
# DATE: July 2020
# NOTES: see relation (36) in JASA paper
# Author: Stefanos Kechagias
#================================================================================#

test_that("Compute the Pred distribution for the Dominick Data", {

library(itsmr)
library(countsFun)

# load the data
data(MySelectedSeries)

# attach the datafrmae
n = 104
Smallsales  = mysales[1:n,]
MOVE = Smallsales$MOVE
Buy = Smallsales$Buy


# regressor variable with intercept
DependentVar   = MOVE
Regressor      = cbind(rep(1,length(Smallsales$Buy)),Smallsales$Buy)
CountDist      = "Negative Binomial"
ARMAModel      = c(2,0)
OptMethod      = "L-BFGS-B"
initialParam   = c(2.1756853 , 1.2048704,0.5, -0.3875602, 0.0603419 )

# call the wrapper function with less arguments
mod = ModelScheme(DependentVar   = DependentVar,
            Regressor      = Regressor,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            OptMethod      = OptMethod,
            initialParam = initialParam)


theta = initialParam

# compute predictive distribution
predDist = PDvalues(theta, mod)

# select two times and also the sum
expect_equal(predDist[1,30], 0.5448965, tolerance = 10^(-5))
expect_equal(predDist[2,30], 0.5594144, tolerance = 10^(-5))
expect_equal(predDist[1,60], 0.3794839, tolerance = 10^(-5))
expect_equal(predDist[2,60], 0.5259663, tolerance = 10^(-5))
expect_equal(sum(predDist), 82.91289, tolerance = 10^(-4))

})

