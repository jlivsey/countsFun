#======================================================================================================#
#Purpose:   Evaluate the Dominick data while providing fewer arguments.
#
# Author:   Stefanos Kechagias
# Date:     August 2024
#=====================================================================================================#

test_that("LGC wrapper works for evaluation MixPois AR(1)", {


# load the data
data(drinksales)

# attach the dataframe
n = 104
Smallsales  = drinksales[1:n,]

# regressor variable with intercept
DependentVar   = Smallsales$MOVE
Regressor      = Smallsales$Buy
CountDist      = "Negative Binomial"
ARMAModel      = c(2,0)
OptMethod      = "L-BFGS-B"
initialParam   = c(2.1756853 , 1.2048704,0.5, -0.3875602, 0.0603419 )
verbose        = FALSE

# save the data in a data frame
df = data.frame(DependentVar, Regressor)

# specify the regression model
RegModel = DependentVar~Regressor

# call the wrapper function with less arguments
mylgc = lgc(RegModel       = RegModel,
            df             = df,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            OptMethod      = OptMethod,
         initialParam      = initialParam,
         verbose           = verbose)

expect_equal(mylgc$FitStatistics[[1]], 392.673, tolerance = 10^(-3))

})
