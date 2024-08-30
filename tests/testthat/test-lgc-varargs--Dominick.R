#======================================================================================================#
#Purpose:   Evaluate the Dominick data while providing fewer arguments.
#
# Author:   Stefanos Kechagias
# Date:     August 2024
#=====================================================================================================#

test_that("LGC wrapper works for evaluation MixPois AR(1)", {

# load libraries
library(optimx)
library(ltsa)
require(countsFun)
library(itsmr)
library(tictoc)
library(devtools)
library(VGAM)
library(iZID)

# load the data
#mysales = read.csv("https://raw.githubusercontent.com/jlivsey/countsFun/master/data/MySelectedSeries.csv")
data(MySelectedSeries)

# attach the dataframe
n = 104
Smallsales  = mysales[1:n,]

# regressor variable with intercept
DependentVar   = Smallsales$MOVE
Regressor      = cbind(rep(1,length(Smallsales$Buy)),Smallsales$Buy)
CountDist      = "Negative Binomial"
ARMAModel      = c(2,0)
OptMethod      = "L-BFGS-B"
initialParam   = c(2.1756853 , 1.2048704,0.5, -0.3875602, 0.0603419 )


# call the wrapper function with less arguments
mylgc = lgc(DependentVar   = DependentVar,
            Regressor      = Regressor,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            OptMethod      = OptMethod,
              initialParam = initialParam)

expect_equal(mylgc$FitStatistics[[1]], 392.673, tolerance = 10^(-3))

})
