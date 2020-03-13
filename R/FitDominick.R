#======================================================================================================#
#Purpose:   Load a series from the Dominick retail data that satisfies
#           1) Has negative lag 1 acf
#           2) Has small mean (<30)
#           3) Is cross correlated with a Sale dummy variable
#           4) After a regression with the dummy variable the residuals have significant lag1 acf.
#
# Notes:    The following steps briefly describe the procedure we followed
#
#           Step 1: Filter intermittent and retired series
#           Step 2: Focus on series with significant lag 1 negative acf
#           Step 3: Select series with mean < 30
#           Step 4: Run a regression on a Buy one get on free dummy variable and examine residuals
#           Step 5: Trial and error to select a series whose residuals have significant lag 1 neg acf
#
# Author:   Stefanos Kechagias
# Team:     Vladas Pipiras, James Livsey, Stefanos Kechagias, Robert Lund, Yisu Jia
# Date:     March 2020
#=====================================================================================================#

# load libraries
library(gcmr)
library(ggplot2)



# set the directory
#setwd("C:/Users/stefa/Dropbox/MVcopula/data/KaggleSalesTransactions")
setwd("C:/Users/statha/Desktop/Dominick")


# load the data
mysales = read.csv("C:/Users/statha/Desktop/Dominick/MySelectedSeries.csv")

# attach the datafrmae
attach(mysales)

# plot the series
#plot.ts(MOVE, ylab = "sales", xlab = "time")

# use GCMR to fit the data
#mod <- gcmr(MOVE~Buy, marginal = negbin.marg, cormat = arma.cormat(3, 0), no.se = FALSE)

# plot residuals
#plot(mod)


b0 = 2.390835
b1 = 0.683351
b = c(b0,b1)
r  = 0.593666
phi = c(-0.1696836, 0.2796241, 0.2267073)

theta = c(b,r,phi)
X = cbind(rep(1,length(Buy)),Buy)
FitGaussianLikNB_2(theta, MOVE, X, 3)

# X = cbind(rep(1,length(Buy)),Buy)
#
# theta0 = c(-0.417066419898571, -0.87759943351099, 9.0907108753284,
#            -0.20906710018397, 0.507553591887987, 0.188611090043294)
# GaussLogLikNB_2(theta0, MOVE, X,3 )





