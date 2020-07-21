#======================================================================================================#
#Purpose:   Load series selected from the Dominick retail data.
#           Look at different windows of data for correlation structure.
#
# Author:   James Livsey
# Team:     Vladas Pipiras, James Livsey, Stefanos Kechagias, Robert Lund, Yisu Jia
# Date:     June 2020
#=====================================================================================================#

# load libraries
library(gcmr)
#library(ggplot2)
library(tictoc)
library(lavaSearch2)
symmetrize <- lavaSearch2:::symmetrize


# load the data
#mysales = read.csv("C:/Users/Stef/Desktop/countsFun/data/MySelectedSeries.csv")
# mysales = read.csv("/Users/stef/Desktop/countsFun/data/MySelectedSeries.csv")
mysales = read.csv("~/github/countsFun/data/MySelectedSeries.csv")


# attach the datafrmae
attach(mysales)

# plot the series
plot.ts(MOVE[1:200], ylab = "sales", xlab = "time")

# use GCMR to fit the data
tic()
mod <- gcmr(MOVE[1:200]~Buy[1:200],
            marginal = negbin.marg,
            cormat = arma.cormat(3, 0),
            no.se = FALSE)
toc()

mod
sqrt(diag(solve(-1*mod$hessian)))

names(mod)

e <- residuals(mod)

acf(e)
pacf(e)

e2 <- residuals(lm(MOVE~Buy))
acf(e2)
pacf(e2)

library(tsbox)
ts_plot(cbind(ts(e), ts(e2)))


length(MOVE)

# break sample up from observations 30 thru 361 into 1st and 2nd part.
# Look lag1 autocorrelation of each part
lag1 <- rep(0, 391)
lag1.2 <- rep(0, 391)
for(i in 30:361){
  x <- MOVE[1:i]
  x2 <- MOVE[(i+1):391]
  lag1[i] <- acf(x, plot = FALSE, lag.max = 2)$acf[2]
  lag1.2[i] <- acf(x2, plot = FALSE, lag.max = 2)$acf[2]
}

plot(c(0, 400), range(lag1, lag1.2), type = 'n', main = "original data")
lines(lag1)
lines(lag1.2, col = 2)
abline(h = 0 , lty = 'dotted')


# break sample up into moving window of length WinSize
# Look lag1 autocorrelation
WinSize <- 100
lag1 <- rep(0, 391)
for(i in 1:(391-WinSize)){
  x <- MOVE[i:(i+WinSize)]
  lag1[i] <- acf(x, plot = FALSE, lag.max = 2)$acf[2]
}
plot.ts(lag1)
abline(h = 0, lty = 'dotted')


# residuals - break sample up from observations 30 thru 361 into 1st and 2nd part.
# Look lag1 autocorrelation of each part
lag1 <- rep(0, 391)
lag1.2 <- rep(0, 391)
for(i in 30:361){
  x <- residuals(lm(MOVE[1:i]~Buy[1:i]))
  x2 <- residuals(lm(MOVE[(i+1):391]~Buy[(i+1):391]))
  lag1[i] <- acf(x, plot = FALSE, lag.max = 2)$acf[2]
  lag1.2[i] <- acf(x2, plot = FALSE, lag.max = 2)$acf[2]
}

plot(c(0, 400), range(lag1, lag1.2), type = 'n', main = "residuals")
lines(lag1)
lines(lag1.2, col = 2)
abline(h = 0 , lty = 'dotted')



# ADF test of residuals - break sample up from observations 30 thru 361 into 1st and 2nd part.
# Look lag1 autocorrelation of each part
adf1 <- rep(0, 391)
adf1.2 <- rep(0, 391)

for(i in 30:361){
  x <- residuals(lm(MOVE[1:i]~Buy[1:i]))
  x2 <- residuals(lm(MOVE[(i+1):391]~Buy[(i+1):391]))
  adf1[i] <- adf.test(x)$p.value
  adf1.2[i] <- adf.test(x2)$p.value
}

plot(c(0, 400), range(adf1, adf1.2), type = 'n', main = "residuals")
lines(adf1)
# lines(adf1.2, col = 2)
abline(h = 0 , lty = 'dotted')

ts.plot(adf1)
abline(h = .05, lty = 'dotted')






