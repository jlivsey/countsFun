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
#library(ggplot2)
library(tictoc)
library(lavaSearch2)
symmetrize <- lavaSearch2:::symmetrize


# load the data
#mysales = read.csv("C:/Users/Stef/Desktop/countsFun/data/MySelectedSeries.csv")
mysales = read.csv("/Users/stef/Desktop/countsFun/data/MySelectedSeries.csv")


# attach the datafrmae
attach(mysales)

# plot the series
#plot.ts(MOVE, ylab = "sales", xlab = "time")

# use GCMR to fit the data
# tic()
# mod <- gcmr(MOVE~Buy, marginal = negbin.marg, cormat = arma.cormat(3, 0), no.se = FALSE)
# toc()

# plot residuals
#plot(mod)


b0 = 2.390835
b1 = 0.683351
b = c(b0,b1)
k  = 1/4
phi = c(-0.1696836, 0.2796241, 0.2267073)
#(b0, b1, k, phi1, phi2, phi3])
LB = c(-Inf, -Inf, 0.001, -Inf, -Inf, -Inf)
UB = c(Inf, Inf, Inf, Inf, Inf, Inf)
MaxCdf = 5000
nHC = 30
theta = c(b,k,phi)
Regressor = cbind(rep(1,length(Buy)),Buy)
ARMAorder = c(3,0)
theta = c(2.9713351, 1.3705352,0.6812780,-0.1979509,0.3350090,  0.2572088  )
data = MOVE


# tic()
# GaussLogLikNB_Reg(theta, MOVE, Regressor, ARMAorder, MaxCdf, nHC)
# toc()

X = MOVE
x0 = theta
ARMAorder = c(3,0)

optim.output <- optim(par       = x0,
                      fn        = GaussLogLikNB_Reg,
                      data      = X,
                      Regressor = Regressor,
                      ARMAorder = ARMAorder,
                      MaxCdf    = MaxCdf,
                      nHC       = nHC,
                      method    = "L-BFGS-B",
                      hessian   = TRUE,
                      lower     = LB,
                      upper     = UB
)

# tic()
# mod2 = FitGaussianLikNB_Reg(x0, X, Regressor, LB, UB, ARMAorder, MaxCdf, nHC)
# toc()


