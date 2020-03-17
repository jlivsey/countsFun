# libraries
library(countsFun)
library(tictoc)
library(georob)
library(lavaSearch2)
symmetrize <- lavaSearch2:::symmetrize
#source('C:/Users/stefa/Desktop/countsFun/R/negBinom-functions.R')

# sample size
n = 200

# generate a regressor dummy variable
X  = cbind(rep(1,n),rbinom(n, 1 ,0.4))

# specify parameters
r = 4
b0 = 1
b1 = 2
b = rbind(b0,b1)
m = exp(X%*%b)
phi = 0.4

p = m/(m+r)


tic()
G5=CovarNegBinAR_Reg(n, r, m, phi)
toc()

tic()
G = CovarNegBinAR(n,r, p[2], phi)
toc()
