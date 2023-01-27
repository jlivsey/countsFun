#-------------------------------------------------------------------------------------------------------------#
# PURPOSE: Compare the derivative of the incomplete Gamma function computed using automatic differentiation
#          versus the standard numerical differentiation computation.
#
#
#
# Notes:   The reason we want to do this is that the incomplete gamma function is the main ingredient in the
#          poisson cdf. If we get the derivative of the incomplete gamma function then we can easily get the
#          derivative of the poisson cdf. In addition the automatic differentiation is is significantly faster
#          so this will be a great speed boost for us.
#
#-------------------------------------------------------------------------------------------------------------#
source("C:/Users/statha/OneDrive - SAS/Documents/countsFun/tests/mygammainc.R")
library(numDeriv)
# fix a point x and take several values of a
x = 4
# a = 3.1
a = seq(2,4.2,0.0005)
x_d = 0
a_d = 1

# allocate memory for the derivatives w.r.t the first variable
AutoDiff_1 = rep(0,length(a))
NumDeriv_1 = rep(0,length(a))
Diff_1 = rep(0,length(a))

# allocate memory for the derivatives w.r.t the second variable
AutoDiff_2 = rep(0,length(a))
NumDeriv_2 = rep(0,length(a))
Diff_2 = rep(0,length(a))

# for a fixed x compute the derivative of the function for sevferla values of a
for(k in 1:length(a)){
  # using the automatic differentiation approach
  AutoDiff_1[k] = NewGammainc_d_2(x, a[k], 0, 1)[2]
  AutoDiff_2[k] = NewGammainc_d_2(x, a[k], 1, 0)[2]

  # using finite differences (grad is a function in the numDeriv package)
  X = c(x,a[k])
  NumDeriv_1[k] = grad(func = gammainc,x=X)[2]
  NumDeriv_2[k] = grad(func = gammainc,x=X)[1]

  # compute the difference of the derivatives for the two approaches
  Diff_1[k] = AutoDiff_1[k]-NumDeriv_1[k]
  Diff_2[k] = AutoDiff_2[k]-NumDeriv_2[k]
}

# plot the derivatives
par(mfrow=c(2,3))
plot(a,AutoDiff_1, type = "l")
plot(a,NumDeriv_1, type = "l")
plot(a,Diff_1, type = "l")
plot(a,AutoDiff_2, type = "l")
plot(a,NumDeriv_2, type = "l")
plot(a,Diff_2, type = "l")



