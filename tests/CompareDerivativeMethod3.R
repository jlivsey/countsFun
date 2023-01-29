#-------------------------------------------------------------------------------------------------------------#
# PURPOSE: Compare the derivative of the poisson cdf computed using automatic differentiation
#          versus the standard numerical differentiation computation. See also CompareDerivativeMethod.R
#
# N
#-------------------------------------------------------------------------------------------------------------#
source("C:/Users/statha/OneDrive - SAS/Documents/countsFun/tests/mygammainc.R")
library(numDeriv)
library(tictoc)

#-------------------------------------------- TEST 1  -----------------------------------------------
# check that the function myppois returns the same value as the ppois in R
x      = 2
lambda = 4
nsim   = 100000
X      = c(x,lambda)

myppois(x,lambda)
c(ppois(x,lambda), grad(function(X){ppois(X[1],X[2])},X)[2])


#-------------------------------------------- TEST 2  -----------------------------------------------
# compare speeds of myppois versus ppois
t1 = tic()
for(i in 1:nsim){
  a1   = ppois(X[1],X[2])
  a1_d = grad(function(X){ppois(X[1],X[2])},X)[2]
}
t1 = tic() - t1


t2 = tic()
for(i in 1:nsim){
  b1_d = myppois(x,lambda)[2]
}
t2 = tic() - t2

sprintf("Finite Differences derivative of cdf is equal to %.5f takes %.2f seconds",b1_d, t1)
sprintf("Autodiff derivative of cdf is equal to %.5f takes %.5f seconds",a1_d, t2)



#-------------------------------------------- TEST 3  -----------------------------------------------
# fix a point x and take several values of lambda
x      = 4
lambda = seq(2,17,0.005)


# allocate memory for the derivatives w.r.t the first variable
MyPoisCDF_Deriv = rep(0,length(lambda))
RPoisCDF_Deriv  = rep(0,length(lambda))
Diff            = rep(0,length(lambda))

# for lambda fixed x compute the derivative of the function for sevferla values of lambda
for(k in 1:length(lambda)){
  # using the automatic differentiation approach
  MyPoisCDF_Deriv[k] = myppois(x, lambda[k])[2]

  # using finite differences (grad is lambda function in the numDeriv package)
  X = c(x,lambda[k])
  RPoisCDF_Deriv[k] = grad(function(X){ppois(X[1],X[2])},X)[2]

  # compute the difference of the derivatives for the two approaches
  Diff[k] = MyPoisCDF_Deriv[k]-RPoisCDF_Deriv[k]
}

# plot the derivatives
par(mfrow=c(1,3))
plot(lambda,MyPoisCDF_Deriv, type = "l")
plot(lambda,RPoisCDF_Deriv, type = "l")
plot(lambda,Diff, type = "l")





