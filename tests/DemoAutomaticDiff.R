#-------------------------------------------------------------------------------------------
# PURPOSE: Demo for the speed difference of using automatic differentiation versus
# numerical differentiation.
#
# NOTES: See the dropbox notes AD-1 in ~\Dropbox\JSS Counts\AD notes
#
# STEPS: 1. Specify a function as we usually do
#        2. Specify a function using Automatic Differentation setup
#        3. Evaluate the two functions at the same point 1000 times and compare time
#
#
#-------------------------------------------------------------------------------------------
library(numDeriv)

# specify a function as we would usuall do
myfun = function(x){

  z = (sin(x[1]/x[2]) + x[1]/x[2] - exp(x[2]))*((x[1]/x[2]) - exp(x[2]))

  return(z)
}


# specify the function above in manner suitable for automatic differentiation
myfun_d = function(x1,x2, x1_d, x2_d){
  # The idea of automatic differentiation is to take each term of our function and name it
  # an intermediate variable vj. For example the term x[1]/x[2] above is denoted here as v1.
  # The arguments x1_d and x2_d take the values zero and 1 depending on which derivatiave
  # I am computing. If for example I call the function with (x1_d, x2_d) = (1,0) then
  # the derivative with respect to x1 will be computed. After each term v_j is computed I
  # will also compute the derivative of that term (and denote it as vj_d) using simple
  # derivative rules and chain rule. For example the term v2 = sin(v1) has derivative
  # v2_d = cos(v1_d)*v1_d.

  # start by "initializing the actual arguement iof the function using the "v" notation
  v_neg1   = x1
  v0       = x2
  v_neg1_d = x1_d
  v0_d     = x2_d

  # first term I will write using the "v" variables is the ratio x[1]/x[2]
  v1   = v_neg1/v0
  v1_d = (v0*v_neg1_d - v_neg1*v0_d)/v0^2

  # sin(x[1]/x[2])
  v2 = sin(v1)
  v2_d = cos(v1)*v1_d

  # exp(x[2]))
  v3 = exp(v0)
  v3_d = v3*v0_d

  # x[1]/x[2] - exp(x[2])
  v4 = v1 - v3
  v4_d = v1_d - v3_d

  # (sin(x[1]/x[2]) + x[1]/x[2] - exp(x[2]))
  v5 = v2+v4
  v5_d = v2_d +v4_d

  # (sin(x[1]/x[2]) + x[1]/x[2] - exp(x[2]))*((x[1]/x[2]) - exp(x[2]))
  v6 = v5*v4
  v6_d = v5_d*v4 + v4_d*v5

  z = c(v6,v6_d)
  return(z)
}

# select a point and a number of simulations
x  = c(1,2)
nsim = 10000

# for each simulation compute the gradient using finite differences
tic()
for(i in 1:nsim){
  a1_d = grad(myfun,x)
}
toc()

# for each simulation compute the gradient using automatic differentiation
tic()
for(i in 1:nsim){
  b1_d = myfun_d(x[1],x[2], 1, 0)
  b2_d = myfun_d(x[1],x[2], 0, 1)
}
toc()




































