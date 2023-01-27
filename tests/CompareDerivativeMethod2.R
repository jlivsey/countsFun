# PURPOSE: Write a simple implementation of the poisson cdf that also computes
#          its derivative with respect to lambda. Ideally this would be done using
#          the incomplete gamma function. I tried to do that but autodiff seems to have some
#          issues (or I had bugs I couldnt find). Even in this simple implementation with
#          the for loop the cdf evaluation is only a little bit slower than the standard R function
#          However, the autodiff implementation of the derivative significantly outperforms the
#          finite differences.
#----------------------------------------------------------------------------------------------------#


ppois_d = function(x, a, x_d, a_d){
  # compute the cdf of the poisson and its derivative with respct to a (lamda)
  # the code is only correct for x_d=0 and a_d = 1
  if ( x < 0 ){
    cdf = 0.0
    cdf_d = 0.0
  }else{

    new = exp(-a)
    new_d = -exp(-a)*a_d

    sum2 = new;
    sum2_d = new_d

    for (i in 1 : x){
      last = new
      last_d = new_d

      new = last * a / i
      new_d = (last_d*a + last*a_d)/i

      sum2 = sum2 + new
      sum2_d = sum2_d + new_d

    }

    cdf = sum2
    cdf_d = sum2_d

  }

  return(c(cdf,cdf_d))
}

# fix a point x and take several values of a
x = c(2,4)
nsim = 10000


# Evaluate derivative using numerical differentiation
t1 = tic()
for(i in 1:nsim){
  x1 = ppois(x[1], x[2])
  x_d = grad(function(x)ppois(x[1],x[2]),x)[2]
}
t1 = tic() - t1

# Evaluate derivative using autodiff
t2 = tic()
for(i in 1:nsim){
  y   = ppois2(x[1], x[2],0,1)
  y   = y[1]
  y_d = y[2]
}
t2 = tic() - t2

t1
t2
