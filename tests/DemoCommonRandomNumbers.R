# Demo to check the effect of Common random numbers when computing the variance of the difference of two random
# variables.

DiffUsingMC = function(AllSim, CRN){

s = rep(1,AllSim)

for (j in 1:AllSim){
  m = 7
  nsim = 2
  Y1bar = rep(1,nsim)
  Y2bar = rep(1,nsim)
  Dmean = rep(1,nsim)
  for (k in 1:nsim){

    # generate random variables using the same CRN
    if(CRN){
      a = runif(m)
      Y1 = qgamma(a, 3, 4)
      Y2 = qgamma(a, 4, 5)
    }else{
      # generate random variables using different CRN
      Y1 = qgamma(runif(m), 3, 4)
      Y2 = qgamma(runif(m), 4, 5)
    }
    D  = Y1-Y2
    Y1bar[k] = mean(Y1)
    Y2bar[k] = mean(Y2)

    Dmean[k] = Y1bar[k] - Y2bar[k]
  }

  s[j] = var(Dmean)

}
return(s)
}

AllSim = 100
Diff = DiffUsingMC(AllSim, 0) - DiffUsingMC(AllSim, 1)
hist(Diff,10)

















