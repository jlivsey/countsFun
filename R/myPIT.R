
myPIT <- function(x, H, theta,phi){
  cdf = function(x, theta){ ppois(x, lambda=theta[1]) }
  pdf = function(x, theta){ dpois(x, lambda=theta[1]) }



  PDvalues = function(theta, phi, data){
    #set.seed(1)
    xt = data
    T1 = length(xt)
    N = 10000 # number of particles
    preddist = matrix(0,2,T1-1) # to collect the values of predictive distribution of interest

    a = qnorm(cdf(xt[1]-1,theta),0,1)
    b = qnorm(cdf(xt[1],theta),0,1)
    a = rep(a,N)
    b = rep(b,N)
    zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
    zhat = phi*zprev

    wprev = rep(1,N)

    for (t in 2:T1){
      temp = rep(0,(xt[t]+1))
      for (x in 0:xt[t]){
        rt = sqrt(1-phi^2)
        a = (qnorm(cdf(x-1,theta),0,1) - zhat)/rt
        b = (qnorm(cdf(x,theta),0,1) - zhat)/rt
        temp[x+1] = mean(wprev*(pnorm(b,0,1) - pnorm(a,0,1)))/mean(wprev)
      }
      err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
      znew = phi*zprev + rt*err
      zhat = phi*znew

      wprev = wprev*(pnorm(b,0,1) - pnorm(a,0,1))

      if (xt[t]==0){
        preddist[,t-1] = c(0,temp[1])
      }else{
        preddist[,t-1] = cumsum(temp)[xt[t]:(xt[t]+1)]
      }
    }
    return(preddist)
  }



  PITvalues = rep(0,H)
  predd = PDvalues(theta, phi, x)




  predd1 = predd[1,]
  predd2 = predd[2,]
  Tr = length(predd1)

  for (h in 1:H){
    id1 = (predd1 < h/H)*(h/H < predd2)
    id2 = (h/H >= predd2)
    tmp1 = (h/H-predd1)/(predd2-predd1)
    tmp1[!id1] = 0
    tmp2 = rep(0,Tr)
    tmp2[id2] = 1
    PITvalues[h] = mean(tmp1+tmp2)
  }
  PITvalues = c(0,PITvalues)

  return(diff(PITvalues))

}
