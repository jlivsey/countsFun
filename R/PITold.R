
PIT <- function(x, count.family = c("Poisson", "mixed-Poisson","neg-Binomial", "gen-Poisson"),
                gauss.series = c("AR","FARIMA","MA"),
                estim.method = c("particlesSIS","particlesSISR"),
                H = NULL,
                p=NULL, d=NULL, q=NULL, n.mix=NULL, n=NULL, theta=NULL, phi=NULL, tht=NULL, ...)
{


  # DEFINE a cdf function from count.family input. Make sure it accepts a vector
  #        valued input.
  if(count.family=="Poisson"){
    cdf = function(x, theta){ ppois(q=x, lambda=theta[1]) }
    pdf = function(x, theta){ dpois(x, lambda=theta[1]) }
  }  else if(count.family=="mixed-Poisson"){
    if(is.null(n.mix)) stop("you must specify the number of Poissons to mix,
                            n.mix, to use count.family=mixed-Poisson")
    if(n.mix==1) stop("n.mix must be greater than 1")
    if(n.mix==2){
      cdf = function(x, theta){
        theta[1]*ppois(x, theta[2]) + (1-theta[1])*ppois(x, theta[3])
      }
      pdf = function(x, theta){
        theta[1]*dpois(x, theta[2]) + (1-theta[1])*dpois(x, theta[3])
      }
    }
    if(n.mix>2) stop("mixed-Poisson is not currently coded for n.mix>2")

  } else if(count.family=="Binomial"){
    if(is.null(n)) stop("you must specify the number of trials,
                        n, to use count.family=Binomial")
    cdf = function(x, theta){ pbinom(q = x, size = n, prob = theta) }
    pdf = function(x, theta){ dbinom(x = x, size = n, prob = theta) }
  } else if(count.family=="neg-Binomial"){
    cdf = function(x, theta){ pnbinom(q = x, size = theta[1], prob = theta[2]) }
    pdf = function(x, theta){ dnbinom(x = x, size = theta[1], prob = theta[2]) }
  } else if(count.family=="gen-Poisson"){ # VP addition
    cdf = function(x, theta1){ cdf.vec <- rep(-99,length(x))
    for (i in 1:length(x)) {
      if (x[i]>=0){
        cdf.vec[i] <- sum( exp(-(theta1[1]+(0:x[i])*theta1[2]))*theta1[1]*(theta1[1]+(0:x[i])*theta1[2])^((0:x[i])-1)/factorial((0:x[i])))
      }else{cdf.vec[i] <-0}
    }
    return(cdf.vec)}
    pdf = function(x, theta1){ exp(-(theta1[1]+x*theta1[2]))*theta1[1]*(theta1[1]+x*theta1[2])^(x-1)/factorial(x) }
  } else{ stop("please specify a valid count.family") }



  if ((gauss.series=="AR") & (estim.method=="particlesSIS")){
    if(is.null(p)) stop("you must specify the AR order, p, to use
                        gauss.series=AR")
    if(p>2) stop("the ACVF is not coded for AR models of order higher than 2
                 currently")
    if(p==1){
      #set.seed(1)
      z.rest = function(a,b){
        # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
        qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
      }
      if (is.vector(x)){
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
          zprev = z.rest(a,b)
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
            err = z.rest(a,b)
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
      }else{ # this is written just for one covariate, excluding intercept
        PDvalues = function(theta, phi, data){

          xt = data[,1]
          cvt = data[,-1]

          T1 = length(xt)
          dd = 1
          theta0 = theta[1:dd]
          theta0 = exp(rowSums(matrix(rep(theta0,each=T1),ncol=dd)*cvt))
          theta1 = cbind(theta0,rep(theta[dd+1],T1))

          N = 5000 # number of particles
          preddist = matrix(0,2,T1-1) # to collect the values of predictive distribution of interest

          a = qnorm(cdf(xt[1]-1,theta1[1,]),0,1)
          b = qnorm(cdf(xt[1],theta1[1,]),0,1)
          a = rep(a,N)
          b = rep(b,N)
          zprev = z.rest(a,b)
          zhat = phi*zprev

          wprev = rep(1,N)

          for (t in 2:T1){
            temp = rep(0,(xt[t]+1))
            for (x in 0:xt[t]){
              rt = sqrt(1-phi^2)
              a = (qnorm(cdf(x-1,theta1[t,]),0,1) - zhat)/rt
              b = (qnorm(cdf(x,theta1[t,]),0,1) - zhat)/rt
              temp[x+1] = mean(wprev*(pnorm(b,0,1) - pnorm(a,0,1)))/mean(wprev)
            }
            err = z.rest(a,b)
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
      }
    } else if(p==2){
      #set.seed(1)
      z.rest = function(a,b){
        # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
        qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
      }
      likSIS = function(theta, data){
        # set.seed(1)
        theta1 = theta[theta1.idx]
        n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
        theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 2)
        phi = theta[theta2.idx]

        if ( (phi[2]+phi[1]<1) & (phi[2]-phi[1]<1) & (abs(phi[2])<1)){
          xt = data
          T1 = length(xt)
          N = 2000 # number of particles
          prt = matrix(0,N,T1) # to collect all particles
          wgh = matrix(0,N,T1) # to collect all particle weights

          a1 = qnorm(cdf(xt[1]-1,theta1),0,1)
          b1 = qnorm(cdf(xt[1],theta1),0,1)
          a1 = rep(a1,N)
          b1 = rep(b1,N)
          zprev1 = z.rest(a1,b1)
          zhat1 = phi[1]*zprev1
          prt[,1] = zhat1

          wprev = rep(1,N)
          wgh[,1] = wprev

          phi0 = if(abs(phi[1])<0.9){phi[1]}else{0.9*sign(phi[1])}
          rt2 = sqrt(1-phi0^2)
          a2 = (qnorm(cdf(xt[2]-1,theta1),0,1) - phi[1]*zprev1)/rt2
          b2 = (qnorm(cdf(xt[2],theta1),0,1) - phi[1]*zprev1)/rt2
          err2 = z.rest(a2,b2)
          znew2 = phi0*zprev1 + rt2*err2
          znew = c(znew2,zprev1)
          zhat2 = sum(phi*znew)
          prt[,2] = zhat2

          wgh[,2] = wprev*(pnorm(b2,0,1) - pnorm(a2,0,1))
          wprev = wgh[,2]

          zprev = znew

          for (t in 3:T1)
          {
            rt = 1/sqrt( ((1-phi[2])/(1+phi[2])) / ((1-phi[2])^2-phi[1]^2) )
            a = (qnorm(cdf(xt[t]-1,theta1),0,1) - sum(phi*zprev))/rt
            b = (qnorm(cdf(xt[t],theta1),0,1) - sum(phi*zprev))/rt
            err = z.rest(a,b)
            znew = sum(phi*zprev) + rt*err
            znew = c(znew,zprev[1])
            zhat = sum(phi*znew)
            prt[,t] = zhat
            zprev = znew

            wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
            wprev = wgh[,t]
          }

          lik = pdf(xt[1],theta1)*mean(wgh[,T1])
          nloglik = (-2)*log(lik)

          out = if (is.na(nloglik)) Inf else nloglik
          return(out)
        }else{
          out = Inf
          return(out)
        }
      }
    } else{ stop("the p specified is not valid") }
  }

  if ((gauss.series=="MA") & (estim.method=="particlesSIS")){ # VP change
    if(q==1){
      #set.seed(1)
      z.rest = function(a,b){
        # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
        qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
      }
      PDvalues = function(theta, tht, data){
        #set.seed(1)
        xt = data
        T1 = length(xt)
        N = 1000 # number of particles

        a = qnorm(cdf(xt[1]-1,theta1),0,1)
        b = qnorm(cdf(xt[1],theta1),0,1)
        a = rep(a,N)
        b = rep(b,N)
        zprev = z.rest(a,b)
        rt0 = 1+tht^2
        zhat = tht*zprev/rt0
        prt[,1] = zhat

        wprev = rep(1,N)
        wgh[,1] = wprev

        for (t in 2:T1)
        {
          rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
          rt = sqrt(rt0/(1+tht^2))
          a = (qnorm(cdf(xt[t]-1,theta1),0,1) - zhat)/rt
          b = (qnorm(cdf(xt[t],theta1),0,1) - zhat)/rt
          err = z.rest(a,b)
          znew = zhat + rt*err
          zhat = tht*(znew-zhat)/rt0
          prt[,t] = zhat

          wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
          wprev = wgh[,t]
        }

        lik = pdf(xt[1],theta1)*mean(wgh[,T1])
        nloglik = (-2)*log(lik)

        out = if (is.na(nloglik)) Inf else nloglik
        return(out)
      }
    }else if(q==2){ # VP addition
      z.rest = function(a,b){
        # Generates N(0,1) variables restricted to (ai,bi),i=1,...n
        qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
      }
      if (is.vector(x)){
        PDvalues = function(theta, tht, data){
            #set.seed(1)

            xt = data
            T1 = length(xt)
            N = 1000 # number of particles

            preddist = matrix(0,2,T1-1) # to collect the values of predictive distribution of interest

            gam1 = tht[1]*(1+tht[2])/(1+tht[1]^2+tht[2]^2)
            gam2 = tht[2]/(1+tht[1]^2+tht[2]^2)


            a = qnorm(cdf(xt[1]-1,theta),0,1)
            b = qnorm(cdf(xt[1],theta),0,1)
            a = rep(a,N)
            b = rep(b,N)
            zprev1 = z.rest(a,b)
            tht11 = gam1
            zhat1 = tht11*zprev1

            rt1 = ( 1 - tht11^2 )^(1/2)
            wprev = rep(1,N)

            temp = rep(0,(xt[2]+1))
            for (x in 0:xt[2]){
              a2 = (qnorm(cdf(x-1,theta),0,1) - zhat1)/rt1
              b2 = (qnorm(cdf(x,theta),0,1) - zhat1)/rt1
              temp[x+1] = mean(wprev*(pnorm(b2,0,1) - pnorm(a2,0,1)))/mean(wprev)
            }
            err2 = z.rest(a2,b2)
            znew2 = zhat1 + rt1*err2
            znew_all = cbind(znew2,zprev1)

            tht22 = gam2
            tht21 = (gam1 - tht11*tht22)/rt1^2
            tht_all = cbind(rep(tht21,N),rep(tht22,N))
            zhat_all = cbind(zhat1,rep(0,N))
            zhat2 = rowSums(tht_all*(znew_all-zhat_all))

            rt2 = (1 - tht21^2*rt1^2 - tht22^2)^(1/2)

            wprev = wprev*(pnorm(b2,0,1) - pnorm(a2,0,1))

            rt_all = c(rt2,rt1)
            zhat_all = cbind(zhat2,zhat1)

            if (xt[2]==0){
              preddist[,2-1] = c(0,temp[1])
            }else{
              preddist[,2-1] = cumsum(temp)[xt[2]:(xt[2]+1)]
            }

            for (t in 3:T1)
            {

              temp = rep(0,(xt[t]+1))
              for (x in 0:xt[t]){
                a = (qnorm(cdf(x-1,theta),0,1) - zhat_all[,1])/rt_all[1]
                b = (qnorm(cdf(x,theta),0,1) - zhat_all[,1])/rt_all[1]
                temp[x+1] = mean(wprev*(pnorm(b,0,1) - pnorm(a,0,1)))/mean(wprev)
              }
              err = z.rest(a,b)
              znew = zhat_all[,1] + rt_all[1]*err
              znew_all = cbind(znew,znew_all[,1])

              thtt2 = gam2/rt_all[2]^2
              thtt1 = (gam1 - tht_all[1]*thtt2*rt_all[2]^2)/rt_all[1]^2
              tht_all = cbind(rep(thtt1,N),rep(thtt2,N))
              zhat2 = rowSums(tht_all*(znew_all-zhat_all))

              rt = (1 - thtt1^2*rt_all[1]^2 - thtt2^2*rt_all[2]^2)^(1/2)

              wprev =  wprev*(pnorm(b,0,1) - pnorm(a,0,1))

              rt_all = c(rt,rt_all[1])
              zhat_all = cbind(zhat2,zhat_all[,1])

              if (xt[t]==0){
                preddist[,t-1] = c(0,temp[1])
              }else{
                preddist[,t-1] = cumsum(temp)[xt[t]:(xt[t]+1)]
              }

            }

            return(preddist)

        }
      }else{
        PDvalues = function(theta, tht, data){
          #set.seed(1)

          xt = data[,1]
          cvt = data[,-1]
          cvt = cbind(rep(1,length(xt)),cvt)
          T1 = length(xt)
          dd = dim(cvt)[2]
          theta0 = theta[1:dd]
          theta0 = exp(rowSums(matrix(rep(theta0,each=T1),ncol=dd)*cvt))
          theta1 = cbind(theta0,rep(theta[dd+1],T1))

          N = 1000 # number of particles
          preddist = matrix(0,2,T1-1) # to collect the values of predictive distribution of interest

          gam1 = tht[1]*(1+tht[2])/(1+tht[1]^2+tht[2]^2)
          gam2 = tht[2]/(1+tht[1]^2+tht[2]^2)

          a = qnorm(cdf(xt[1]-1,theta1[1,]),0,1)
          b = qnorm(cdf(xt[1],theta1[1,]),0,1)
          a = rep(a,N)
          b = rep(b,N)
          zprev1 = z.rest(a,b)
          tht11 = gam1
          zhat1 = tht11*zprev1

          rt1 = ( 1 - tht11^2 )^(1/2)
          wprev = rep(1,N)

          temp = rep(0,(xt[2]+1))
          for (x in 0:xt[2]){
            a2 = (qnorm(cdf(x-1,theta1[2,]),0,1) - zhat1)/rt1
            b2 = (qnorm(cdf(x,theta1[2,]),0,1) - zhat1)/rt1
            temp[x+1] = mean(wprev*(pnorm(b2,0,1) - pnorm(a2,0,1)))/mean(wprev)
          }
          err2 = z.rest(a2,b2)
          znew2 = zhat1 + rt1*err2
          znew_all = cbind(znew2,zprev1)

          tht22 = gam2
          tht21 = (gam1 - tht11*tht22)/rt1^2
          tht_all = cbind(rep(tht21,N),rep(tht22,N))
          zhat_all = cbind(zhat1,rep(0,N))
          zhat2 = rowSums(tht_all*(znew_all-zhat_all))

          rt2 = (1 - tht21^2*rt1^2 - tht22^2)^(1/2)

          wprev = wprev*(pnorm(b2,0,1) - pnorm(a2,0,1))

          rt_all = c(rt2,rt1)
          zhat_all = cbind(zhat2,zhat1)

          if (xt[2]==0){
            preddist[,2-1] = c(0,temp[1])
          }else{
            preddist[,2-1] = cumsum(temp)[xt[2]:(xt[2]+1)]
          }

          for (t in 3:T1)
          {

            temp = rep(0,(xt[t]+1))
            for (x in 0:xt[t]){
              a = (qnorm(cdf(x-1,theta1[t,]),0,1) - zhat_all[,1])/rt_all[1]
              b = (qnorm(cdf(x,theta1[t,]),0,1) - zhat_all[,1])/rt_all[1]
              temp[x+1] = mean(wprev*(pnorm(b,0,1) - pnorm(a,0,1)))/mean(wprev)
            }
            err = z.rest(a,b)
            znew = zhat_all[,1] + rt_all[1]*err
            znew_all = cbind(znew,znew_all[,1])

            thtt2 = gam2/rt_all[2]^2
            thtt1 = (gam1 - tht_all[1]*thtt2*rt_all[2]^2)/rt_all[1]^2
            tht_all = cbind(rep(thtt1,N),rep(thtt2,N))
            zhat2 = rowSums(tht_all*(znew_all-zhat_all))

            rt = (1 - thtt1^2*rt_all[1]^2 - thtt2^2*rt_all[2]^2)^(1/2)

            wprev = wprev*(pnorm(b,0,1) - pnorm(a,0,1))

            rt_all = c(rt,rt_all[1])
            zhat_all = cbind(zhat2,zhat_all[,1])

            if (xt[t]==0){
              preddist[,t-1] = c(0,temp[1])
            }else{
              preddist[,t-1] = cumsum(temp)[xt[t]:(xt[t]+1)]
            }

          }

          return(preddist)

        }
      }



    }
  }




  PITvalues = rep(0,H)

  if (gauss.series == "MA"){
    predd = PDvalues(theta, tht, x)
  }

  if (gauss.series == "AR"){
    predd = PDvalues(theta, phi, x)
  }



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
