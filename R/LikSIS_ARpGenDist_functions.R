# Add some functions that I will need in particle filtering approximation of 
# likelihood. See file LikSIS_ARpGenDist.R

#####------------------------------------------------------###################
# Add some functions that I will need
get_rand_state <- function() {
  # Using `get0()` here to have `NULL` output in case object doesn't exist.
  # Also using `inherits = FALSE` to get value exactly from global environment
  # and not from one of its parent.
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

set_rand_state <- function(state) {
  # Assigning `NULL` state might lead to unwanted consequences
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
  }
}

# Generate AR series
sim_pois_ar = function(n, phi, lam){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qpois(pnorm(z), lam)
  return(x)
}

sim_negbin_ar = function(n, phi, r, prob){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qnbinom(pnorm(z), r,prob)
  return(x)
}

sim_MixPois_ar <- function(n,phi,lam1,lam2, prob){
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x  = qMixedPoisson(pnorm(z),lam1,lam2, prob)
  return(x)
}


# mixed poisson cdf
pMixedPoisson = function(q, lam1, lam2, prob){
  # theta[1] is the mixing probability
  # theta[2], theta[3] are the lambda parameters
  prob*ppois(q, lam1) + (1-prob)*ppois(q, lam2)
}

dMixedPoisson = function(x, lam1, lam2, prob){
  prob*dpois(x, lam1) + (1-prob)*dpois(x, lam2)
}

qMixedPoisson = function(y, lam1, lam2, prob){
  yl = length(y)
  x = rep(0,yl)
  for (n in 1:yl){
    while(pMixedPoisson(x[n], lam1, lam2, prob) < y[n]){ # R qpois would use <y; this choice makes the function right-continuous; this does not really matter for our model
      x[n] = x[n]+1
    } 
  }
  return(x)
}


# Generalized Poisson cdf, pdf
pGenPoisson = function(q, theta, lam){ cdf.vec <- rep(-99,length(q))
for (i in 1:length(q)) {
  if (q[i]>=0){
    cdf.vec[i] <- sum( exp(-(theta+(0:q[i])*lam))*theta*(theta+(0:q[i])*lam)^((0:q[i])-1)/factorial((0:q[i])))
  }else{cdf.vec[i] <-0}
}
return(cdf.vec)}

dGenPoisson = function(x, theta, lam){ exp(-(theta+x*lam))*theta*(theta+x*lam)^(x-1)/factorial(x) }



FitMultiplePF = function(initialParam, data, CountDist, nfit, ParticleSchemes){
  # Let nparts by the length of the vector ParticleSchemes.
  # This function fits maximizes the PF likelihood, nfit manys times for nparts many choices
  # of particle numbers, thus yielding a total of nfit*nparts many estimates.

  # how many choices for the number of particles
  nparts = length(ParticleSchemes)
  nparms = length(initialParam)
  
  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst = matrix(0,nrow=nfit*nparts,ncol=nparms)
  se =  matrix(NA,nrow=nfit*nparts,ncol=nparms)
  loglik = rep(0,nfit*nparts)
  

  # Each realization will be fitted nfit*nparts many times
  for (j in 1:nfit){
    set.seed(j)
    # for each fit repeat for different number of particles
    for (k in 1:nparts){
      # number of particles to be used
      ParticleNumber = ParticleSchemes[k]
      
      # remove the ParticleNumber from the likelihood function arguments
      myfun = function(theta,data)LikSISGenDist_ARp_Res(theta,data,ParticleNumber, CountDist)
      
      # run optimization for our model
      optim.output <- optim(par = initialParam, fn = myfun,
                            data=data,
                            hessian=TRUE, method = "BFGS")
      
      # save estimates, loglik value and diagonal hessian
      ParmEst[nfit*(k-1)+j,]  = optim.output$par
      loglik[nfit*(k-1) +j]   = optim.output$value
      se[nfit*(k-1)+j,]       = sqrt(abs(diag(solve(optim.output$hessian))))
    }
  }
  
  All = cbind(ParmEst, se, loglik)
  return(All)
} 



#####------------------------------------------------------###################