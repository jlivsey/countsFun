#--------------------------------------------------------------------------------------------------
# PURPOSE: Compare different implementations of Innov Algorithms that I found online
#
# 1. innovations.algorithm found at http://faculty.washington.edu/dbp/s519/R-code/innovations-algorithm.R
# 2. InnovAlgMA same as a bove but with while loop implemented by me
# 3. InnovKernel this is the implementaitons in itsmr package with the output v = v/v[1]
# 4. InnovKernel2, InnovKernel3, adaptations of InnovKernel to run faster
#--------------------------------------------------------------------------------------------------

# innovations algorithm code
innovations.algorithm <- function(acvf,n.max=length(acvf)-1){
  # I found this implementation of IA online, I need to check it
  # http://faculty.washington.edu/dbp/s519/R-code/innovations-algorithm.R
  # if(acvf[1]!=1) acvf = acvf/acvf[1]
  thetas <- vector(mode="list",length=n.max)
  vs <- rep(acvf[1],n.max+1)
  for(n in 1:n.max){
    thetas[[n]] <- rep(0,n)
    thetas[[n]][n] <- acvf[n+1]/vs[1]
    if(n>1){
      for(k in 1:(n-1)){
        js <- 0:(k-1)
        thetas[[n]][n-k] <- (acvf[n-k+1] - sum(thetas[[k]][k-js]*thetas[[n]][n-js]*vs[js+1]))/vs[k+1]
      }
    }
    js <- 0:(n-1)
    vs[n+1] <- vs[n+1] - sum(thetas[[n]][n-js]^2*vs[js+1])
  }
  return(structure(list(vs=vs,thetas=thetas)))
}

# innovations algorithm code for MA series
InnovAlgMA = function(acvf,n.max=length(acvf)-1,q,maxdiff=0.000000001){
  # In the case of MA series the innovation alogrithm coefficients
  # converege fast to the MA coefficients
  if(acvf[1]!=1) acvf = acvf/acvf[1]
  thetas   = list()
  vs   = rep(acvf[1],n.max+1)
  Diff = rep(1,q)
  n=1
  while( abs(max(Diff))>maxdiff && n<n.max ){

    thetas[[n]] <- rep(0,q)
    thetas[[n]][n] = ifelse(n<=q, acvf[n+1]/vs[1],0)
    if(n>1){
      for(k in 1:(n-1)){
        js <- 0:(k-1)
        thetas[[n]][n-k] <- (acvf[n-k+1] - sum(thetas[[k]][k-js]*thetas[[n]][n-js]*vs[js+1]))/vs[k+1]
        Diff[k] = thetas[[n]][n-k] - thetas[[n-1]][n-k]
      }
      #Diff[n] = thetas[[n]][1:min(n,q)] - thetas[[n-1]][1:min(n,q)]
    }
    js <- 0:(n-1)
    vs[n+1] <- vs[n+1] - sum(thetas[[n]][n-js]^2*vs[js+1])
    n = n+1
  }

  return(structure(list(vs=vs[1:(n-1)],thetas=lapply(thetas[ ][1:(n-1)], function(x) {x[1:q]}))))
}



InnovKernel = function(a, N) {

  # Compute autocovariance kappa(i,j) per equation (3.3.3)

  # Optimized for i >= j and j > 0

  kappa = function(i,j) {
    if (j > m)
      return(sum(theta_r[1:(q+1)] * theta_r[(i-j+1):(i-j+q+1)]))
    else if (i > 2*m)
      return(0)
    else if (i > m)
      return((gamma[i-j+1] - sum(phi * gamma[abs(seq(1-i+j,p-i+j))+1]))/sigma2)
    else
      return(gamma[i-j+1]/sigma2)
  }

  phi = a$phi
  theta = a$theta
  sigma2 = a$sigma2


  theta_r = c(1,theta,numeric(N))

  # Autocovariance of the model

  gamma = aacvf(a,N-1)
  #gamma = ARMAacf(ma = theta,lag.max = N)
  # Innovations algorithm

  p = length(phi)
  q = length(theta)
  m = max(p,q)

  Theta = matrix(0,N-1,N-1)
  v = numeric(N)
  v[1] = kappa(1,1)
  for (n in 1:(N-1)) {
    for (k in 0:(n-1)) {
      u = kappa(n+1,k+1)
      s = 0
      if (k > 0) s = sum(Theta[k,k:1] * Theta[n,n:(n-k+1)] * v[1:k])
      Theta[n,n-k] = (u - s) / v[k+1]
    }
    s = sum(Theta[n,n:1]^2 * v[1:n])
    v[n+1] = kappa(n+1,n+1) - s
  }
  v = v/v[1]

  # Compute xhat per equation (5.2.7)
  # xhat = numeric(N)
  # if (m > 1)
  #   for (n in 1:(m-1))
  #     xhat[n+1] = sum(Theta[n,1:n]*(x[n:1]-xhat[n:1]))
  # for (n in m:(N-1)) {
  #   A = sum(phi*x[n:(n-p+1)])
  #   B = sum(Theta[n,1:q]*(x[n:(n-q+1)]-xhat[n:(n-q+1)]))
  #   xhat[n+1] = A + B
  # }

  return(list(theta=Theta,v=v))
}



InnovKernel2 = function(a, N) {

  # Compute autocovariance kappa(i,j) per equation (3.3.3)

  # Optimized for i >= j and j > 0

  kappa = function(i,j) {
    if (j > m)
      return(sum(theta_r[1:(q+1)] * theta_r[(i-j+1):(i-j+q+1)]))
    else if (i > 2*m)
      return(0)
    else if (i > m)
      return((gamma[i-j+1] - sum(phi * gamma[abs(seq(1-i+j,p-i+j))+1]))/sigma2)
    else
      return(gamma[i-j+1]/sigma2)
  }

  phi = a$phi
  theta = a$theta
  sigma2 = a$sigma2


  theta_r = c(1,theta,numeric(N))

  # Autocovariance of the model

  gamma = aacvf(a,N-1)
  #gamma = ARMAacf(ma = theta,lag.max = N)
  # Innovations algorithm

  p = length(phi)
  q = length(theta)
  m = max(p,q)

  Theta = matrix(0,N-1,N-1)
  v = numeric(N)
  Diff  = 1
  v[1] = kappa(1,1)
  # n=1 k=0 case
  Theta[1,1] = kappa(2,1)/v[1]
  v[2] = kappa(2,2) - Theta[1,1]^2 * v[1]

  n=2
  maxdiff=0.000000001
  while(abs(Diff)>maxdiff && n<N ) {
    for (k in 0:(n-1)) {
      u = kappa(n+1,k+1)
      s = 0
      if (k > 0) s = sum(Theta[k,k:1] * Theta[n,n:(n-k+1)] * v[1:k])
      Theta[n,n-k] = (u - s) / v[k+1]
    }
    s = sum(Theta[n,n:1]^2 * v[1:n])
    v[n+1] = kappa(n+1,n+1) - s
    Diff = v[n+1]-v[n]
    n = n+1
  }
  v = v/v[1]

  return(list(theta=Theta[1:(n-1),],v=v[1:(n-1)]))
}



InnovKernel3 = function(gamma,phi,theta,sigma2, N) {

  # Compute autocovariance kappa(i,j) per equation (3.3.3)

  # Optimized for i >= j and j > 0

  kappa = function(i,j) {
    if (j > m)
      return(sum(theta_r[1:(q+1)] * theta_r[(i-j+1):(i-j+q+1)]))
    else if (i > 2*m)
      return(0)
    else if (i > m)
      return((gamma[i-j+1] - sum(phi * gamma[abs(seq(1-i+j,p-i+j))+1]))/sigma2)
    else
      return(gamma[i-j+1]/sigma2)
  }

  # phi = a$phi
  # theta = a$theta
  # sigma2 = a$sigma2


  theta_r = c(1,theta,numeric(N))

  # Autocovariance of the model

  # gamma = aacvf(a,N-1)
  #gamma = ARMAacf(ma = theta,lag.max = N)
  # Innovations algorithm

  p = length(phi)
  q = length(theta)
  m = max(p,q)

  # if(acvf[1]!=1) acvf = acvf/acvf[1]
  thetas   = list()
  vs   = rep(NA,N+1)
  vs[1] = kappa(1,1)
  Diff = rep(1,q)
  n=1
  maxdiff=0.000000001
  while( abs(max(Diff))>maxdiff && n<N ){

    thetas[[n]] <- rep(0,q)
    thetas[[n]][n] = ifelse(n<=q, kappa(2,1)/vs[1],0)
    if(n>1){
      for(k in 1:(n-1)){
        js <- 0:(k-1)
        thetas[[n]][n-k] <- (kappa(n+1,k+1) - sum(thetas[[k]][k-js]*thetas[[n]][n-js]*vs[js+1]))/vs[k+1]
        Diff[k] = thetas[[n]][n-k] - thetas[[n-1]][n-k]
      }
      #Diff[n] = thetas[[n]][1:min(n,q)] - thetas[[n-1]][1:min(n,q)]
    }
    js <- 0:(n-1)
    vs[n+1] <- kappa(n+1,n+1) - sum(thetas[[n]][n-js]^2*vs[js+1])
    vs
    n = n+1
  }
  vs = vs/vs[1]
  return(structure(list(vs=vs[1:(n-1)],thetas=lapply(thetas[ ][1:(n-1)], function(x) {x[1:q]})) ))
}



theta    = NULL
phi      = c(0.8, -0.25)
sigma2   = 1
N        = 5
a        = list()
a$phi    = phi
a$theta  = theta
a$sigma2 = sigma2
acvf2    = aacvf(a,N-1)
acvf     = ARMAacf(ar = phi,lag.max = N-1)
names(acvf) = NULL
ia0      = innovations.algorithm(acvf)
ia4      = InnovKernel(a, N)


# ia2      = innovations.algorithm(acvf2)
# ia3      = InnovKernel(a,N)
# #ia4      = InnovAlgMA(acvf2,n.max=length(acvf)-1,q,maxdiff=0.000000001)
# ia5      = InnovKernel2(a,N)


nsim = 1
t1 = tic()
for (i in 1:nsim){
  ia1      = InnovAlgMA(acvf2,n.max=length(acvf)-1,q=length(theta),maxdiff=0.000000001)
}
t1 = tic()-t1

t2 = tic()
for (i in 1:nsim){
  ia2      = InnovKernel2(a,N)
}
t2 = tic()-t2

t3 = tic()
for (i in 1:nsim){
  ia3      = InnovKernel3(acvf2,phi,theta,sigma2, N)
}
t3 = tic()-t3

t4 = tic()
for (i in 1:nsim){
  ia4      = InnovKernel(a, N)
}
t4 = tic()-t4
ia1
ia2
ia3
ia4








