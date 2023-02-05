#--------------------------------------------------------------------------------------------------
# PURPOSE: Test the prediction of the Innovation Algorithm

# 4. InnovKernel4 is mild adaptation of the innovation.kernel funciton in the itsmr package
#--------------------------------------------------------------------------------------------------


InnovKernel4 = function(Parms,gamma, x) {

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

  phi = Parms$AR
  theta = Parms$MA
  sigma2 = 1
  theta_r = c(1,theta,numeric(N))
  N = length(gamma)

  # Innovations algorithm

  p = ifelse(is.null(Parms$AR),0,length(Parms$AR))
  q = ifelse(is.null(Parms$MA),0,length(Parms$MA))
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

  # Compute xhat per equation (5.2.7)
  xhat = numeric(N)
  if (m > 1)
    for (n in 1:(m-1))
      xhat[n+1] = sum(Theta[n,1:n]*(x[n:1]-xhat[n:1]))
  for (n in m:(N-1)) {
    A = sum(Parms$AR*x[n:(n-p+1)])
    B = sum(Theta[n,1:q]*(x[n:(n-q+1)]-xhat[n:(n-q+1)]))
    xhat[n+1] = A + B
  }

  return(list(xhat = xhat, thetas=Theta[1:(n-1),],vs=v[1:(n-1)]))
}


InnovKernel4Xhat = function (Theta,x,Parms){

  N = length(x)+1
  p = ifelse(is.null(Parms$AR),0,length(Parms$AR))
  q = ifelse(is.null(Parms$MA),0,length(Parms$MA))
  m = max(p,q)

  # Compute xhat per equation (5.2.7)
  xhat = numeric(N)
  if (m > 1)
    for (n in 1:(m-1))
      xhat[n+1] = sum(Theta[n,1:n]*(x[n:1]-xhat[n:1]))
  for (n in m:(N-1)) {
    A = sum(Parms$AR*x[n:(n-p+1)])
    B = sum(Theta[n,1:q]*(x[n:(n-q+1)]-xhat[n:(n-q+1)]))
    xhat[n+1] = A + B
  }

  return(xhat)
}


theta    = 0.5
phi      = 0
sigma2   = 1
Parms    = list()
Parms$AR = phi
Parms$MA = theta
N        = 10

a        = list()
a$phi    = ifelse(is.null(Parms$AR), 0,Parms$AR)
a$theta  = ifelse(is.null(Parms$MA), 0,Parms$MA)
a$sigma2 = sigma2


gamma = aacvf(a,N-1)

set.seed(2)
x  = arima.sim(model = list( ma=theta  ), 4)

ia1 = InnovKernel4(Parms, gamma, x[1])
xhat1 = InnovKernel4Xhat(ia1$thetas, x[1], Parms)[2]

ia2 = InnovKernel4(Parms, gamma, x[1:2])
xhat2 = InnovKernel4Xhat(ia2$thetas, x[1:2], Parms)[3]

ia3 = InnovKernel4(Parms, gamma, x[1:3])
xhat3 = InnovKernel4Xhat(ia3$thetas, x[1:3], Parms)[4]

xhat1
xhat2
xhat3

















