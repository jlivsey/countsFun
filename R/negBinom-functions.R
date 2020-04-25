
#---------Covariance matrix---------#
CovarNegBinAR = function(n,r, p, phi, N){
  #====================================================================================#
  # PURPOSE    Compute the covariance matrix of a NegBin AR series.
  #
  # INPUT
  #   r,p      Marginal parameters
  #   phi      AR parameter
  #   n        size of the matrix
  #   N        truncation of relation (21)
  #
  # Output
  #   GAMMA    covariance matrix ofcount series
  #
  # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date       April 2020
  # Version    3.6.3
  #====================================================================================#

  # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
  HC = HermCoefNegBin(r,p)

  # ARMA autocorrelation function
  ar.acf <- ARMAacf(ar = phi, lag.max = n)

  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  gamma_x = CountACVF(h = 0:(n-1), myacf = ar.acf, g = HC)

  # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = toeplitz(gamma_x)
  return(GAMMA)
}

#----------Link coefficients with Regressor---------#
LinkCoef_Reg = function(r, p, M){
  #====================================================================================#
  # PURPOSE    Compute the product  factorial(k)*g_{t1,k}*g_{t2,k} in relation (67) in
  #            https://arxiv.org/pdf/1811.00203.pdf when the regressor variable is only
  #            one dummy variable.
  #
  # INPUT
  #   r, p     NB Marginal parameters
  #
  # Output
  #   l        An MX3 matrix of link coefficients. M is the default truncatuon of relation
  #            (67) and 3 is the number of different combinations of products
  #            g_{t1,k}*g_{t2,k} (since for each t I will either have 1 or 0 for the
  #            regressor variable).
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       April 2020
  # Version    3.6.3
  #====================================================================================#


  # compute Hermite coefficients from relation (21)
  g1 = HermCoefNegBin_2(r, unique(p)[1], M)
  g2 = HermCoefNegBin_2(r, unique(p)[2], M)
  HC = cbind(g1, g2)

  # Compute the products g1^2, g1*g2 and g2^2 of the HCs
  HCprod = matrix(NA, M, 3)
  for(i in 0:(length(unique(p))-1) ){
    for(j in i: (length(unique(p))-1) ){
      HCprod[,i+j+1] = HC[,i+1]*HC[,j+1]
    }
  }

  return(HCprod)
}


#---------Covariance matrix with regressor---------#
CovarNegBinAR_Reg = function(n, r, p, phi, M){
  #====================================================================================#
  # PURPOSE    Compute the covariance matrix of a NegBin AR series that
  #            includes one dummy variable as a regressor. Here p depends on each
  #            observation. See relation (67) in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   r, p     NB MArginal parameters
  #   phi      AR parameter
  #   n        size of the matrix
  #   M        truncation of relation (67)
  #
  # Notes      Since X is a dummy there are only two values of m=exp(X%*%beta),
  #            and hence only two values of the NeG Bin probability of p, and hence
  #            only two differenet Hermite Coeffcients I need to compute (for each k).
  #
  # Output
  #   GAMMA    covariance matrix ofcount series with a dummy regressor
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       March 2020
  # Version    3.6.3
  #====================================================================================#

  # Compute ARMA autocorrelation function
  All.ar = apply(as.matrix(ARMAacf(ar = phi, lag.max = n)), 1,function(x)x^(1:M))

  # Compute the link coefficients l_k = factorial(k)*g_{t1,k}*g_{t2,k}
  linkCoef = LinkCoef_Reg(r, p, M)

  # keep track of which indices each unique HC is located at
  index = cbind(unique(p)[1]==p, unique(p)[2]==p)

  # STEP 3: Implement relation 67
  k = 1:M
  kfac = factorial(k)
  G = matrix(NA,n,n)
  for(t1 in 0:(n-1)){
    for(t2 in 0:t1 ){
      h = abs(t1-t2)+1
      G[t1+1,t2+1]= sum(kfac*linkCoef[,index[t1+1,2]+index[t2+1,2]+1]*All.ar[,h])
    }
  }
  G = symmetrize(G, update.upper=TRUE)
  return(G)
}


#---------Hermitte Coefficients for all k---------#
HermCoefNegBin <- function(r, p, N, nHC){
  #====================================================================================#
  # PURPOSE    Compute all Hermite Coefficients. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   r,p      Marginal parameters
  #   maxCoef  number of coefficients to return. Default = 20
  #   N        truncation of relation (21)
  #   nHC      number of Hermitte coefficients
  #
  # Output
  #   HC       All Hermite coeficients
  #
  # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date       April 2020
  # Version    3.6.3
  #====================================================================================#

  h = 1:nHC #check me
  HC = rep(NA, length(h)) # storage
  for(i in h) {
    HC[i] <- HermCoefNegBin_k(r, p , k = i, N)
  }

  return(HC)

}


#---------Hermitte Coefficients for all k---------#
HermCoefNegBin_2 <- function(r, p, nHC){
  #====================================================================================#
  # PURPOSE    Compute all Hermite Coefficients. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   r,p      Marginal parameters
  #   nHC      number of coefficients to return
  #
  # Output
  #   HC       All Hermite coeficients
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #====================================================================================#

  # Compute truncation of relation (21) in arxiv link above
  N = sapply(unique(p),function(x)which(pnbinom(1:1000, r,1-x)>=1-1e-17)[1])-1
  N[is.na(N)] = 1000

  h = 1:nHC # 20 is truncation of (67)
  HC = rep(NA, length(h)) # storage
  for(i in h) {
    HC[i] <- HermCoefNegBin_k(r, p , k = i, N)
  }

  return(HC)

}


#---------Hermitte Coefficients for on k---------#
HermCoefNegBin_k <- function(r, p, k, N){
  #====================================================================================#
  # PURPOSE    Compute kth Hermite Coefficient. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   r,p      Marginal parameters
  #   k        index of Hermite coefficient
  #
  # Output
  #   HC_k     kth hermite coeficient
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #====================================================================================#

  # function for kth Hermite Polynomial
  her <- function(x){
    evalHermPolynomial(k-1, x)
  }

  # # truncation numbe: check me
  # N <- which(round(pnbinom(1:1000, r,p), 7) == 1)[1]
  # if(is.na(N)){
  #   N=1000
  # }
  # N = max(sapply(unique(p),function(x)which(round(pnbinom(1:1000, r,x), 7) == 1)[1]))

  # compute terms in the sum of relation (21) in
  terms <- exp((-qnorm(pnbinom(0:max(N), r,1-p, lower.tail= TRUE))^2)/2) *
    her(qnorm(pnbinom(0:max(N), r,1-p, lower.tail = TRUE)))

  # take the sum of all terms
  HC_k <- sum(terms) / (sqrt(2*pi) *  factorial(k))
  return(HC_k)
}


#---------simulate negbin series (r,p) parametrization---------#
sim_negbin_ar = function(n, phi, r,p){
  #====================================================================================#
  # PURPOSE   Simulate Poisson series with AR structure. See relation (1)
  #           in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   n        series length
  #   phi      AR parameter
  #   r,p      Marginal Parameters
  #
  # Output
  #   x        Poisson series
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #====================================================================================#
  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  x = qnbinom(pnorm(z), r,1-p)
  return(x)
}


#---------simulate negbin series (r,m) parametrization---------#
sim_negbin_ar_2 = function(n, phi, r, m){
  #====================================================================================#
  # PURPOSE   Simulate Neg Bin series with AR structure. See relation (1)
  #           in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   n        series length
  #   phi      AR parameter
  #   r,m      Marginal Parameters (mean parametrization)
  #
  # Output
  #   x        Neg Bin series
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       March 2020
  # Version    3.6.1
  #====================================================================================#

  z = arima.sim(model = list(ar=phi), n = n); z = z/sd(z) # standardized
  y = qnbinom(pnorm(z), size = r,mu = m)
  return(y)
}


#---------Covariance matrix version 2---------#
CovarNegBinAR_2 = function(n,r, m, phi){
  #====================================================================================#
  # PURPOSE    Compute the covariance matrix of a NegBin AR series.
  #            This is the parametrization through the mean used in GLM
  #
  # INPUT
  #   r,m      Marginal parameters (m is the mean)
  #   phi      AR parameter
  #   n        size of the matrix
  #
  # Output
  #   GAMMA    covariance matrix ofcount series
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #====================================================================================#


  # retrieve the NegBin probability parameter
  p = m/(m+r)

  # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
  HC = apply(as.matrix(p), 1, HermCoefNegBin, r=r)

  # ARMA autocorrelation function
  ar.acf <- ARMAacf(ar = phi, lag.max = n)


  # create a grid for nonstationary covariance arguments
  arg_df <- expand.grid(0:(n-1), 0:(n-1))

  # compute the covariance for all values of t1,t2 = 0,...,n-1
  gamma = CountACVF_2(arg_df[,1], arg_df[,2], myacf = ar.acf, g = HC)
  GAMMA = matrix(gamma,n,n)

  # GAMMA = matrix(NA,n,n)
  # # Final toeplitz covariance matrix when covariates are present
  # for (t2 in 0:(n-1)){
  #    GAMMA[t2+1,] = CountACVF_2(t1=0:(n-1), t2, myacf = ar.acf, g = HC)
  # }
  return(GAMMA)
}


#---------Covariance matrix version 3---------#
CovarNegBinAR_3 = function(n, r, m, phi){
  # retrieve the NegBin probability parameter
  p = m/(m+r)

  # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
  HC = apply(as.matrix(p), 1, HermCoefNegBin, r=r)

  # ARMA autocorrelation function
  ar.acf <- ARMAacf(ar = phi, lag.max = n)
  All.ar = apply(as.matrix(ar.acf), 1,function(x)x^(1:20))

  G = matrix(NA,n,n)
  k = 1:20
  kfac = factorial(k)

  for(t1 in 0:(n-1)){
    for(t2 in 0:t1 ){
      h = abs(t1-t2)+1
      G[t1+1,t2+1]= sum(kfac*HC[,t1+1]*HC[,t2+1]*All.ar[,h])
    }
  }
  G = symmetrize(G, update.upper=TRUE)
  return(G)
}


#---------Covariance matrix version 4---------#
CovarNegBinAR_4 = function(n, r, m, phi){
  # retrieve the NegBin probability parameter
  p = m/(m+r)

  # truncation numbe: check me
  N = sapply(unique(p),function(x)which(pnbinom(1:1000, r,1-x)>1-1e-7)[1])
  #N = sapply(unique(p),function(x)which(round(pnbinom(1:1000, r,x), 7) == 1)[1] )
  N[is.na(N)] = 1000

  # Nall = rep(NA,n)
  # Nall[unique(p)[1]==p] = N[1]
  # Nall[unique(p)[2]==p] = N[2]

  # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
  HC = matrix(NA,20,n)
  for(i in  1:length(unique(p))){
    HC[,unique(p)[i]==p] = HermCoefNegBin_2(r, unique(p)[i], N[i])
  }

  # HCsmall = mapply(function(x,y)HermCoefNegBin_2(r,x,y),unique(p),N)
  #
  # HC[,unique(p)[1]==p] = HCsmall[,1]
  # HC[,unique(p)[2]==p] = HCsmall[,2]



  # ARMA autocorrelation function
  ar.acf <- ARMAacf(ar = phi, lag.max = n)
  All.ar = apply(as.matrix(ar.acf), 1,function(x)x^(1:20))

  G = matrix(NA,n,n)
  k = 1:20
  kfac = factorial(k)

  for(t1 in 0:(n-1)){
    for(t2 in 0:t1 ){
      h = abs(t1-t2)+1
      G[t1+1,t2+1]= sum(kfac*HC[,t1+1]*HC[,t2+1]*All.ar[,h])
    }
  }
  G = symmetrize(G, update.upper=TRUE)
  return(G)
}


#---------autoCovariance function---------#
CountACVF_t1t2 = function(t1,t2, myacf, g){
  #====================================================================================#
  # PURPOSE    Compute the autocovariance matrix of the count series.
  #            See relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  # INPUT
  #   h        acvf lag
  #   myacf    autocorrelation of Gaussian series: rho_z
  #   g        Hermitte Coefficients
  # Output
  #   gamma_x  count acvf
  #
  # Notes:     See also CountACVF function--vectorized version of CountACVF_h
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #====================================================================================#

  k = nrow(g) #check me
  gamma_x = sum(g[,t1+1]*g[,t2+1] *  factorial(1:k) * (myacf[abs(t1-t2)+1])^(1:k))
  return(gamma_x)

}


#---------vectorized autoCovariance function---------#
CountACVF_2 <- Vectorize(CountACVF_t1t2, vectorize.args = c("t1", "t2"))





#---------simulate negbin series (r,p) parametrization---------#
sim_negbin_ma = function(n, theta, r,p){
  #====================================================================================#
  # PURPOSE   Simulate Poisson series with MA structure. See relation (1)
  #           in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   n        series length
  #   theta    MA parameter
  #   r,p      Marginal Parameters
  #
  # Output
  #   x        Poisson series
  #
  # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date       April 2020
  # Version    3.6.3
  #====================================================================================#
  z = arima.sim(model = list(order = c(0,0,1), ma=theta), n = n); z = z/sd(z) # standardized
  x = qnbinom(pnorm(z), r, 1-p)
  return(x)
}

#---------Covariance matrix---------#
CovarNegBinMA = function(n,r, p, theta){
  #====================================================================================#
  # PURPOSE    Compute the covariance matrix of a NegBin MA series.
  #
  # INPUT
  #   r,p      Marginal parameters
  #   phi      MA parameter
  #   n        size of the matrix
  #
  # Output
  #   GAMMA    covariance matrix ofcount series
  #
  # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date       April 2020
  # Version    3.6.3
  #====================================================================================#

  # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
  HC = HermCoefNegBin(r,p)

  # ARMA autocorrelation function
  ma.acf <- ARMAacf(ma = theta, lag.max = n)

  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  gamma_x = CountACVF(h = 0:(n-1), myacf = ma.acf, g = HC)

  # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = toeplitz(gamma_x)
  return(GAMMA)
}

#---------Likelihood function---------#
GaussLogLikNB_MA = function(theta, data){
  #====================================================================================#
  # PURPOSE    Compute Gaussian log-likelihood for NegBin MA series
  #
  # INPUT
  #   theta    parameter vector containing the marginal and MA parameters
  #   data     count series
  #
  # Output
  #   loglik   Gaussian log-likelihood
  #
  # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date       April 2020
  # Version    3.6.3
  #====================================================================================#

  # retrieve parameters and sample size
  r = theta[1]
  p = theta[2]
  MAParm = theta[-c(1,2)]
  n = length(data)

  #Select the mean value used to demean--sample or true?
  MeanValue = r*p/(1-p)

  # assign large likelihood value if not causal or if meanm outside range
  if(abs(MAParm) > 0.99 || r<0.0001 || p<0.0001 || p>0.999){
    return(10^16) #check me
  }


  # assign large likelihood value if not causal or if meanm outside range
  # if(abs(MAParm)>0.999 || MeanValue<0 ){
  #   return(NA) #check me
  # }

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = CovarNegBinMA(n, r, p, MAParm)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}

#---------Fit Gaussian Likelihood function---------#
FitGaussianLikNB_MA = function(initialParam, x){
  #====================================================================================#
  # PURPOSE    Fit the Gaussian log-likelihood for NegBin MA series
  #
  # INPUT
  #   initialParam       parameter vector containing the marginal and MA parameters
  #   x                  count series
  #
  # Output
  #   optim.output$par   parameter estimates
  #
  # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date       April 2020
  # Version    3.6.3
  #====================================================================================#
  optim.output <- optim(par     = initialParam,
                        fn      = GaussLogLikNB_MA,
                        data    = x,
                        method  = "BFGS",
                        hessian = TRUE)

  nparms  = length(initialParam)
  ParmEst = matrix(0,nrow=1,ncol=nparms)
  se      = matrix(NA,nrow=1,ncol=nparms)
  loglik  = rep(0,1)

  # save estimates, loglik and standard errors
  ParmEst[,1:nparms]   = optim.output$par
  loglik               = optim.output$value
  #se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))

  All      = cbind(ParmEst, se, loglik)
  return(All)

}




#---------simulate negbin series (r,p) parametrization---------#
sim_negbin = function(n, ARMAmodel, r,p){
  #====================================================================================#
  # PURPOSE       Simulate NegBin series with ARMA structure. See relation (1)
  #               in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   n           series length
  #   ARMAmodel   list with ARMA parameters
  #   r,p         Marginal Parameters
  #
  # Output
  #   X           Neg bin series
  #
  # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date          April 2020
  # Version       3.6.3
  #====================================================================================#

  z = arima.sim(model = list(ar=ARMAmodel[[1]], ma=ARMAmodel[[2]]), n = n); z = z/sd(z)
  X = qnbinom(pnorm(z), r, 1-p)
  return(X)
}

#---------Covariance matrix---------#
CovarNegBin = function(n, r, p, AR, MA, N, nHC){
  #====================================================================================#
  # PURPOSE     Compute the covariance matrix of a NegBin series.
  #
  # INPUT
  #   r,p       Marginal parameters
  #   AR,MA     ARMA parameters
  #   n         size of the matrix
  #   N         truncation for relation (21)
  #   nHC       number of Hermitte coefficents to be computed
  #
  # Output
  #   GAMMA     covariance matrix ofcount series
  #
  # Authors     Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date        April 2020
  # Version     3.6.3
  #====================================================================================#

  # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
  HC = HermCoefNegBin(r,p,N, nHC)

  # ARMA autocorrelation function
  if(is.na(AR)){arma.acf <- ARMAacf(ma = MA, lag.max = n)}
  if(is.na(MA)){arma.acf <- ARMAacf(ar = AR, lag.max = n)}
  if(!is.na(AR) & !is.na(MA)){arma.acf <- ARMAacf(ar = AR, ma = MA, lag.max = n)}

  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  gamma_x = CountACVF(h = 0:(n-1), myacf = arma.acf, g = HC)

  # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = toeplitz(gamma_x)
  return(GAMMA)
}

#---------Likelihood function---------#
GaussLogLikNB = function(theta, data, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE      Compute Gaussian log-likelihood for NegBin series
  #
  # INPUT
  #   theta      parameter vector containing the marginal and ARMA parameters
  #   data       count series
  #   ARMAorder  ordeer of ARMA model
  #   MaxCdf     cdf will be computed up to this number (for light tails cdf=1 fast)
  #   nHC        number of HC to be computed
  #
  #
  # Output
  #   loglik     Gaussian log-likelihood
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         April 2020
  # Version      3.6.3
  #====================================================================================#

  # retrieve parameters and sample size
  r = theta[1]
  p = theta[2]
  AR = ifelse(ARMAorder[1]>0, theta[3: (2+ARMAorder[1])], NA)
  MA = ifelse(ARMAorder[2]>0, theta[(2+ARMAorder[1]+1) : length(theta)], NA)
  n = length(data)

  # compute truncation of relation (21)
  N <- which(round(pnbinom(1:MaxCdf, r,1-p), 7) == 1)[1]
  # if(length(N)==0 |is.na(N) ){
  #   cat(sprintf("The max cdf value is %f and N=%f", max(round(pnbinom(1:MaxCdf, r,1-p), 7)),N))
  #   stop("Haven't reached upper limit for cdf")
  # }
  if(length(N)==0 |is.na(N) ){
    N =MaxCdf
  }

  #Select the mean value used to demean--sample or true?
  MeanValue = r*p/(1-p)


  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = CovarNegBin(n, r, p, AR, MA, N, nHC)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}

#---------Fit Gaussian Likelihood function---------#
FitGaussianLikNB = function(x0, X, LB, UB, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE       Fit the Gaussian log-likelihood for NegBin series
  #
  # INPUT
  #   x0          initial parameters
  #   X           count series
  #   LB          parameter lower bounds
  #   UB          parameter upper bounds
  #   ARMAorder   order of the udnerlying ARMA model
  #   MaxCdf      cdf will be computed up to this number (for light tails cdf=1 fast)
  #   nHC         number of HC to be computed
  #
  # OUTPUT
  #   All         parameter estimates, standard errors, likelihood value
  #
  # NOTES         I may comment out se in cases where maximum is achieved at the boundary
  #
  # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date          April 2020
  # Version       3.6.3
  #====================================================================================#
  optim.output <- optim(par       = x0,
                        fn        = GaussLogLikNB,
                        data      = X,
                        ARMAorder = ARMAorder,
                        MaxCdf    = MaxCdf,
                        nHC       = nHC,
                        method    = "L-BFGS-B",
                        hessian   = TRUE,
                        lower     = LB,
                        upper     = UB
                        )

  nparms  = length(x0)
  ParmEst = matrix(0,nrow=1,ncol=nparms)
  se      = matrix(NA,nrow=1,ncol=nparms)
  loglik  = rep(0,1)

  # save estimates, loglik and standard errors
  ParmEst[,1:nparms]   = optim.output$par
  loglik               = optim.output$value
  #se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))

  All      = cbind(ParmEst, se, loglik)
  return(All)

}


#---------Likelihood function---------#
GaussLogLikNB_Reg = function(theta, data, Regressor, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE      Compute Gaussian log-likelihood for NegBin AR series
  #
  #
  # NOTES        Here I use the following parametrization:
  #
  #              log(mu) = b0 + b1X_i, with
  #              E(Y/X) = mu, Var(Y/X) = mu + kmu^2,
  #              The parameters that enter the likelihood are
  #              b0, b1, k, and the AR parameters. The parameters k and mu are relted with
  #              the classic Neg Bin parameters r and p via:
  #              k =1/r and p = k*mu/(1+k*mu). Note that mu and p depend on each obs but k
  #              and r do not.
  #
  # INPUT
  #   theta      parameter vector containing the marginal and AR parameters
  #   data       count series
  #   Regressor  regressors (first column is always 1--intercept)
  #   ARorder    AR order
  #   M          truncation ofrelation (67) in https://arxiv.org/pdf/1811.00203.pdf
  #
  # Output
  #   loglik     Gaussian log-likelihood
  #
  # Authors      Stefanos Kechagias, James Livsey
  # Date         April 2020
  # Version      3.6.3
  #====================================================================================#

  # retrieve parameters and sample size
  nparms = length(theta)
  nMargParms = nparms - sum(ARMAorder)
  beta       = theta[1:(nparms-sum(ARMAorder)-1)]
  k          = theta[nparms-sum(ARMAorder)]
  AR = ifelse(ARMAorder[1]>0, theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ], NA)
  MA = ifelse(ARMAorder[2]>0, theta[ (length(theta) - ARMAorder[2]) : length(theta)], NA)
  n  = length(data)

  # retrieve mean
  m = exp(Regressor%*%beta)

  # retrieve neg binomial parameters
  r = 1/k
  p = k*m/(1+k*m)

  # Compute truncation of relation (21) in arxiv
  N = sapply(unique(p),function(x)which(pnbinom(1:MaxCdf, r,1-x)>=1-1e-7)[1])-1
  if(length(N)==0 |is.na(N) ){
    N =MaxCdf
  }

  #Select the mean value used to demean--sample or true?
  MeanValue = r*p/(1-p)

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = CovarNegBin(n, r, p, AR, MA, N, nHC)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}


#---------Fit Gaussian Likelihood function---------#
FitGaussianLikNB_Reg = function(x0, X, Regrsessor, LB, UB, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE       Fit the Gaussian log-likelihood for NegBin series
  #
  # INPUT
  #   x0          initial parameters
  #   X           count series
  #   LB          parameter lower bounds
  #   UB          parameter upper bounds
  #   ARMAorder   order of the udnerlying ARMA model
  #   MaxCdf      cdf will be computed up to this number (for light tails cdf=1 fast)
  #   nHC         number of HC to be computed
  #
  # OUTPUT
  #   All         parameter estimates, standard errors, likelihood value
  #
  # NOTES         I may comment out se in cases where maximum is achieved at the boundary
  #
  # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date          April 2020
  # Version       3.6.3
  #====================================================================================#
  optim.output <- optim(par       = x0,
                        fn        = GaussLogLikNB,
                        data      = X,
                        Regressor = Regressor,
                        ARMAorder = ARMAorder,
                        MaxCdf    = MaxCdf,
                        nHC       = nHC,
                        method    = "L-BFGS-B",
                        hessian   = TRUE,
                        lower     = LB,
                        upper     = UB
  )

  nparms  = length(x0)
  ParmEst = matrix(0,nrow=1,ncol=nparms)
  se      = matrix(NA,nrow=1,ncol=nparms)
  loglik  = rep(0,1)

  # save estimates, loglik and standard errors
  ParmEst[,1:nparms]   = optim.output$par
  loglik               = optim.output$value
  #se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))

  All      = cbind(ParmEst, se, loglik)
  return(All)

}






# #---------Likelihood function with regressor---------#
# GaussLogLikNB_Reg = function(theta, Y, X, ARorder, M){
#
#
#   # Number of parameters and sample size
#   nparms = length(theta)
#   n = length(Y)
#
#   # retrieve parameters
#   beta    = theta[1:(nparms-ARorder-1)]
#   k       = theta[nparms-ARorder]
#   phi     = theta[(nparms-ARorder+1):nparms]
#
#
#   # retrieve mean
#   m = exp(X%*%beta)
#
#   # retrieve neg binomial parameters
#   r = 1/k
#   p = k*m/(1+k*m)
#
#   # assign large likelihood value if not causal
#   if(any(abs( polyroot(c(1, -phi))  ) < 1)){
#     return(NA) #check me
#   }
#
#   # Compute the covariance matrix--relation (67) in https://arxiv.org/pdf/1811.00203.pdf
#   GAMMA = CovarNegBinAR_Reg(n, r, p, phi, M)
#
#   # Compute the logdet and the quadratic part
#   logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(Y), m)
#
#   # final loglikelihood value
#   out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]
#
#   # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
#   # out = -2*dmvnorm(as.numeric(Y), rep(lam, n), GAMMA, log = TRUE)
#   return(out)
# }








