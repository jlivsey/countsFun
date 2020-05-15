#---------simulate negbin series with prob of success p---------#
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
  # The third argument of qnbinom is the prob of failure
  X = qnbinom(pnorm(z), r, 1-p)
  return(X)
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


#---------Hermitte Coefficients for one k---------#
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

  terms[is.nan(terms)]=0

  # take the sum of all terms
  HC_k <- sum(terms) / (sqrt(2*pi) *  factorial(k))
  return(HC_k)
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
  if(!length(AR)){arma.acf <- ARMAacf(ma = MA, lag.max = n)}
  if(!length(MA)){arma.acf <- ARMAacf(ar = AR, lag.max = n)}
  if(length(AR) & length(MA)){arma.acf <- ARMAacf(ar = AR, ma = MA, lag.max = n)}

  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  gamma_x = CountACVF(h = 0:(n-1), myacf = arma.acf, g = HC)

  # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = toeplitz(gamma_x)
  return(GAMMA)
}


#---------Gaussian Likelihood function---------#
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
  nparms  = length(theta)
  nMargParms = nparms - sum(ARMAorder)
  if(ARMAorder[1]>0){
    AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
  }else{
    AR = NULL
  }

  if(ARMAorder[2]>0){
    MA = theta[ (nMargParms+ARMAorder[1]+1) : length(theta)]
  }else{
    MA = NULL
  }
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


#---------wrapper to fit Gaussian Likelihood function---------#
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


#----------Link coefficients with one dummy Regressor---------#
LinkCoef_Reg = function(r, p, N, nHC){
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
  g1 = HermCoefNegBin(r, unique(p)[1],N[1], nHC)
  g2 = HermCoefNegBin(r, unique(p)[2],N[2], nHC)
  HC = cbind(g1, g2)

  # Compute the products g1^2, g1*g2 and g2^2 of the HCs
  HCprod = matrix(NA, nHC, 3)
  for(i in 0:(length(unique(p))-1) ){
    for(j in i: (length(unique(p))-1) ){
      HCprod[,i+j+1] = HC[,i+1]*HC[,j+1]
    }
  }

  return(HCprod)
}


#---------Covariance matrix with one dummy Regressor---------#
CovarNegBinAR_Reg = function(n, r, p, phi, N, nHC){
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
  All.ar = apply(as.matrix(ARMAacf(ar = phi, lag.max = n)), 1,function(x)x^(1:nHC))

  # Compute the link coefficients l_k = factorial(k)*g_{t1,k}*g_{t2,k}
  linkCoef = LinkCoef_Reg(r, p, N, nHC)

  # keep track of which indices each unique HC is located at
  index = cbind(unique(p)[1]==p, unique(p)[2]==p)

  # STEP 3: Implement relation 67
  k = 1:nHC
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


#---------Guassian Likelihood function with a dummy Regressor---------#
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
  #   ARMAorder  ARMA order
  #   MaxCdf     cdf will be computed up to this number (for light tails cdf=1 fast)
  #   nHC        number of HC to be computed
  #
  # Output
  #   loglik     Gaussian log-likelihood
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         April 2020
  # Version      3.6.3
  #====================================================================================#

  # retrieve parameters and sample size
  nparms     = length(theta)
  nMargParms = nparms - sum(ARMAorder)
  beta       = theta[1:(nparms-sum(ARMAorder)-1)]
  k          = theta[nparms-sum(ARMAorder)]
  n          = length(data)

  if(ARMAorder[1]>0){
    AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
  }else{
    AR = NULL
  }

  if(ARMAorder[2]>0){
    MA = theta[ (length(theta) - ARMAorder[2]) : length(theta)]
  }else{
    MA = NULL
  }

  # need only only causal
  if(any(abs( polyroot(c(1, -AR))  ) < 1)){
    return(10^6) #check me
  }

  # retrieve mean
  m = exp(Regressor%*%beta)

  # retrieve neg binomial parameters
  r = 1/k
  p = k*m/(1+k*m)

  # Compute truncation of relation (21) in arxiv
  N = sapply(unique(p),function(x)which(pnbinom(1:MaxCdf, r,1-x)>=1-1e-7)[1])-1
  N[is.na(N)] = MaxCdf

  #Select the mean value used to demean--sample or true?
  MeanValue = r*p/(1-p)

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  # GAMMA = CovarNegBin(n, r, p, AR, MA, N, nHC)
  GAMMA = CovarNegBinAR_Reg(n, r, p, AR, N, nHC)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}


#----------Link coefficients with one dummy Regressor---------#
LinkCoef_Reg2 = function(r, p, N, nHC){
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
  g1 = HermCoefNegBin(unique(r)[1], p,N[1], nHC)
  g2 = HermCoefNegBin(unique(r)[1], p,N[2], nHC)
  HC = cbind(g1, g2)

  # Compute the products g1^2, g1*g2 and g2^2 of the HCs
  HCprod = matrix(NA, nHC, 3)
  for(i in 0:(length(unique(r))-1) ){
    for(j in i: (length(unique(r))-1) ){
      HCprod[,i+j+1] = HC[,i+1]*HC[,j+1]
    }
  }

  return(HCprod)
}


#---------Covariance matrix with one dummy Regressor---------#
CovarNegBinAR_Reg2 = function(n, r, p, phi, N, nHC){
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
  All.ar = apply(as.matrix(ARMAacf(ar = phi, lag.max = n)), 1,function(x)x^(1:nHC))

  # Compute the link coefficients l_k = factorial(k)*g_{t1,k}*g_{t2,k}
  linkCoef = LinkCoef_Reg2(r, p, N, nHC)

  # keep track of which indices each unique HC is located at
  index = cbind(unique(r)[1]==r, unique(r)[2]==r)

  # STEP 3: Implement relation 67
  k = 1:nHC
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


#---------Guassian Likelihood function with a dummy Regressor---------#
GaussLogLikNB_Reg2 = function(theta, data, Regressor, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE      Compute Gaussian log-likelihood for NegBin AR series
  #
  #
  # NOTES        Here I use the parametrization given in our paper see Section 5
  #
  #
  # INPUT
  #   theta      parameter vector containing the marginal and AR parameters
  #   data       count series
  #   Regressor  regressors (first column is always 1--intercept)
  #   ARMAorder  ARMA order
  #   MaxCdf     cdf will be computed up to this number (for light tails cdf=1 fast)
  #   nHC        number of HC to be computed
  #
  # Output
  #   loglik     Gaussian log-likelihood
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         April 2020
  # Version      3.6.3
  #====================================================================================#

  # retrieve parameters and sample size
  nparms     = length(theta)
  nMargParms = nparms - sum(ARMAorder)
  beta       = theta[1:(nparms-sum(ARMAorder)-1)]
  p          = theta[nparms-sum(ARMAorder)]
  n          = length(data)

  if(ARMAorder[1]>0){
    AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
  }else{
    AR = NULL
  }

  if(ARMAorder[2]>0){
    MA = theta[ (length(theta) - ARMAorder[2]) : length(theta)]
  }else{
    MA = NULL
  }

  # need only only causal
  if(any(abs( polyroot(c(1, -AR))  ) < 1)){
    return(10^6) #check me
  }

  # retrieve mean
  r = exp(Regressor%*%beta)


  # Compute truncation of relation (21) in arxiv
  N = sapply(unique(r),function(x)which(pnbinom(1:MaxCdf, x,1-p)>=1-1e-7)[1])-1
  N[is.na(N)] = MaxCdf

  #Select the mean value used to demean--sample or true?
  MeanValue = r*p/(1-p)

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  # GAMMA = CovarNegBin(n, r, p, AR, MA, N, nHC)
  GAMMA = CovarNegBinAR_Reg2(n, r, p, AR, N, nHC)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}


#---------wrapper to fit Guassian Likelihood function with a dummy Regressor---------#
FitGaussianLikNB_Reg = function(x0, X, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, Model){
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
  #   Model       =1 then r depends on t, p fixed, =0 p depends on t, r is fixed

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
if(Model){
  optim.output <- optim(par       = x0,
                        fn        = GaussLogLikNB_Reg2,
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
}else{
  optim.output <- optim(par       = x0,
                        fn        = GaussLogLikNB_Reg,
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
}

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



#---------initial estimates via method of moments and reversion
ComputeInitNegBinMA = function(x,n,nHC){
  #---------------------------------#
  # Purpose: Method of Moment Initial estimates for negbin MA(1)
  #
  #
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         April 2020
  # Version      3.6.3
  #---------------------------------#

  xbar = mean(x)
  sSquare = var(x)

  # Method of Moments for negBin
  rEst = xbar^2/(sSquare - xbar)
  pEst = 1 - xbar/sSquare

  # compute thetaEst using reversion as in IYW
  initParms = ComputeInitNegBinMAterm(x, rEst, pEst, n, nHC)

  return(initParms)

  }


#---------link coefficients
link_coefs <- function(g_coefs, gamx0){
  K <- length(g_coefs)
  return(factorial(1:K) * g_coefs^2 / gamx0)
}

#--------obtain initial estimate for MA term using acf and reversion
ComputeInitNegBinMAterm = function(x, rEst, pEst, N, nHC){

  # compute Hermite coefficients
  g.coefs  = HermCoefNegBin(rEst, pEst, N, nHC)

  # Compute acf of count series at lag 0
  NegBinVar = rEst*pEst/(1-pEst)^2

  # compute link coeffcients
  link.coefs <- link_coefs(g.coefs, NegBinVar)

  # compute Inverse Link coefficients of f^-1: gam.z --> gam.x
  inv.link.coefs <- reversion(link.coefs)

  # sample acf of count data
  gam.x <- acf(x = x, lag.max = 30, plot = FALSE, type = "correlation")$acf

  # compute gamma Z thru reversion
  gam.z <- power_series(gam.x[,,1], inv.link.coefs)
  gam.z = gam.z/gam.z[1]

  thetaEst = gam.z[2]
  InitEstimates = c(rEst,pEst,thetaEst)
  return(InitEstimates)
}

