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
HermCoefNegBin_Reg <- function(MargParms, ConstMargParm, DynamMargParm, N, nHC, mycdf, nreg){
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
    HC[i] <- HermCoefNegBin_Reg_k(MargParms, ConstMargParm, DynamMargParm , k = i, N, mycdf, nreg)
  }

  return(HC)

}


#---------Hermitte Coefficients for one k---------#
HermCoefNegBin_Reg_k <- function(MargParms, ConstMargParm, DynamMargParm, k, N, mycdf, nreg){
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
  if(nreg>0){
    terms <- exp((-qnorm(mycdf(0:max(N), ConstMargParm,DynamMargParm))^2)/2) * her(qnorm(mycdf(0:max(N), ConstMargParm,DynamMargParm)))
  }else{
    terms <- exp((-qnorm(mycdf(0:max(N), MargParms))^2)/2) * her(qnorm(mycdf(0:max(N), MargParms)))
  }
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


#----------Link coefficients with one dummy Regressor---------#
LinkCoef_Reg = function(ConstMargParm, DynamMargParm, N, nHC, mycdf,nreg){
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
  g1 = HermCoefNegBin_Reg(MargParms = NULL, ConstMargParm, unique(DynamMargParm)[1],N[1], nHC, mycdf, nreg)
  g2 = HermCoefNegBin_Reg(MargParms = NULL, ConstMargParm, unique(DynamMargParm)[2],N[2], nHC, mycdf, nreg)
  HC = cbind(g1, g2)

  # Compute the products g1^2, g1*g2 and g2^2 of the HCs
  HCprod = matrix(NA, nHC, 3)
  for(i in 0:(length(unique(DynamMargParm))-1) ){
    for(j in i: (length(unique(DynamMargParm))-1) ){
      HCprod[,i+j+1] = HC[,i+1]*HC[,j+1]
    }
  }

  return(HCprod)
}


#---------Covariance matrix with one dummy Regressor---------#
CovarNegBin_Reg = function(n, ConstMargParm, DynamMargParm, AR, MA, N, nHC,mycdf, nreg){
  #====================================================================================#
  # PURPOSE    Compute the covariance matrix of a NegBin series that
  #            includes one dummy variable as a regressor. Here p depends on each
  #            observation. See relation (67) in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   r, p     NB MArginal parameters
  #   AR,MA    AR, MA parameters
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


  if(!length(AR) && !length(MA)){
    All.arma = apply(as.matrix( c(1,rep(0,n-1))), 1, function(x)x^(1:nHC))
  }
  # Compute ARMA autocorrelation function
  if(!length(AR) && length(MA)){
    All.arma = apply(as.matrix(ARMAacf(ma = MA, lag.max = n)), 1,function(x)x^(1:nHC))
  }
  if(!length(MA) && length(AR)){
    All.arma = apply(as.matrix(ARMAacf(ar = AR, lag.max = n)), 1,function(x)x^(1:nHC))
  }
  if(length(AR) & length(MA)){
    All.arma = apply(as.matrix(ARMAacf(ar = AR,ma=MA, lag.max = n)), 1,function(x)x^(1:nHC))
  }


  # Compute the link coefficients l_k = factorial(k)*g_{t1,k}*g_{t2,k}
  linkCoef = LinkCoef_Reg(ConstMargParm, DynamMargParm, N, nHC,mycdf,nreg)

  # keep track of which indices each unique HC is located at
  index = cbind(unique(DynamMargParm)[1]==DynamMargParm, unique(DynamMargParm)[2]==DynamMargParm)

  # STEP 3: Implement relation 67
  k = 1:nHC
  kfac = factorial(k)
  G = matrix(NA,n,n)
  for(t1 in 0:(n-1)){
    for(t2 in 0:t1 ){
      h = abs(t1-t2)+1
      G[t1+1,t2+1]= sum(kfac*linkCoef[,index[t1+1,2]+index[t2+1,2]+1]*All.arma[,h])
    }
  }
  G = symmetrize(G, update.upper=TRUE)
  return(G)
}


#---------Guassian Likelihood function with a dummy Regressor---------#
GaussLogLikNB_Reg = function(theta, data, Regressor, ARMAorder, MaxCdf, nHC, CountDist){
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



  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # retrieve sample size
  n = length(data)

  # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments

  # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1:(1+nreg),
                           "Negative Binomial"   = 1:(2+nreg),
                           "Mixed Poisson"       = 1:(3+nreg),
                           "Generalized Poisson" = 1:(2+nreg),
                           "Binomial"            = 1:(2+nreg))
  if(nreg<1){
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = ppois,
                   "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
                   "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = pGpois,
                   "Binomial"            = pbinom
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = dGpois,
                   "Binomial"            = dbinom
    )
  }else{
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = function(x, ConstMargParm, DynamMargParm){             ppois   (x, DynamMargParm)},
                   "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ pnbinom (x, ConstMargParm, 1-DynamMargParm)},
                   "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ pGpois  (x, ConstMargParm, DynamMargParm)}
    )
    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = function(x, ConstMargParm, DynamMargParm){             dpois   (x, DynamMargParm)},
                   "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ dnbinom (x, ConstMargParm, 1-DynamMargParm)},
                   "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ dGpois  (x, ConstMargParm, DynamMargParm)}
    )
  }



  # retrieve marginal distribution parameters
  MargParms  = theta[MargParmIndices]
  nMargParms = length(MargParms)
  nparms     = length(theta)

  # check if the number of parameters matches the model setting
  if(nMargParms + sum(ARMAorder)!=nparms) stop('The length of theta does not match the model specification.')


  # Add the regressor to the parameters--only works for Negative Binomial
  if(CountDist == "Negative Binomial" && nreg>0){
    beta = MargParms[1:(nreg+1)]
    k = MargParms[nreg+2]
    m = exp(Regressor%*%beta)
    r = 1/k
    p = k*m/(1+k*m)
    ConstMargParm = r
    DynamMargParm = p
  }

  if(CountDist == "Generalized Poisson" && nreg>0){
    beta   = MargParms[1:(nreg+1)]
    ConstMargParm  = MargParms[nreg+2]
    DynamMargParm = exp(Regressor%*%beta)
  }

  if(CountDist == "Poisson" && nreg>0){
    beta           = MargParms[1:(nreg+1)]
    ConstMargParm  = NULL
    DynamMargParm  = exp(Regressor%*%beta)
  }


  # retrieve ARMA parameters
  AR = NULL
  if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]

  MA = NULL
  if(ARMAorder[2]>0) MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]


  if(!is.null(AR) && is.null(MA)){
    if(any(abs( polyroot(c(1, -AR))  ) < 1)){
      return(10^6) #check me
    }
  }

  if(!is.null(MA) && is.null(AR)){
    if(any(abs( polyroot(c(1, -MA))  ) < 1)){
      return(10^6) #check me
    }
  }

  if(!is.null(MA) && !is.null(AR)){
    if(checkPoly(AR,MA)[1]!="Causal" && check(poly)[2]!="Invertible"){
      return(10^6) # check me--do i need inveritibility?
    }
  }


  # retrieve mean
  if(nreg>0){
    MeanValue = exp(Regressor%*%beta)
  }else{
    MeanValue = switch(CountDist,
                       "Poisson"             = MargParms[1],
                       "Negative Binomial"   = MargParms[1]*MargParms[2]/(1-MargParms[2]),
                       "Generalized Poisson" = MargParms[2])
  }

  # Compute truncation of relation (21)
  if(nreg>0){
    N = sapply(unique(DynamMargParm),function(x)which(mycdf(1:MaxCdf, ConstMargParm, x)>=1-1e-7)[1])-1
    N[is.na(N)] = MaxCdf
  }else{
    N <- which(round(mycdf(1:MaxCdf, MargParms), 7) == 1)[1]
    if(length(N)==0 |is.na(N) ){
      N =MaxCdf
    }
  }


  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  # GAMMA = CovarNegBin(n, r, p, AR, MA, N, nHC)
  if (nreg>0){
    GAMMA = CovarNegBin_Reg(n, ConstMargParm, DynamMargParm, AR, MA, N, nHC, mycdf, nreg)
  }else{
    GAMMA = CovarNegBin(n, MargParms, AR, MA, N, nHC, mycdf, nreg)
  }


  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}


#----------Link coefficients with one dummy Regressor---------#
LinkCoef_Reg2 = function(ConstMargParm, DynamMargParm, N, nHC, mycdf){
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
  g1 = HermCoefNegBin_Reg(MargParms = NULL, unique(DynamMargParm)[1], ConstMargParm,N[1], nHC, mycdf, nreg)
  g2 = HermCoefNegBin_Reg(MargParms = NULL, unique(DynamMargParm)[2], ConstMargParm,N[2], nHC, mycdf, nreg)
  HC = cbind(g1, g2)

  # Compute the products g1^2, g1*g2 and g2^2 of the HCs
  HCprod = matrix(NA, nHC, 3)
  for(i in 0:(length(unique(DynamMargParm))-1) ){
    for(j in i: (length(unique(DynamMargParm))-1) ){
      HCprod[,i+j+1] = HC[,i+1]*HC[,j+1]
    }
  }

  return(HCprod)
}


#---------Covariance matrix with one dummy Regressor---------#
CovarNegBinAR_Reg2 = function(n, ConstMargParm, DynamMargParm, phi, N, nHC,mycdf){
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

  if(!length(AR) && !length(MA)){
    All.arma = apply(as.matrix( c(1,rep(0,n-1))), 1, function(x)x^(1:nHC))
  }
  # Compute ARMA autocorrelation function
  if(!length(AR) && length(MA)){
    All.arma = apply(as.matrix(ARMAacf(ma = MA, lag.max = n)), 1,function(x)x^(1:nHC))
  }
  if(!length(MA) && length(AR)){
    All.arma = apply(as.matrix(ARMAacf(ar = AR, lag.max = n)), 1,function(x)x^(1:nHC))
  }
  if(length(AR) & length(MA)){
    All.arma = apply(as.matrix(ARMAacf(ar = AR,ma=MA, lag.max = n)), 1,function(x)x^(1:nHC))
  }

  # Compute the link coefficients l_k = factorial(k)*g_{t1,k}*g_{t2,k}
  linkCoef = LinkCoef_Reg2(ConstMargParm, DynamMargParm, N, nHC, mycdf)

  # keep track of which indices each unique HC is located at
  index = cbind(unique(DynamMargParm)[1]==DynamMargParm, unique(DynamMargParm)[2]==DynamMargParm)

  # STEP 3: Implement relation 67
  k = 1:nHC
  kfac = factorial(k)
  G = matrix(NA,n,n)
  for(t1 in 0:(n-1)){
    for(t2 in 0:t1 ){
      h = abs(t1-t2)+1
      G[t1+1,t2+1]= sum(kfac*linkCoef[,index[t1+1,2]+index[t2+1,2]+1]*All.arma[,h])
    }
  }
  G = symmetrize(G, update.upper=TRUE)
  return(G)
}


#---------Guassian Likelihood function with a dummy Regressor---------#
GaussLogLikNB_Reg2 = function(theta, data, Regressor, ARMAorder, MaxCdf, nHC, CountDist){
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

  # retrieve marginal cdf
  mycdf = switch(CountDist,
                 "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ pnbinom (x, ConstMargParm, 1-DynamMargParm)},
                 "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ pGpois  (x, ConstMargParm, DynamMargParm)}
  )

  # retrieve parameters and sample size
  nparms     = length(theta)
  nreg       = dim(Regressor)[2]-1
  nMargParms = nparms - sum(ARMAorder)
  beta       = theta[1:(nreg+1)]
  k          = theta[nparms-sum(ARMAorder)]
  n          = length(data)


  AR = NULL
  if(ARMAorder[1]>0) AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])]

  MA = NULL
  if(ARMAorder[2]>0){
    MA = theta[ (nMargParms+ARMAorder[1]+1) : (nMargParms + ARMAorder[1] + ARMAorder[2]) ]
  }

  if(!is.null(AR) && is.null(MA)){
    if(any(abs( polyroot(c(1, -AR))  ) < 1)){
      return(10^6) #check me
    }
  }

  if(!is.null(MA) && is.null(AR)){
    if(any(abs( polyroot(c(1, -MA))  ) < 1)){
      return(10^6) #check me
    }
  }

  if(!is.null(MA) && !is.null(AR)){
    if(checkPoly(AR,MA)[1]!="Causal" && check(poly)[2]!="Invertible"){
      return(10^6) # check me--do i need inveritibility?
    }
  }


  # retrieve mean
  r = exp(Regressor%*%beta)


  # Compute truncation of relation (21) in arxiv
  N = sapply(unique(r),function(x)which(mycdf(1:MaxCdf, x,1-p)>=1-1e-7)[1])-1
  N[is.na(N)] = MaxCdf

  #Select the mean value used to demean--sample or true?
  MeanValue = r*p/(1-p)

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  # GAMMA = CovarNegBin(n, r, p, AR, MA, N, nHC)
  GAMMA = CovarNegBinAR_Reg2(n, r, p, AR, N, nHC, mycdf)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}



#---------wrapper to fit Guassian Likelihood function with a dummy Regressor---------#
FitGaussianLikNB_Reg= function(x0, X, Regressor, CountDist, Particles, LB, UB,
                               ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod){
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

  nparms  = length(x0)
  n = length(X)
  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst = matrix(0,nrow=1,ncol=nparms)
  se =  matrix(NA,nrow=1,ncol=nparms)
  loglik = rep(0,1)
  convcode = rep(0,1)
  kkt1 = rep(0,1)
  kkt2 = rep(0,1)

  if(Model){
    optim.output <- optimx(par      = x0,
                           fn        = GaussLogLikNB_Reg2,
                           data      = X,
                           Regressor = Regressor,
                           ARMAorder = ARMAorder,
                           CountDist = CountDist,
                           MaxCdf    = MaxCdf,
                           nHC       = nHC,
                           hessian   = TRUE,
                           lower     = LB,
                           upper     = UB,
                           method    = OptMethod
    )
  }else{
    optim.output <- optimx(par       = x0,
                           fn        = GaussLogLikNB_Reg,
                           data      = X,
                           Regressor = Regressor,
                           ARMAorder = ARMAorder,
                           CountDist = CountDist,
                           MaxCdf    = MaxCdf,
                           nHC       = nHC,
                           hessian   = TRUE,
                           lower     = LB,
                           upper     = UB,
                           method    = OptMethod
    )
  }

  # save estimates, loglik value and diagonal hessian
  ParmEst  = as.numeric(optim.output[1:nparms])
  loglik   = optim.output$value
  convcode = optim.output$convcode
  kkt1     = optim.output$kkt1
  kkt2     = optim.output$kkt2

  if(Model){
    # compute hessian
    H = gHgen(par            = ParmEst,
              fn             = GaussLogLikNB_Reg2,
              data           = X,
              CountDist      = CountDist,
              Regressor      = Regressor,
              ARMAorder      = ARMAorder,
              MaxCdf         = MaxCdf,
              nHC            = nHC
    )
  }else{
    # compute hessian
    H = gHgen(par            = ParmEst,
              fn             = GaussLogLikNB_Reg,
              data           = X,
              CountDist      = CountDist,
              Regressor      = Regressor,
              ARMAorder      = ARMAorder,
              MaxCdf         = MaxCdf,
              nHC            = nHC
    )
  }

  if(H$hessOK && det(H$Hn)>10^(-12)){
    se = sqrt(abs(diag(solve(H$Hn))))
  }else{
    se = rep(NA, nparms)
  }

  # Compute model selection criteria
  Criteria = ComputeCriteria(loglik, nparms, n, Particles)


  # get the names of the final output
  parmnames = colnames(optim.output)
  mynames = c(parmnames[1:nparms],paste("se", parmnames[1:nparms], sep="_"), "loglik", "AIC", "BIC","AICc", "status", "kkt1", "kkt2")



  All = matrix(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)
  colnames(All) = mynames
  return(All)

}

#---------Compute AIC, BIC, AICc
ComputeCriteria = function(loglik, nparms, n, Particles){
  #---------------------------------#
  # Purpose:     Compute AIC, BIC, AICc for our models
  #
  # Inputs:
  #  loglik      log likelihood value at optimum
  #  nparms      number of parametersin the model
  #       n      sample size
  #
  # NOTES:       The Gaussian likelihood we are minimizing has the form:
  #              l0 = 0.5*logdet(G) + 0.5*X'*inv(G)*X. We will need to
  #              bring this tothe classic form before we compute the
  #              criteria (see log of lrealtion (8.6.1) in Brockwell and Davis)
  #
  #              No correction is necessary in the PF case
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         July 2020
  # Version      3.6.3
  #---------------------------------#

  if(!is.null(Particles)){
    l1 = -loglik
  }else{
    l1 = -loglik  - n/2*log(2*pi)
  }

  AIC = 2*nparms - 2*l1
  BIC = log(n)*nparms - 2*l1
  AICc = AIC + (2*nparms^2 + 2*nparms)/(n-nparms-1)

  AllCriteria = c(AIC, BIC, AICc)
}


#---------initial estimates via method of moments and reversion
ComputeInitNegBinMA = function(x,n,nHC,LB,UB, mycdf, nreg){
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
  initParms = ComputeInitNegBinMAterm(x, MargParms, ConstMargParm, DynamMargParm, n, nHC,LB,UB, mycdf, nreg)

  return(initParms)

}


#---------link coefficients
link_coefs <- function(g_coefs, gamx0){
  K <- length(g_coefs)
  return(factorial(1:K) * g_coefs^2 / gamx0)
}

#--------obtain initial estimate for MA term using acf and reversion
ComputeInitNegBinMAterm = function(x, MargParms, ConstMargParm, DynamMargParm, N, nHC, LB, UB, mycdf, nreg){

  # compute Hermite coefficients
  g.coefs  = HermCoefNegBin_Reg(MargParms, ConstMargParm=NULL, DynamMargParm=NULL, N, nHC, mycdf, nreg)

  # Compute acf of count series at lag 0
  NegBinVar = MargParms[1]*MargParms[2]/(1-MargParms[2])^2

  # compute link coeffcients
  link.coefs <- link_coefs(g.coefs, NegBinVar)

  # compute Inverse Link coefficients of f^-1: gam.z --> gam.x
  inv.link.coefs <- reversion(link.coefs)

  # sample acf of count data
  gam.x <- acf(x = x, lag.max = 30, plot = FALSE, type = "correlation")$acf

  # compute gamma Z thru reversion
  gam.z <- power_series(gam.x[,,1], inv.link.coefs)
  #gam.z = gam.z/gam.z[1]

  thetaEst = gam.z[2]

  # correct if I am outside the boundaries
  if(thetaEst<LB[3]){thetaEst = 1.1*LB[3]}
  if(thetaEst>UB[3]){thetaEst = 0.9*UB[3]}

  InitEstimates = c(rEst,pEst,thetaEst)
  return(InitEstimates)
}

#-------- check causality and invertible
checkPoly = function(AR, MA) {
  #---------------------------------#
  # Purpose: Cehck cauality and invertibility
  #
  # Note: Thje function is modified from itsmr fnc check()
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         July 2020
  # Version      3.6.3
  #---------------------------------#
  c1 = "Causal"
  c2 = "Invertible"

  if (!is.null(AR))  c1 = ifelse(all(abs(polyroot( c(1,-AR))) > 1),"Causal", "NonCausal" )

  if (!is.null(MA))  c2 = ifelse(all(abs(polyroot( c(1,MA))) > 1),"Invertible", "NonInvertible" )

  return(c(c1,c2))
}



#---------Covariance matrix---------#
CovarNegBin = function(n, MargParms, AR, MA, N, nHC,mycdf,nreg){
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
  HC = HermCoefNegBin_Reg(MargParms, ConstMargParm=NULL, DynamMargParm=NULL, N, nHC, mycdf, nreg)

  # ARMA autocorrelation function
  if(!length(AR) && length(MA)  ){arma.acf <- ARMAacf(ma = MA, lag.max = n)}
  if(!length(MA) && length(AR)){arma.acf <- ARMAacf(ar = AR, lag.max = n)}
  if(length(AR)  && length(MA)){arma.acf <- ARMAacf(ar = AR, ma = MA, lag.max = n)}


  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  if(is.null(AR) && is.null(MA)){
    gamma_x = c(as.numeric(HC^2%*%factorial(1:nHC)), rep(0,n-1))
  }
  else{
    gamma_x = CountACVF(h = 0:(n-1), myacf = arma.acf, g = HC)
  }
  # Final toeplitz covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = toeplitz(gamma_x)
  return(GAMMA)
}


#---------Gaussian Likelihood function---------#
# GaussLogLikNB = function(theta, data, ARMAorder, MaxCdf, nHC){
#   #====================================================================================#
#   # PURPOSE      Compute Gaussian log-likelihood for NegBin series
#   #
#   # INPUT
#   #   theta      parameter vector containing the marginal and ARMA parameters
#   #   data       count series
#   #   ARMAorder  ordeer of ARMA model
#   #   MaxCdf     cdf will be computed up to this number (for light tails cdf=1 fast)
#   #   nHC        number of HC to be computed
#   #
#   #
#   # Output
#   #   loglik     Gaussian log-likelihood
#   #
#   # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date         April 2020
#   # Version      3.6.3
#   #====================================================================================#
#
#   # retrieve parameters and sample size
#   r = theta[1]
#   p = theta[2]
#   nparms  = length(theta)
#   nMargParms = nparms - sum(ARMAorder)
#   if(ARMAorder[1]>0){
#     AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
#   }else{
#     AR = NULL
#   }
#
#   if(ARMAorder[2]>0){
#     MA = theta[ (nMargParms+ARMAorder[1]+1) : length(theta)]
#   }else{
#     MA = NULL
#   }
#   n = length(data)
#
#   # compute truncation of relation (21)
#   N <- which(round(pnbinom(1:MaxCdf, r,1-p), 7) == 1)[1]
#   # if(length(N)==0 |is.na(N) ){
#   #   cat(sprintf("The max cdf value is %f and N=%f", max(round(pnbinom(1:MaxCdf, r,1-p), 7)),N))
#   #   stop("Haven't reached upper limit for cdf")
#   # }
#   if(length(N)==0 |is.na(N) ){
#     N =MaxCdf
#   }
#
#   #Select the mean value used to demean--sample or true?
#   MeanValue = r*p/(1-p)
#
#   # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
#   GAMMA = CovarNegBin(n, MargParms, AR, MA, N, nHC, mycdf)
#
#   # Compute the logdet and the quadratic part
#   logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)
#
#   # final loglikelihood value
#   out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]
#
#   # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
#   # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
#   return(out)
# }
#

#---------wrapper to fit Gaussian Likelihood function---------#
# FitGaussianLikNB = function(x0, X, LB, UB, ARMAorder, MaxCdf, nHC){
#   #====================================================================================#
#   # PURPOSE       Fit the Gaussian log-likelihood for NegBin series
#   #
#   # INPUT
#   #   x0          initial parameters
#   #   X           count series
#   #   LB          parameter lower bounds
#   #   UB          parameter upper bounds
#   #   ARMAorder   order of the udnerlying ARMA model
#   #   MaxCdf      cdf will be computed up to this number (for light tails cdf=1 fast)
#   #   nHC         number of HC to be computed
#   #
#   # OUTPUT
#   #   All         parameter estimates, standard errors, likelihood value
#   #
#   # NOTES         I may comment out se in cases where maximum is achieved at the boundary
#   #
#   # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
#   # Date          April 2020
#   # Version       3.6.3
#   #====================================================================================#
#   optim.output <- optim(par       = x0,
#                         fn        = GaussLogLikNB,
#                         data      = X,
#                         ARMAorder = ARMAorder,
#                         MaxCdf    = MaxCdf,
#                         nHC       = nHC,
#                         method    = "L-BFGS-B",
#                         hessian   = TRUE,
#                         lower     = LB,
#                         upper     = UB
#                         )
#
#   nparms  = length(x0)
#   ParmEst = matrix(0,nrow=1,ncol=nparms)
#   se      = matrix(NA,nrow=1,ncol=nparms)
#   loglik  = rep(0,1)
#
#   # save estimates, loglik and standard errors
#   ParmEst[,1:nparms]   = optim.output$par
#   loglik               = optim.output$value
#   #se[,1:nparms]        = sqrt(abs(diag(solve(optim.output$hessian))))
#
#   All      = cbind(ParmEst, se, loglik)
#   return(All)
#
# }
