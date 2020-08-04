# PURPOSE: Function for Gen Poisson analysis
#
# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    July 2020
#
# VERSION: 3.6.3


#--------- Generalized Poisson pdfs---------#
dGpois = function(y,a,m){
  k = m/(1+a*m)
  return( k^y * (1+a*y)^(y-1) * exp(-k*(1+a*y)-lgamma(y+1)))
}

#--------- Generalized Poisson cdf for one m---------#
pgpois1 = function(x,a,m){
  M = max(x,0)
  bb = dGpois(0:M,a,m)
  cc = cumsum(bb)
  g = rep(0,length(x))
  g[x>=0] = cc[x+1]
  return(g)
}

#--------- Generalized Poisson pdf for multiple m---------#
pGpois = function(x,a,m){

  if (length(m)>1){
    r = mapply(pgpois1, x=x, a, m = m)
  }else{
    r = pgpois1(x, a, m)
  }
  return(r)
}

#--------- Generalized Poisson icdf for 1 m---------#
qgpois1 <- function(p,a,m) {
  check.list <- pGpois(0:100,a,m)
  quantile.vec <- rep(-99,length(p))

  for (i in 1:length(p)) {
    x <- which(check.list>=p[i],arr.ind = T)[1]
    quantile.vec[i] <- x-1
  }
  return(quantile.vec)
}

#--------- Generalized Poisson icdf for many m---------#
qGpois = function(p,a,m){

  if (length(m)>1){
    r = mapply(qgpois1, p=p, a, m = m)
  }else{
    r = qgpois1(p,a,m)
  }
  return(r)
}

#--------- Generate Gen Poisson datas---------#
rGpois = function(n, a,m){
  u = runif(n)
  x = qGpois(u,a, m)
  return(x)
}







#---------simulateGen Poisson series---------#
sim_genpois = function(n, ARMAmodel,a, m){
  #====================================================================================#
  # PURPOSE       Simulate GenPois series with ARMA structure. See relation (1)
  #               in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   n           series length
  #   ARMAmodel   list with ARMA parameters
  #   r,p         Marginal Parameters
  #
  # Output
  #   X           GenPois series
  #
  # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date          April 2020
  # Version       3.6.3
  #====================================================================================#

  z = arima.sim(model = list(ar=ARMAmodel[[1]], ma=ARMAmodel[[2]]), n = n); z = z/sd(z)
  X = qGpois(pnorm(z),a,m)
  return(X)
}


#---------Hermitte Coefficients for all k---------#
HermCoefGenPois <- function(a,m, N, nHC, mycdf){
  #====================================================================================#
  # PURPOSE    Compute all Hermite Coefficients. See relation (21) in
  #            https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   a,m      Marginal parameters
  #   maxCoef  number of coefficients to return. Default = 20
  #   N        truncation of relation (21)
  #   nHC      number of Hermitte coefficients
  #
  # Output
  #   HC       All Hermite coeficients
  #
  # Authors    Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date       July 2020
  # Version    3.6.3
  #====================================================================================#

  h = 1:nHC #check me
  HC = rep(NA, length(h)) # storage
  for(i in h) {
    HC[i] <- HermCoefGenPois_k(a,m , k = i, N, mycdf)
  }

  return(HC)

}


#---------Hermitte Coefficients for one k---------#
HermCoefGenPois_k <- function(a,m, k, N, mycdf){
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
  # Date       July 2020
  # Version    3.6.3
  #====================================================================================#

  # function for kth Hermite Polynomial
  her <- function(x){
    evalHermPolynomial(k-1, x)
  }


  # compute terms in the sum of relation (21) in
  terms <- exp((-qnorm(mycdf(0:max(N), a,m))^2)/2) *
    her(qnorm(mycdf(0:max(N), a,m)))

  terms[is.nan(terms)]=0

  # take the sum of all terms
  HC_k <- sum(terms) / (sqrt(2*pi) *  factorial(k))
  return(HC_k)
}



#----------Link coefficients with one dummy Regressor---------#
LinkCoefGP_Reg = function(a,m, N, nHC, mycdf){
  #====================================================================================#
  # PURPOSE    Compute the product  factorial(k)*g_{t1,k}*g_{t2,k} in relation (67) in
  #            https://arxiv.org/pdf/1811.00203.pdf when the regressor variable is only
  #            one dummy variable.
  #
  # INPUT
  #   a
  #   m        GP Marginal parameters
  #
  # Output
  #   l        An MX3 matrix of link coefficients. M is the default truncatuon of relation
  #            (67) and 3 is the number of different combinations of products
  #            g_{t1,k}*g_{t2,k} (since for each t I will either have 1 or 0 for the
  #            regressor variable).
  #
  # Authors    Stefanos Kechagias, James Livsey,  Vladas Pipiras
  # Date       July 2020
  # Version    3.6.3
  #====================================================================================#


  # compute Hermite coefficients from relation (21)
  g1 = HermCoefGenPois(a, unique(m)[1], N[1], nHC, mycdf)
  g2 = HermCoefGenPois(a, unique(m)[2], N[2], nHC, mycdf)
  HC = cbind(g1, g2)

  # Compute the products g1^2, g1*g2 and g2^2 of the HCs
  HCprod = matrix(NA, nHC, 3)
  for(i in 0:(length(unique(m))-1) ){
    for(j in i: (length(unique(m))-1) ){
      HCprod[,i+j+1] = HC[,i+1]*HC[,j+1]
    }
  }

  return(HCprod)
}



#---------Covariance matrix with one dummy Regressor---------#
CovarGenPois_Reg = function(n, a,m, AR, MA, N, nHC, mycdf){
  #====================================================================================#
  # PURPOSE    Compute the covariance matrix of a GenPois series that
  #            includes one dummy variable as a regressor. Here p depends on each
  #            observation. See relation (67) in https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   a
  #   m        GP MArginal parameters
  #   AR,MA    AR, MA parameter
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
    All.arma = apply(as.matrix(ARMAacf(ar = AR, ma=MA, lag.max = n)), 1,function(x)x^(1:nHC))
  }
  # Compute the link coefficients l_k = factorial(k)*g_{t1,k}*g_{t2,k}
  linkCoef = LinkCoefGP_Reg(a,m, N, nHC, mycdf)

  # keep track of which indices each unique HC is located at
  index = cbind(unique(m)[1]==m, unique(m)[2]==m)

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


#---------Covariance matrix---------#
CovarGenPois = function(n,a,m, AR, MA, N, nHC, mycdf){
  #====================================================================================#
  # PURPOSE     Compute the covariance matrix of a GenPois series.
  #
  # INPUT
  #   a,m       Marginal parameters
  #   AR,MA     ARMA parameters
  #   n         size of the matrix
  #   N         truncation for relation (21)
  #   nHC       number of Hermitte coefficents to be computed
  #
  # Output
  #   GAMMA     covariance matrix of count series
  #
  # Authors     Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date        July 2020
  # Version     3.6.3
  #====================================================================================#

  # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
  HC = HermCoefGenPois(a,m,N, nHC, mycdf)

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
GaussLogLikGP = function(theta, data, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE      Compute Gaussian log-likelihood fo Gen Pois series
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
  # Date         July 2020
  # Version      3.6.3
  #====================================================================================#


  # retrieve marginal cdf
  mycdf = switch(CountDist,
                 "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
                 "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
  )



  # retrieve parameters and sample size
  a   = theta[1]
  m  = theta[2]
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
  N <- which(round(pGpois(1:MaxCdf, a,m), 7) == 1)[1]
  # if(length(N)==0 |is.na(N) ){
  #   cat(sprintf("The max cdf value is %f and N=%f", max(round(pnbinom(1:MaxCdf, r,1-p), 7)),N))
  #   stop("Haven't reached upper limit for cdf")
  # }
  if(length(N)==0 |is.na(N) ){
    N =MaxCdf
  }

  #Select the mean value used to demean--sample or true?
  MeanValue = m

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA = CovarGenPois(n, a,m, AR, MA, N, nHC, mycdf)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}

#---------Guassian Likelihood function with a dummy Regressor---------#
GaussLogLikGP_Reg = function(theta, data, Regressor, ARMAorder, MaxCdf, nHC, CountDist){
  #====================================================================================#
  # PURPOSE      Compute Gaussian log-likelihood for Gen Pois series
  #
  #
  # NOTES        Here I use the following parametrization:
  #
  #              log(mu) = b0 + b1X_i, with
  #              E(Y/X) = m,
  #              The parameters that enter the likelihood are
  #              b0, b1, k, and the AR parameters.

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
  # Date         July 2020
  # Version      3.6.3
  #====================================================================================#


  # retrieve marginal cdf
  mycdf = switch(CountDist,
                 "Negative Binomial"   = function(x, ConstTheta, DynamTheta){ pnbinom (x, ConstTheta, 1-DynamTheta)},
                 "Generalized Poisson" = function(x, ConstTheta, DynamTheta){ pGpois  (x, ConstTheta, DynamTheta)}
  )


  # retrieve parameters and sample size
  nparms     = length(theta)
  nreg       = dim(Regressor)[2]-1
  nMargParms = nparms - sum(ARMAorder)
  beta       = theta[1:(nparms-sum(ARMAorder)-1)]
  a          = theta[nparms-sum(ARMAorder)]
  n          = length(data)

  if(ARMAorder[1]>0){
    AR = theta[(nMargParms+1):(nMargParms + ARMAorder[1])  ]
  }else{
    AR = NULL
  }

  if(ARMAorder[2]>0){
    MA = theta[ (nMargParms+ARMAorder[1]+1) : length(theta)]
  }else{
    MA = NULL
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
  m = exp(Regressor%*%beta)

  # Compute truncation of relation (21) in arxiv
  N = sapply(unique(m),function(x)which(mycdf(1:MaxCdf,a,x)>=1-1e-7)[1])-1
  N[is.na(N)] = MaxCdf

  #Select the mean value used to demean--sample or true?
  MeanValue = m

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  # GAMMA = CovarNegBin(n, r, p, AR, MA, N, nHC)
  GAMMA = CovarGenPois_Reg(n, a,m, AR,MA, N, nHC,mycdf)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}


#-----------Fit the gaussian likelihood----------#
FitGaussianLikGP = function(x0, X, LB, UB, ARMAorder, MaxCdf, nHC){
  #====================================================================================#
  # PURPOSE       Fit the Gaussian log-likelihood for Gen Pois series
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
                        fn        = GaussLogLikGP,
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


#---------wrapper to fit Guassian Likelihood function with a dummy Regressor---------#
FitGaussianLikGP_Reg = function(x0, X, Regressor, LB, UB, ARMAorder, MaxCdf, nHC, OptMethod, CountDist){
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


  optim.output <- optimx(par            = x0,
                         fn             = GaussLogLikGP_Reg,
                         data           = X,
                         Regressor      = Regressor,
                         ARMAorder      = ARMAorder,
                         MaxCdf         = MaxCdf,
                         CountDist      = CountDist,
                         nHC            = nHC,
                         lower          = LB,
                         upper          = UB,
                         hessian        = TRUE,
                         method         = OptMethod)

  # save estimates, loglik value and diagonal hessian
  ParmEst  = as.numeric(optim.output[1:nparms])
  loglik   = optim.output$value
  convcode = optim.output$convcode
  kkt1     = optim.output$kkt1
  kkt2     = optim.output$kkt2

  # compute hessian
  H = gHgen(par            = ParmEst,
            fn             = GaussLogLikGP_Reg,
            data           = X,
            Regressor      = Regressor,
            ARMAorder      = ARMAorder,
            CountDist      = CountDist,
            MaxCdf         = MaxCdf,
            nHC            = nHC
  )

  if(H$hessOK){
    se = sqrt(abs(diag(solve(H$Hn))))
  }else{
    se = rep(NA, nparms)
  }



  # Compute model selection criteria
  Criteria = ComputeCriteria(loglik, nparms, n)


  # get the names of the final output
  parmnames = colnames(optim.output)
  mynames = c(parmnames[1:nparms],paste("se", parmnames[1:nparms], sep="_"), "loglik", "AIC", "BIC","AICc", "status", "kkt1", "kkt2")



  All = matrix(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)
  colnames(All) = mynames
  return(All)

}


