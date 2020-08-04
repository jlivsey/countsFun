#---------retrieve the model scheme
ModelScheme = function(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon ){



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
                   "Generalized Poisson" = function(x,theta) { pGpois  (x, theta[1], theta[2])},
                   "Binomial"            = pbinom
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "Generalized Poisson" = function(x,theta) { dGpois  (x, theta[1], theta[2])},
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
  nMargParms = length(MargParmIndices)
  nparms     = nMargParms + sum(ARMAorder)


  out = list(mycdf       = mycdf,
       mypdf           = mypdf,
       MargParmIndices = MargParmIndices,
       nMargParms      = nMargParms,
       n               = n,
       nreg            = nreg,
       MaxCdf          = MaxCdf,
       nHC             = nHC,
       CountDist       = CountDist,
       ARMAorder       = ARMAorder,
       ParticleNumber  = ParticleNumber,
       epsilon         = epsilon
       )
  return(out)

}

#---------Hermitte Coefficients for all k ---------#
HermCoef <- function(MargParms, ConstMargParm, DynamMargParm, N, nHC, mycdf, nreg){
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
    HC[i] <- HermCoef_k(MargParms, ConstMargParm, DynamMargParm , k = i, N, mycdf, nreg)
  }

  return(HC)

}


#---------Hermitte Coefficients for one k ---------#
HermCoef_k <- function(MargParms, ConstMargParm, DynamMargParm, k, N, mycdf, nreg){
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


#----------Link coefficients with one dummy Regressor ---------#
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
  g1 = HermCoef(MargParms = NULL, ConstMargParm, unique(DynamMargParm)[1],N[1], nHC, mycdf, nreg)
  g2 = HermCoef(MargParms = NULL, ConstMargParm, unique(DynamMargParm)[2],N[2], nHC, mycdf, nreg)
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


#---------Covariance matrix ---------#
CountCovariance = function(n, MargParms, ConstMargParm, DynamMargParm, AR, MA, N, nHC, mycdf, nreg){
  #====================================================================================#
  # PURPOSE
  # Authors    Stefanos Kechagias, James Livsey
  # Date       March 2020
  # Version    3.6.3
  #====================================================================================#

  if(nreg>0){
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
  }else{
    # Hermite coeficients--relation (21) in https://arxiv.org/pdf/1811.00203.pdf
    HC = HermCoef(MargParms, ConstMargParm=NULL, DynamMargParm=NULL, N, nHC, mycdf, nreg)

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
    G = toeplitz(gamma_x)
  }
  return(G)
}


#---------Guassian Likelihood function ---------#
GaussianLogLik = function(theta, data, Regressor, mod){
  #====================================================================================#
  # PURPOSE      Compute Gaussian log-likelihood for Count series. This will lead to a
  #              pseudo MLE.
  #
  #
  # Output
  #   loglik     Gaussian log-likelihood
  #
  # Authors      Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date         April 2020
  # Version      3.6.3
  #====================================================================================#

  # retrieve marginal parameters
  MargParms        = theta[mod$MargParmIndices]

  # retrieve regressor parameters
  if(mod$nreg>0){
    beta  = MargParms[1:(mod$nreg+1)]
    m     = exp(Regressor%*%beta)
  }

  # retrieve GLM type  parameters
  if(mod$CountDist == "Negative Binomial" && mod$nreg>0){
    ConstMargParm  = 1/MargParms[mod$nreg+2]
    DynamMargParm  = MargParms[mod$nreg+2]*m/(1+MargParms[mod$nreg+2]*m)
  }

  if(mod$CountDist == "Generalized Poisson" && mod$nreg>0){
    ConstMargParm  = MargParms[mod$nreg+2]
    DynamMargParm  = m
  }

  if(mod$CountDist == "Poisson" && mod$nreg>0){
    ConstMargParm  = NULL
    DynamMargParm  = m
  }

  # retrieve ARMA parameters
  AR = NULL
  if(mod$ARMAorder[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAorder[1])]

  MA = NULL
  if(mod$ARMAorder[2]>0) MA = theta[(mod$nMargParms+mod$ARMAorder[1]+1) : (mod$nMargParms + mod$ARMAorder[1] + mod$ARMAorder[2]) ]

  # check for causality
  if( CheckStability(AR,MA) ) return(10^(-6))


  # retrieve mean
  if(mod$nreg>0){
    MeanValue = m
  }else{
    MeanValue = switch(mod$CountDist,
                       "Poisson"             = MargParms[1],
                       "Negative Binomial"   = MargParms[1]*MargParms[2]/(1-MargParms[2]),
                       "Generalized Poisson" = MargParms[2])
  }

  # Compute truncation of relation (21)
  if(mod$nreg>0){
    N = sapply(unique(DynamMargParm),function(x)which(mod$mycdf(1:mod$MaxCdf, ConstMargParm, x)>=1-1e-7)[1])-1
    N[is.na(N)] = mod$MaxCdf
  }else{
    N <- which(round(mod$mycdf(1:mod$MaxCdf, MargParms), 7) == 1)[1]
    if(length(N)==0 |is.na(N) ){
      N =mod$MaxCdf
    }
  }

  # Compute the covariance matrix--relation (56) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA =  CountCovariance(mod$n, MargParms, ConstMargParm, DynamMargParm, AR, MA, N, mod$nHC, mod$mycdf, mod$nreg)

  # Compute the logdet and the quadratic part
  logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(data), MeanValue)

  # final loglikelihood value
  out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

  # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
  # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
  return(out)
}


#---------check causality and invertibility
CheckStability = function(AR,MA){
  if (is.null(AR) && is.null(MA)) return(0)
  # return 1 if model is not satble (causal and invertible)
  if(!is.null(AR) && is.null(MA)){
    rc = ifelse(any(abs( polyroot(c(1, -AR))  ) < 1), 1,0)
    }

  if(!is.null(MA) && is.null(AR)){
    rc = ifelse(any(abs( polyroot(c(1, -MA))  ) < 1),1,0)
  }

  if(!is.null(MA) && !is.null(AR)){
    rc = ifelse(checkPoly(AR,MA)[1]!="Causal" && check(poly)[2]!="Invertible", 1,0)
  }

  return(rc)
}


# compute arbitrary polynomial
evalPolynomial_scalarX <- function(coefs, x){
  #=======================================================================#
  # PURPOSE    compute arbitrary polynomial given coefficients and input
  #
  # INPUT
  #   coefs    vector of coefficients
  #   x        input
  #
  # Output
  #   value   returns f(x) = coefs[1] * x^0 +
  #                          coefs[2] * x^1 +
  #                                     ... +
  #                          coefs[n] * x^(n-1)
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #=======================================================================#
  if(length(coefs) == 1){
    return(coefs) # handle scalar case
  }
  n <- length(coefs) - 1
  out <- sum(coefs * x^(0:n))
  return(out)
}


# Evaluate k^th Hermite Polynomial at scalar x
evalHermPolynomial <- function(k, x){
  #=======================================================================#
  # PURPOSE    Evaluate k^th Hermite Polynomial at scaler input.
  #           See relation (7) https://arxiv.org/pdf/1811.00203.pdf
  #
  # INPUT
  #   k        index of Hermite coefficient
  #   x        input to Hermite Polynomial
  #
  #
  # Output
  #   value   returns H_k(x). See relation (7)
  #           in https://arxiv.org/pdf/1811.00203.pdf
  #
  # Authors    Stefanos Kechagias, James Livsey
  # Date       January 2020
  # Version    3.6.1
  #=======================================================================#
  if(k < 0){
    return(0)
  }
  coefs <- HermPolyCoefs[[k+1]]
  out <- evalPolynomial(coefs, x)
  return(out)
}


# Evaluate Hermite Polynomial at a vector x
evalPolynomial <- Vectorize(evalPolynomial_scalarX, vectorize.args = "x")


# count acvf for each h
CountACVF_h = function(h, myacf, g){
  #=======================================================================#
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
  #=======================================================================#

  k = length(g) #check me
  gamma_x = sum(g^2 *  factorial(1:k) * (myacf[h+1])^(1:k))
  return(gamma_x)

}


# vectorized count acvf
CountACVF <- Vectorize(CountACVF_h, vectorize.args = "h")


# fnc to evaluate Gaussial lik components
EvalInvQuadForm = function(A, v, DataMean ){
  #=======================================================================#
  # Evaluate quadrative form v`*inv(A)*v where
  # A is symmetric positive definite
  # Want   QuadForm = v` * inv(A) * v = v` * w  where A*w = v
  # ==> (U`U)*w = v             Cholesky decomp of A
  # ==>  First solve U` z = v   for z,
  # then solve   U w = z   for w */
  #=======================================================================#

  U = chol(A)
  z  =  forwardsolve(t(U), v - DataMean)
  w  = backsolve(U,z)
  QuadForm = t(v-DataMean)%*%w

  logdetA = 2*sum(log(diag(U)))

  logLikComponents = c(logdetA, QuadForm)
  return(logLikComponents)
}


#---------wrapper to fit Guassian Likelihood function with a dummy Regressor---------#
FitGaussianLogLik = function(theta, xt, Regressor, mod, LB, UB, OptMethod){
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
  nparms  = length(theta)
  n = length(xt)
  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst = matrix(0,nrow=1,ncol=nparms)
  se =  matrix(NA,nrow=1,ncol=nparms)
  loglik = rep(0,1)
  convcode = rep(0,1)
  kkt1 = rep(0,1)
  kkt2 = rep(0,1)


  optim.output <- optimx(par            = theta,
                         fn             = GaussianLogLik,
                         data           = xt,
                         Regressor      = Regressor,
                         mod            = mod,
                         lower          = LB,
                         upper          = UB,
                         method         = OptMethod,
                         hessian        = TRUE)

  # save estimates, loglik value and diagonal hessian
  ParmEst  = as.numeric(optim.output[1:nparms])
  loglik   = optim.output$value
  convcode = optim.output$convcode
  kkt1     = optim.output$kkt1
  kkt2     = optim.output$kkt2

  # compute hessian
  H = gHgen(par            = ParmEst,
            fn             = GaussianLogLik,
            data           = xt,
            Regressor      = Regressor,
            mod            = mod)

  # get standard errors
  if(H$hessOK){
    se = sqrt(abs(diag(solve(H$Hn))))
  }else{
    se = rep(NA, nparms)
  }


  # Compute model selection criteria
  Criteria = ComputeCriteria(loglik, nparms, n, mod$ParticleNumber)


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
