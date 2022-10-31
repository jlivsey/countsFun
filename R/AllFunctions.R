#---------retrieve the model scheme
ModelScheme = function(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon, initialParam, EstMethod, maxit){

  error = 0
  errorMsg = NULL

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

  if(!is.null(initialParam) && length(initialParam)!=nparms) {
    error = 1
    errorMsg = "The length of the initial parameter doesn't match the model specifications."
  }

  # create the constraints
  if(CountDist =="Poisson"){
    if(nreg==0){
      LB = c(0.001, rep(-Inf, sum(ARMAorder)))
      UB = rep(Inf, sum(ARMAorder)+1)
    }else{
      LB = rep(-Inf, sum(ARMAorder)+nreg+1)
      UB = rep(Inf, sum(ARMAorder)+nreg+1)
    }
  }

  if(CountDist == "Negative Binomial"){
    if(nreg==0){
      LB = c(0.01, 0.01, rep(-Inf, sum(ARMAorder)))
      UB = c(Inf, 0.99,   rep( Inf, sum(ARMAorder)))
    }else{
      LB = c(rep(-Inf, nreg+1), 0.001, rep(-Inf, sum(ARMAorder)))
      UB = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAorder)))
    }
  }


  if(CountDist == "Generalized Poisson"){
    if(nreg==0){
      LB = c(0.001, 0.001, rep(-Inf, sum(ARMAorder)))
      UB = c(Inf, Inf,     rep( Inf, sum(ARMAorder)))
    }else{
      LB = c(rep(-Inf, nreg+1), 0.001, rep(-Inf, sum(ARMAorder)))
      UB = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAorder)))
    }
  }


  out = list(
    mycdf           = mycdf,
    mypdf           = mypdf,
    MargParmIndices = MargParmIndices,
    initialParam    = initialParam,
    nMargParms      = nMargParms,
    n               = n,
    nreg            = nreg,
    MaxCdf          = MaxCdf,
    nHC             = nHC,
    CountDist       = CountDist,
    ARMAorder       = ARMAorder,
    ParticleNumber  = ParticleNumber,
    epsilon         = epsilon,
    nparms          = nparms,
    UB              = UB,
    LB              = LB,
    error           = error,
    errorMsg        = errorMsg,
    EstMethod       = EstMethod,
    maxit           = maxit,
    data            = data,
    Regressor       = Regressor
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
CountCovariance = function(n, MargParms, ConstMargParm, DynamMargParm, AR, MA, N, nHC, mycdf, nreg, format){
  #====================================================================================#
  # PURPOSE
  # INPUT
  #   n               sample size
  #   MargParms       marginal parameters if no regressor
  #   ConstMargParm   contant marginal parameters if regressor is present
  #   DynamMargParm   marginal parameters depending on regressor
  #   AR              AR parameters
  #   MA              MA parameters
  #   N               cdf truncation
  #   nHC             number of Hermitte coefficients
  #   mycdf           cdf of marginal distribution
  #   nreg            number of regressors
  #   format          flag, 0 return first row, 1 return entire matrix
  # OUTPUT
  #   Covariance      matrix or first row depending on flag
  #
  # Authors           Stefanos Kechagias, James Livsey
  # Date              March 2020
  # Version           3.6.3
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
    if (format==1){
      G = toeplitz(gamma_x)
    }else{
      G = gamma_x
    }
  }
  return(G)
}


#---------Guassian Likelihood function ---------#
GaussianLogLik = function(theta, mod){
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
    m     = exp(mod$Regressor%*%beta)
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

  # check for causality: FIX ME, this introduces a jump when we do not have causal models, I guess it is a heuristic to avoid causal models
  # We need to discuss how this is treated.
  if( CheckStability(AR,MA) ){
    return(10^(8))
  }else{
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
    GAMMA =  CountCovariance(mod$n, MargParms, ConstMargParm, DynamMargParm, AR, MA, N, mod$nHC, mod$mycdf, mod$nreg, 1)

    # Compute the logdet and the quadratic part
    logLikComponents = EvalInvQuadForm(GAMMA, as.numeric(mod$data), MeanValue)

    # final loglikelihood value
    out = 0.5*logLikComponents[1] + 0.5*logLikComponents[2]

    # the following will match the above if you subtract N/2*log(2*pi) and don't multiply with 2
    # out = -2*dmvnorm(as.numeric(data), rep(lam, n), GAMMA, log = TRUE)
    return(out)
  }
}


#---------check causality and invertibility
CheckStability = function(AR,MA){
  if (is.null(AR) && is.null(MA)) return(0)

  # return 1 if model is not stable (causal and invertible)
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
FitGaussianLogLik = function(theta, mod){
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
  nparms    = length(mod$initialParam)
  n         = length(mod$data)

  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst   = matrix(0,nrow=1,ncol=nparms)
  se        = matrix(NA,nrow=1,ncol=nparms)
  loglik    = rep(0,1)
  convcode  = rep(0,1)
  kkt1      = rep(0,1)
  kkt2      = rep(0,1)


  optim.output <- optimx(
    par     = theta,
    fn      = GaussianLogLik,
    mod     = mod,
    lower   = mod$LB,
    upper   = mod$UB,
    method  = mod$OptMethod,
    hessian = TRUE)

  # save estimates, loglik value and diagonal hessian
  ParmEst   = as.numeric(optim.output[1:nparms])
  loglik    = optim.output$value
  convcode  = optim.output$convcode
  kkt1      = optim.output$kkt1
  kkt2      = optim.output$kkt2

  # compute sandwich standard errors
  se        = sand(ParmEst, data, Regressor, mod)

  # Compute model selection criteria
  Criteria  = ComputeCriteria(loglik, mod)


  # get the names of the final output
  parmnames = colnames(optim.output)
  mynames   = c(parmnames[1:nparms],paste("se", parmnames[1:nparms], sep="_"), "loglik", "AIC", "BIC","AICc", "status", "kkt1", "kkt2")

  # gather results
  All       = matrix(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)
  colnames(All) = mynames
  return(All)

}



#---------Compute AIC, BIC, AICc
ComputeCriteria = function(loglik, mod){
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

  if(mod$EstMethod!="GL"){
    l1 = -loglik
  }else{
    l1 = -loglik  - mod$n/2*log(2*pi)
  }

  AIC = 2*mod$nparms - 2*l1
  BIC = log(mod$n)*mod$nparms - 2*l1
  AICc = AIC + (2*mod$nparms^2 + 2*mod$nparms)/(mod$n-mod$nparms-1)

  AllCriteria = c(AIC, BIC, AICc)
}


#----------standard erros with sandwhich method following notation of Freedman (2006)
sand <- function(theta, data, Regressor, mod){

  # Calulate numerical Hessian
  h <- gHgen(fn         = GaussianLogLik,
             par       = theta,
             data      = data,           # additional arg for GaussLogLik
             Regressor = Regressor,      # additional arg for GaussLogLik
             mod       = mod)            # additional arg for GaussLogLik


  if(!(h$hessOK && det(h$Hn)>10^(-8))){
    SE.sand = rep(NA, mod$nparms)
  }else{
    gi <- matrix(NA, length(data), length(theta))

    for(k in 1:(length(data)-1)){
      gi[k, ] <- grad(func =  logf_i, x = theta, mod = mod, i = k)
    }
    gi <- gi[-length(data), ]

    # Calculate Cov matrix
    A <- h$Hn
    B <- t(gi) %*% gi
    V.sand  <- solve(-A) %*% B %*% solve(-A)
    SE.sand <- sqrt(diag(V.sand))

    # standard errors usual way using hessian
    #SE.hess <- sqrt(diag(solve(h)))
  }

  return(SE.sand)
}


#----------log of density at observation i, Freedman (2006)
logf_i <- function(theta, mod, i){

  # retrieve ARMA parameters
  AR = NULL
  if(mod$ARMAorder[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAorder[1])]

  MA = NULL
  if(mod$ARMAorder[2]>0) MA = theta[(mod$nMargParms+mod$ARMAorder[1]+1) : (mod$nMargParms + mod$ARMAorder[1] + mod$ARMAorder[2]) ]

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


  # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
  GAMMA =  CountCovariance(mod$n, MargParms, ConstMargParm, DynamMargParm, AR, MA, N, mod$nHC, mod$mycdf, mod$nreg, 0)


  # DL algorithm
  if(mod$nreg==0){
    DLout <- DLalg(data, GAMMA, MeanValue)
    ei <- DLout$e[i]
    vi <- DLout$v[i]
  }else{
    IAout <- Innalg(data, GAMMA)
    ei <- IAout$e[i]
    vi <- IAout$v[i]
  }



  # else{
  #   INAlgout = innovations.algorithm(gamma)
  #   Theta    = INAlgout$thetas
  #   vi <- INAlgout$vs[i]
  #
  #   ei <- sum(data - Theta*data)
  #
  #
  #
  #   Theta = ia$thetas
  #   # first stage of Innovations
  #   v0    = ia$vs[1]
  #   zhat = -Theta[[1]][1]*zprev
  #
  # }


  return(-(log(vi) + ei^2/vi)/2)
}






# PF likelihood with resampling for AR(p)
ParticleFilter_Res_AR = function(theta, mod){
  #--------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling
  #           to approximate the likelihood of the
  #           a specified count time series model with an underlying AR(p)
  #           dependence structure. A singloe dummy regression is added here.
  #
  # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
  #           details. A first version of the paer can be found at:
  #           https://arxiv.org/abs/1811.00203
  #           2. This function is very similar to LikSISGenDist_ARp but here
  #           I have a resampling step.
  #
  # INPUTS:
  #    theta:            parameter vector
  #    data:             data
  #    ParticleNumber:   number of particles to be used.
  #    CountDist:        count marginal distribution
  #    epsilon           resampling when ESS<epsilon*N
  #
  # OUTPUT:
  #    loglik:           approximate log-likelihood
  #
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    July  2020
  #--------------------------------------------------------------------------#

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # retrieve marginal parameters
  MargParms        = theta[mod$MargParmIndices]

  # retrieve regressor parameters
  if(mod$nreg>0){
    beta  = MargParms[1:(mod$nreg+1)]
    m     = exp(mod$Regressor%*%beta)
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
  if( CheckStability(AR,MA) ) return(10^(8))


  T1 = length(data)
  N = mod$ParticleNumber          # number of particles

  wgh = matrix(0,T1,N)        # to collect all particle weights

  # allocate memory for zprev
  ZprevAll = matrix(0,mod$ARMAorder[1],N)

  if(mod$nreg==0){
    # Compute integral limits
    a = rep( qnorm(mod$mycdf(data[1]-1,t(MargParms)),0,1), N)
    b = rep( qnorm(mod$mycdf(data[1],t(MargParms)),0,1), N)
  }else{
    a = rep( qnorm(mod$mycdf(data[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
    b = rep( qnorm(mod$mycdf(data[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
  }
  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

  # save the currrent normal variables
  ZprevAll[1,] = zprev

  # initial estimate of first AR coefficient as Gamma(1)/Gamma(0) and corresponding error
  # phit = TacvfAR(AR)[2]/TacvfAR(AR)[1] - this FitAR package is obsolete
  # FIX ME: check if the code below is correct
  phit = ARMAacf(ar = 0.2)[2]/ARMAacf(ar = 0.2)[1]
  rt = as.numeric(sqrt(1-phit^2))

  # particle filter weights
  wprev = rep(1,N)
  wgh[1,] = wprev
  nloglik = 0 # initialize likelihood

  # First p steps:
  if (mod$ARMAorder[1]>=2){
    for (t in 2:mod$ARMAorder[1]){

      # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
      if (t==2) {
        zhat = ZprevAll[1:(t-1),]*phit
      } else{
        zhat = colSums(ZprevAll[1:(t-1),]*phit)
      }

      # Recompute integral limits
      if(mod$nreg==0){
        a = (qnorm(mod$mycdf(data[t]-1,t(MargParms)),0,1) - zhat)/rt
        b = (qnorm(mod$mycdf(data[t],t(MargParms)),0,1) - zhat)/rt
      }else{
        a = (qnorm(mod$mycdf(data[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
        b = (qnorm(mod$mycdf(data[t],ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
      }

      # compute random errors from truncated normal
      err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

      # compute the new Z and add it to the previous ones
      znew = rbind(zhat + rt*err, ZprevAll[1:(t-1),])
      ZprevAll[1:t,] = znew

      # recompute weights
      wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
      wprev = wgh[t,]

      # use YW equation to compute estimates of phi and of the errors
      # FIX ME: Check if the code below is correct i nterms of the ARMAacf
      Gt    = toeplitz(ARMAacf(AR)[1:t])
      gt    = ARMAacf(AR)[2:(t+1)]
      phit  = as.numeric(solve(Gt) %*% gt)
      rt    =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/ARMAacf(AR)[1]))

    }
  }


  # From p to T1 I dont need to estimate phi anymore
  for (t in (mod$ARMAorder[1]+1):T1){

    # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
    if(mod$ARMAorder[1]>1){# colsums doesnt work for 1-dimensional matrix
      zhat = colSums(ZprevAll*AR)
    }else{
      zhat=ZprevAll*AR
    }

    # compute limits of truncated normal distribution
    if(mod$nreg==0){
      a = as.numeric(qnorm(mod$mycdf(mod$data[t]-1,MargParms),0,1) - zhat)/rt
      b = as.numeric(qnorm(mod$mycdf(mod$data[t],  MargParms),0,1) - zhat)/rt
    }else{
      a = as.numeric(qnorm(mod$mycdf(mod$data[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
      b = as.numeric(qnorm(mod$mycdf(mod$data[t],  ConstMargParm, DynamMargParm[t]),0,1) - zhat)/rt
    }

    # draw errors from truncated normal
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

    # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
    znew = zhat + rt*err

    # Resampling Step--here the function differs from LikSISGenDist_ARp

    # compute unnormalized weights
    wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)

    # break if I got NA weight
    if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
      # print(theta)
      # print(t)
      nloglik = 10^8
      # FIX ME: Do I need break here or return?
      break
      # I am fitting a poisson with initialParam= NULL in the Sales data and the initial estimated
      # mean form GLM, returns parameters that yield a, and b above that lead to Zprev = -Inf. I will comment
      # out the break and return, also change the value from 10^8 to 10^18.
      # return(nloglik)
    }

    # normalized weights
    wghn = wgh[t,]/sum(wgh[t,])

    old_state1 <- get_rand_state()

    # sample indices from multinomial distribution-see Step 4 of SISR in paper
    ESS = 1/sum(wghn^2)
    if(ESS<mod$epsilon*N){
      ind = rmultinom(1,N,wghn)
      # sample particles
      znew = rep(znew,ind)

      # use low variance resampling
      #znew = lowVarianceRS(znew, wghn, N)
    }
    set_rand_state(old_state1)


    # save particles
    if (mod$ARMAorder[1]>1){
      ZprevAll = rbind(znew, ZprevAll[1:(mod$ARMAorder[1]-1),])
    }else {
      ZprevAll[1,]=znew
    }
    # update likelihood
    nloglik = nloglik - log(mean(wgh[t,]))
  }

  # likelihood approximation
  if(mod$nreg==0){
    nloglik = nloglik - log(mod$mypdf(mod$data[1],MargParms))
  }else{
    # FIX ME: check this. what happens if the pdf is zero at a given point? Perhaps it is numerical zero
    # for now I ll ignore it but this is wrong!!
    # if (mod$mypdf(mod$data[1], ConstMargParm, DynamMargParm[1])<10^(-12)){
    #   nloglik = nloglik
    # }else{
    #   nloglik = nloglik - log(mod$mypdf(mod$data[1], ConstMargParm, DynamMargParm[1]))
    # }
    nloglik = nloglik - log(mod$mypdf(mod$data[1], ConstMargParm, DynamMargParm[1]))

  }

  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  nloglik = nloglik
  #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  return(nloglik)
}


# PF likelihood with resampling for MA(q)
ParticleFilter_Res_MA = function(theta, mod){
  #------------------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
  #           of the a specified count time series model with an underlying MA(1)
  #           dependence structure.
  #
  # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
  #           details. A first version of the paper can be found at:
  #           https://arxiv.org/abs/1811.00203
  #           2. This function is very similar to LikSISGenDist_ARp but here
  #           I have a resampling step.
  #
  # INPUTS:
  #    theta:            parameter vector
  #    data:             data
  #    ParticleNumber:   number of particles to be used.
  #    Regressor:        independent variable
  #    CountDist:        count marginal distribution
  #    epsilon           resampling when ESS<epsilon*N
  #
  # OUTPUT:
  #    loglik:           approximate log-likelihood
  #
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    July 2020
  #------------------------------------------------------------------------------------#

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # retrieve marginal parameters
  MargParms        = theta[mod$MargParmIndices]

  # retrieve regressor parameters
  if(mod$nreg>0){
    beta  = MargParms[1:(mod$nreg+1)]
    m     = exp(mod$Regressor%*%beta)
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


  T1 = length(mod$data)
  N = mod$ParticleNumber          # number of particles


  # allocate matrix to collect all particle weights
  wgh = matrix(0,length(mod$data),N)

  # Compute integral limits
  if(mod$nreg==0){
    a = rep( qnorm(mod$mycdf(mod$data[1]-1,t(MargParms)),0,1), N)
    b = rep( qnorm(mod$mycdf(mod$data[1],t(MargParms)),0,1), N)
  }else{
    a = rep( qnorm(mod$mycdf(mod$data[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
    b = rep( qnorm(mod$mycdf(mod$data[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
  }

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)


  # run innovations Algorithm for MA models that are not WN
  if(mod$ARMAorder[2]>0) Inn = matrix(0,N,mod$ARMAorder[2])   # I will save here the q many innovations (Z - Zhat) --see (5.3.9) BD book
  if (is.null(MA) && is.null(AR)){
    v0   = 1
    zhat = 0
  }else{
    # FIX ME: Check if the code below is correct in terms of the ARMAacf
    MA.acvf <- as.vector(ARMAacf(ma = MA, lag.max=T1))
    ia = innovations.algorithm(MA.acvf)
    Theta = ia$thetas
    # first stage of Innovations
    v0    = ia$vs[1]
    zhat = -Theta[[1]][1]*zprev
    Inn[,ARMAorder[2]] = zprev
  }

  # particle filter weights
  wprev   = rep(1,N)
  wgh[1,] = wprev
  nloglik = 0 # initialize likelihood

  for (t in 2:T1){

    # update innovations quantities if not White noise
    if (is.null(MA) && is.null(AR)){
      vt=1
    }else{
      vt0 = ia$vs[t]
      vt  = sqrt(vt0/v0)
    }

    # roll the old Innovations to earlier columns
    if(mod$ARMAorder[2]>1) Inn[,1:(mod$ARMAorder[2]-1)] = Inn[,2:(mod$ARMAorder[2])]

    # compute limits of truncated normal distribution
    if(mod$nreg==0){
      a = as.numeric(qnorm(mod$mycdf(mod$data[t]-1,MargParms),0,1) - zhat)/vt
      b = as.numeric(qnorm(mod$mycdf(mod$data[t],MargParms),0,1) -   zhat)/vt
    }else{
      a = as.numeric(qnorm(mod$mycdf(mod$data[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/vt
      b = as.numeric(qnorm(mod$mycdf(mod$data[t],ConstMargParm, DynamMargParm[t]),0,1) -   zhat)/vt
    }

    # draw errors from truncated normal
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

    # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
    znew = zhat + vt*err

    # compute new innovation
    Inn[,mod$ARMAorder[2]] = (znew-zhat)

    # compute unnormalized weights
    wgh[t,] = pnorm(b,0,1) - pnorm(a,0,1)

    # break if I got NA weight
    if (any(is.na(wgh[t,]))| sum(wgh[t,])<10^(-8) ){
      nloglik = 10^8
      break
    }

    # normalized weights
    wghn = wgh[t,]/sum(wgh[t,])


    # Resampling: sample indices from multinomial distribution-see Step 4 of SISR in paper
    ESS = 1/sum(wghn^2)
    old_state1 <- get_rand_state()
    if(ESS<mod$epsilon*N){
      ind = rmultinom(1,N,wghn)
      # sample particles
      znew = rep(znew,ind)

      # use low variance resampling
      #znew = lowVarianceRS(znew, wghn, N)
    }

    # update zhat--fix me can probably be vectorized
    if (is.null(MA) && is.null(AR)){
      zhat = 0
    }else{
      S = 0
      for(j in 1:min(t,mod$ARMAorder[2])){
        S = S-Theta[[t]][j]*Inn[,mod$ARMAorder[2]-j+1]
      }
      zhat = S
    }

    set_rand_state(old_state1)

    # update likelihood
    nloglik = nloglik - log(mean(wgh[t,]))
  }

  # likelihood approximation
  if(mod$nreg<1){
    nloglik = nloglik - log(mod$mypdf(mod$data[1],MargParms))
  }else{
    nloglik = nloglik - log(mod$mypdf(mod$data[1], ConstMargParm, DynamMargParm[1]))
  }

  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))


  if (nloglik==Inf | is.na(nloglik)){
    nloglik = 10^8
  }


  return(nloglik)
}


# PF likelihood with resampling
ParticleFilter_Res = function(theta, mod){
  #--------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling
  #           to approximate the likelihood of the
  #           a specified count time series model with an underlying AR(p)
  #           dependence structure or MA(q) structure.
  #
  # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
  #           details. A first version of the paer can be found at:
  #           https://arxiv.org/abs/1811.00203

  #
  # INPUTS:
  #    theta:            parameter vector
  #    data:             dependent variable
  #    Regressor:        independent variables
  #    ParticleNumber:   number of particles to be used in likelihood approximation
  #    CountDist:        count marginal distribution
  #    epsilon           resampling when ESS<epsilon*N
  #
  # OUTPUT:
  #    loglik:           approximate log-likelihood
  #
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
  # DATE:    July 2020
  #--------------------------------------------------------------------------#


  # Pure AR model
  if(mod$ARMAorder[1]>0 && mod$ARMAorder[2]==0) loglik = ParticleFilter_Res_AR(theta, mod)
  # Pure MA model or White noise
  if(mod$ARMAorder[1]==0&& mod$ARMAorder[2]>=0) loglik = ParticleFilter_Res_MA(theta, mod)
  return(loglik)
}


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


# innovations algorithm code
innovations.algorithm <- function(acvf,n.max=length(acvf)-1){
  # Found this onlinbe need to check it
  # http://faculty.washington.edu/dbp/s519/R-code/innovations-algorithm.R
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


# Nonstationary innovations algorithm
Innalg <- function(data, GAMMA){
  N <- length(data)
  x.hat <- numeric(N)
  v <- numeric(N)
  e <- numeric(N)
  theta <- matrix(0, N, N)

  x.hat[1] <- 0
  v[1] <- GAMMA[1, 1]
  e[1] <- data[1]

  for (n in 1:(N-1)){
    for (k in 0:(n-1)){
      a <- 0
      if (k > 0) {
        a <- sum(theta[k, 1:k] * theta[n, 1:k] * v[1:k])
      }

      theta[n, k+1] <- (1/v[k+1]) * (GAMMA[n+1, k+1] - a)
    }
    if(n<N){
      x.hat[n+1] <- sum(theta[n, 1:n] * (data[1:n] - x.hat[1:n]))
      v[n+1] <- GAMMA[n+1, n+1] - sum(theta[n, 1:n]^2 * v[1:n])
      e[n+1] <- data[n+1] - x.hat[n+1]
    }
  }

  return(list(x.hat = x.hat,
              theta = theta,
              v     = v    ,
              e     = e     ))
}


# Optimization wrapper to fit PF likelihood with resamplinbg
FitMultiplePF_Res = function(theta, mod){
  #====================================================================================#
  # PURPOSE       Fit the Particle Filter log-likelihood with resampling.
  #               This function maximizes the PF likelihood, nfit manys times for nparts
  #               many choices of particle numbers, thus yielding a total of nfit*nparts
  #               many estimates.
  #
  # INPUT
  #   theta       initial parameters
  #   data        count series
  #   CountDist   prescribed count distribution
  #   Particles   vector with different choices for number of particles
  #   ARMAorder   order of the udnerlying ARMA model
  #   epsilon     resampling when ESS<epsilon*N
  #   LB          parameter lower bounds
  #   UB          parameter upper bounds
  #   OptMethod
  #
  #
  # OUTPUT
  #   All         parameter estimates, standard errors, likelihood value, AIC, etc
  #
  # Authors       Stefanos Kechagias, James Livsey, Vladas Pipiras
  # Date          April 2020
  # Version       3.6.3
  #====================================================================================#

  # retrieve parameter, sample size etc
  nparts   = length(mod$ParticleNumber)
  nparms   = length(theta)
  nfit     = 1
  n        = length(mod$data)

  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst  = matrix(0,nrow=nfit*nparts,ncol=nparms)
  se       = matrix(NA,nrow=nfit*nparts,ncol=nparms)
  loglik   = rep(0,nfit*nparts)
  convcode = rep(0,nfit*nparts)
  kkt1     = rep(0,nfit*nparts)
  kkt2     = rep(0,nfit*nparts)


  # Each realization will be fitted nfit*nparts many times
  for (j in 1:nfit){
    set.seed(j)
    # for each fit repeat for different number of particles
    for (k in 1:nparts){
      # FIX ME: I need to somehow update this in mod. (number of particles to be used). I t now works only for 1 choice of ParticleNumber
      ParticleNumber = mod$ParticleNumber[k]

      # run optimization for our model --no ARMA model allowed
      optim.output <- optimx(
        par = theta,
        fn  = ParticleFilter_Res,
        lower = mod$LB,
        upper = mod$UB,
        hessian = TRUE,
        method  = mod$OptMethod,
        mod   = mod)



      # save estimates, loglik value and diagonal hessian
      ParmEst[nfit*(k-1)+j,]  = as.numeric(optim.output[1:nparms])
      loglik[nfit*(k-1) +j]   = optim.output$value
      convcode[nfit*(k-1) +j] = optim.output$convcode
      kkt1[nfit*(k-1) +j]     = optim.output$kkt1
      kkt2[nfit*(k-1) +j]     = optim.output$kkt2


      # compute Hessian
      H = gHgen(par            = ParmEst[nfit*(k-1)+j,],
                fn             = ParticleFilter_Res,
                mod            = mod)

      # if I get all na for one row and one col of H remove it
      # H$Hn[rowSums(is.na(H$Hn)) != ncol(H$Hn), colSums(is.na(H$Hn)) != nrow(H$Hn)]

      if (!any(is.na(rowSums(H$Hn)))){
        # save standard errors from Hessian
        if(H$hessOK && det(H$Hn)>10^(-8)){
          se[nfit*(k-1)+j,]   = sqrt(abs(diag(solve(H$Hn))))
        }else{
          se[nfit*(k-1)+j,] = rep(NA, nparms)
        }
      }else{
        # remove the NA rows and columns from H
        Hnew = H$Hn[rowSums(is.na(H$Hn)) != ncol(H$Hn), colSums(is.na(H$Hn)) != nrow(H$Hn)]

        # find which rows are missing and which are not
        NAIndex = which(colSums(is.na(H$Hn))==nparms)
        NonNAIndex = which(colSums(is.na(H$Hn))==1)

        #repeat the previous ifelse for the reduced H matrix
        if(det(Hnew)>10^(-8)){
          se[nfit*(k-1)+j,NonNAIndex]   = sqrt(abs(diag(solve(Hnew))))
        }else{
          se[nfit*(k-1)+j,NAIndex] = rep(NA, length(NAindex))
        }

      }


    }
  }

  # Compute model selection criteria (assuming one fit)
  Criteria = ComputeCriteria(loglik, mod)


  # get the names of the final output
  parmnames = colnames(optim.output)
  mynames = c(parmnames[1:nparms],paste("se", parmnames[1:nparms], sep="_"), "loglik", "AIC", "BIC","AICc", "status", "kkt1", "kkt2")


  All = matrix(c(ParmEst, se, loglik, Criteria, convcode, kkt1, kkt2),nrow=1)
  colnames(All) = mynames

  return(All)
}




# compute initial estimates
InitialEstimates = function(mod){

  est  = rep(NA, mod$nMargParms+sum(mod$ARMAorder))

  #-----Poisson case
  if(mod$CountDist=="Poisson"){
    if(mod$nreg==0){
      est[1] = mean(mod$data)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      # GLM for the mean that depends on time
      # CHECK ME: If I fit a Poisson AR(3) in the the data example of the JASA paper, but the code below doesn't specify poisson family (it would pick up the default distribution that glm function has) then there will be a numerical error in the likelihood. Check it!
      glmPoisson            = glm(mod$data~mod$Regressor[,2:(mod$nreg+1)], family = "poisson")
      est[1:mod$nMargParms] = as.numeric(glmPoisson[1]$coef)

      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }

  #-----Neg Binomial case
  if(mod$CountDist=="Negative Binomial"){
    if(mod$nreg==0){
      xbar = mean(mod$data)
      sSquare = var(mod$data)

      # Method of Moments for negBin
      rEst = xbar^2/(sSquare - xbar)
      pEst = 1 - xbar/sSquare
      est[1:2] = c(rEst, pEst)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      # GLM for the mean that depends on time
      glmNegBin                 = glm.nb(mod$data~mod$Regressor[,2:(mod$nreg+1)])
      est[1:(mod$nMargParms-1)] = as.numeric(glmNegBin[1]$coef)
      # Mom on constant variance
      est[mod$nMargParms]       = NegBinMoM(mod$data,glmNegBin$fitted.values)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }



  if(mod$CountDist=="Mixed Poisson"){
    if(mod$nreg==0){
      # pmle for marginal parameters
      MixPois_PMLE <- pmle.pois(x,2)

      pEst  = MixPois_PMLE[[1]][1]
      l1Est = MixPois_PMLE[[2]][1]
      l2Est = MixPois_PMLE[[2]][2]


      # correct estimates if they are outside the feasible region
      if(pEst<LB[1]){pEst = 1.1*mod$LB[1]}
      if(pEst>UB[1]){pEst = 0.9*mod$UB[1]}

      if(l1Est<LB[2]){l1Est = 1.1*mod$LB[2]}
      if(l2Est<LB[3]){l2Est = 1.1*mod$LB[3]}

      est[1:3] = c(l1Est, l1Est, pEst)

      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(mod$data~mod$Regressor)[1]$coef)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }



  #-----Generalized Poisson case
  if(mod$CountDist=="Generalized Poisson"){
    if(mod$nreg==0){
      xbar = mean(mod$data)
      sSquare = var(mod$data)

      # Method of Moments for negBin
      rEst = xbar^2/(sSquare - xbar)
      pEst = 1 - xbar/sSquare
      est[1:2] = c(rEst, pEst)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$nMargParms+mod$ARMAorder[1]):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(mod$data~mod$Regressor)[1]$coef)
      if(mod$ARMAorder[1]) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$ARMAorder[1])] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$phi
      if(mod$ARMAorder[2]) est[(1+mod$ARMAorder[1]+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAorder))] = itsmr::arma(mod$data,mod$ARMAorder[1],mod$ARMAorder[2])$theta
    }
  }





  return(est)
}

NegBinMoM = function(data, GLMMeanEst){
  # the GLMMeanEst is the GLM estimate of the standard log-link
  # th formula below is standard MoM for the over dispersion param in NegBin2 parametrization
  PhiMomEst = sum(GLMMeanEst^2)/(sum((data-GLMMeanEst)^2-GLMMeanEst))
  return(PhiMomEst)
}







# # check validity of input arguments
# CheckInputSpecs = function(data, Regressor, CountDist, EstMethod, ARMAorder,
#                            nHC, MaxCdf, ParticleNumber, epsilon, initialParam,OptMethod ){
#
#   rc = 0
#
#   # Truncation of cdf
#   if (MaxCdf<0) rc = 1
#
#   # check distributions
#   if ( !(EstMethod %in%  c("PFR", "GL", "IYW")))  rc = 2
#
#   # check Particle number
#   if (EstMethod=="PFR" && (epsilon > 1 || epsilon<0)) rc = 3
#
#   # check Particle number
#   if (EstMethod=="PFR" && ParticleNumber<1) rc = 4
#
#   # check distributions
#   if ( !(CountDist %in%  c("Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial" ))) rc = 5
#
#   # check ARMAorder
#   if (prod(ARMAorder)<0 || length(ARMAorder)!= 2) rc = 6
#
#   # Mixed ARMA model
#   if (mod$ARMAorder[1]>0 && mod$ARMAorder[2]>0) rc = 7
#
#   # check data
#   if (is.null(data) ||  length(data)<3) rc = 8
#
#   errors = list(
#     e1 = 'Please specify a nonnegative MaxCdf.',
#     e2 = 'The argument EstMethod must take one of the following values: "GL", IYW","PFR".',
#     e3 = 'Please select a value between 0 and 1 for epsilon.',
#     e4 = 'Please select a nonegative value for the argument ParticleNumber.',
#     e5 = 'The argument CountDist must take one of the following values:
#          "Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial".',
#     e6 = 'The argument ARMAorder must have length 2 and can not take negative values.',
#     e7 = 'Please specify a pure AR or a pure MA model. ARMA(p,q) models with p>0 and q>0 have not yet been implemented.',
#     e8 = 'Data must be non null with sample size greater than 3.'
#   )
#
#
#   out = list(
#     rc  = rc,
#     e   = errors[[rc]]
#   )
#   return(out)
#
# }
#



