#---------retrieve the model scheme
ModelScheme = function(DependentVar, Regressor=NULL, EstMethod="PFR", ARMAModel=c(0,0), CountDist="Poisson",
                       ParticleNumber = 5, epsilon = 0.5, initialParam=NULL, TrueParam=NULL, Task="Optimization", SampleSize=NULL,
                       OptMethod="bobyqa", OutputType="list", ParamScheme=1, maxdiff=10^(-8)){

  # retrieve sample size
  n = ifelse(!is.null(DependentVar), length(DependentVar), SampleSize)

  # Distribution list
  if( !(CountDist %in% c("Poisson", "Negative Binomial", "Generalized Poisson", "Mixed Poisson", "ZIP", "Binomial")))
    stop("The specified distribution in not supported.")
  # Task
  if( !(Task %in% c("Evaluation", "Optimization", "Simulation")))
    stop("The specified distribution in not supported.")

  # Estimation Method
  if( !(EstMethod %in% c("PFR")))
    stop("The specified estimation method in not supported.")

  # check that the ARMA order has dimension 2, the ARMA orders are integers
  if(length(ARMAModel)!=2 | !prod(ARMAModel %% 1 == c(0,0)) )
    stop("The specified ARMA model is not supported.")


  # number of regressors assuming there is an intercept
  nreg = ifelse(is.null(Regressor), 0,dim(Regressor)[2]-1)

  # retrieve indices of marginal distribution parameters-the regressor is assumed to have an intercept
  MargParmIndices = switch(CountDist,
                           "Poisson"             = 1:(1+nreg),
                           "Negative Binomial"   = 1:(2+nreg),
                           "Generalized Poisson" = 1:(2+nreg),
                           "Binomial"            = 1:(2+nreg),
                           "Mixed Poisson"       = 1:(3+nreg*2),
                           "ZIP"                 = 1:(2+nreg*2)
  )

  # retrieve marginal distribution parameters
  nMargParms = length(MargParmIndices)
  nparms     = nMargParms + sum(ARMAModel)
  nAR        = ARMAModel[1]
  nMA        = ARMAModel[2]

  # check if initial param length is wrong
  if(!is.null(initialParam) & length(initialParam)!=nparms)
    stop("The specified initial parameter has wrong length.")

  # parse all information in the case without Regressors or in the case with Regressors
  if(nreg<1){
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"             = ppois,
                   "Negative Binomial"   = function(x, theta){ pnbinom (x, theta[1], 1-theta[2])},
                   "Generalized Poisson" = function(x, theta) { pGpois  (x, theta[1], theta[2])},
                   "Binomial"            = pbinom,
                   "Mixed Poisson"       = function(x, theta){ pmixpois(x, theta[1], theta[2], theta[3])},
                   "ZIP"                 = function(x, theta){ pzipois(x, theta[1], theta[2])}
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"             = dpois,
                   "Negative Binomial"   = function(x, theta){ dnbinom (x, theta[1], 1-theta[2]) },
                   "Generalized Poisson" = function(x, theta){ dGpois  (x, theta[1], theta[2])},
                   "Binomial"            = dbinom,
                   "Mixed Poisson"       = function(x, theta){ dmixpois(x, theta[1], theta[2], theta[3])},
                   "ZIP"                 = function(x, theta){ dzipois(x, theta[1], theta[2]) }
    )

    # retrieve marginal inverse cdf
    myinvcdf = switch(CountDist,
                      "Poisson"             = qpois,
                      "Negative Binomial"   = function(x, theta){ qnbinom (x, theta[1], 1-theta[2]) },
                      "Generalized Poisson" = function(x, theta){ qGpois  (x, theta[1], theta[2])},
                      "Binomial"            = qbinom,
                      "Mixed Poisson"       = function(x, theta){ qmixpois(x, theta[1], theta[2], theta[3])},
                      "ZIP"                 = function(x, theta){ qzipois(x, theta[1], theta[2]) }
    )

    # lower bound constraints
    LB = switch(CountDist,
                "Poisson"              = c(0.01,              rep(-Inf, sum(ARMAModel))),
                "Negative Binomial"    = c(0.01,  0.01,       rep(-Inf, sum(ARMAModel))),
                "Generalized Poisson"  = c(0.01, 0.01,        rep(-Inf, sum(ARMAModel))),
                "Binomial"             = c(0.01,  0.01,       rep(-Inf, sum(ARMAModel))),
                "Mixed Poisson"        = c(0.01, 0.01, 0.01,  rep(-Inf, sum(ARMAModel))),
                "ZIP"                  = c(0.01, 0.01,        rep(-Inf, sum(ARMAModel)))
    )
    # upper bound constraints
    UB = switch(CountDist,
                "Poisson"              = c(Inf,            rep( Inf, sum(ARMAModel))),
                "Negative Binomial"    = c(Inf, 0.99,      rep( Inf, sum(ARMAModel))),
                "Generalized Poisson"  = c(Inf, Inf,       rep( Inf, sum(ARMAModel))),
                "Binomial"             = c(Inf, 0.99,      rep( Inf, sum(ARMAModel))),
                "Mixed Poisson"        = c(Inf, Inf, 0.99, rep( Inf, sum(ARMAModel))),
                "ZIP"                  = c(Inf, 0.99,      rep( Inf, sum(ARMAModel)))
    )
    # names of marginal parameters
    MargParmsNames = switch(CountDist,
                            "Poisson"              = c("lambda"),
                            "Negative Binomial"    = c("r","p"),
                            "Mixed Poisson"        = c("lambda_1", "lambda_2", "p"),
                            "Generalized Poisson"  = c("lambda", "a"),
                            "Binomial"             = c("n", "p"),
                            "ZIP"                  = c("lambda", "p")
    )
  }else{
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"              = function(x, ConstMargParm, DynamMargParm){ ppois   (x, DynamMargParm)},
                   "Negative Binomial"    = function(x, ConstMargParm, DynamMargParm){ pnbinom (x, ConstMargParm, 1-DynamMargParm)},
                   "Generalized Poisson"  = function(x, ConstMargParm, DynamMargParm){ pGpois  (x, ConstMargParm, DynamMargParm)},
                   "Binomial"             = function(x, ConstMargParm, DynamMargParm){ pbinom  (x, ConstMargParm, DynamMargParm)},
                   "Mixed Poisson"        = function(x, ConstMargParm, DynamMargParm){ pmixpois(x, DynamMargParm, ConstMargParm)},
                   "ZIP"                  = function(x, ConstMargParm, DynamMargParm){ pzipois (x, DynamMargParm[1], DynamMargParm[2]) }
    )
    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"              = function(x, ConstMargParm, DynamMargParm){ dpois   (x, DynamMargParm)},
                   "Negative Binomial"    = function(x, ConstMargParm, DynamMargParm){ dnbinom (x, ConstMargParm, 1-DynamMargParm)},
                   "Generalized Poisson"  = function(x, ConstMargParm, DynamMargParm){ dGpois  (x, ConstMargParm, DynamMargParm)},
                   "Binomial"             = function(x, ConstMargParm, DynamMargParm){ dbinom  (x, ConstMargParm, DynamMargParm)},
                   "Mixed Poisson"        = function(x, ConstMargParm, DynamMargParm){ dmixpois(x, DynamMargParm, ConstMargParm)},
                   "ZIP"                  = function(x, ConstMargParm, DynamMargParm){ dzipois (x, DynamMargParm[1], DynamMargParm[2]) }
    )
    # retrieve marginal inverse cdf
    myinvcdf = switch(CountDist,
                      "Poisson"             = function(x, ConstMargParm, DynamMargParm){ qpois   (x, DynamMargParm)},
                      "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ qnbinom (x, ConstMargParm, 1-DynamMargParm)},
                      "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ qGpois  (x, ConstMargParm, DynamMargParm)},
                      "Binomial"            = function(x, ConstMargParm, DynamMargParm){ qbinom  (x, ConstMargParm, DynamMargParm)},
                      "Mixed Poisson"       = function(x, ConstMargParm, DynamMargParm){ qmixpois(x, DynamMargParm, ConstMargParm)},
                      "ZIP"                 = function(x, ConstMargParm, DynamMargParm){ qzipois (x, DynamMargParm[1], DynamMargParm[2]) }
    )
    # lower bound contraints
    LB = switch(CountDist,
                "Poisson"             = rep(-Inf, sum(ARMAModel)+nreg+1),
                "Negative Binomial"   = c(rep(-Inf, nreg+1), 0.001, rep(-Inf, sum(ARMAModel))),
                "Generalized Poisson" = c(rep(-Inf, nreg+1), 0.001, rep(-Inf, sum(ARMAModel))),
                "Binomial"            = c(rep(-Inf, nreg+1), 0.01, rep(-Inf, sum(ARMAModel))),
                "Mixed Poisson"       = c(rep(-Inf, 2*nreg+2), 0.001, rep(-Inf, sum(ARMAModel))),
                "ZIP"                 = c(rep(-Inf, 2*nreg+2), rep(-Inf, sum(ARMAModel)))
    )
    # upper bound constraints
    UB = switch(CountDist,
                "Poisson"             = rep(Inf, sum(ARMAModel)+nreg+1),
                "Negative Binomial"   = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAModel))),
                "Generalized Poisson" = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAModel))),
                "Binomial"            = c(rep(Inf, nreg+1), Inf, rep(Inf, sum(ARMAModel))),
                "Mixed Poisson"       = c(rep(Inf, 2*nreg+2), 0.99, rep(Inf, sum(ARMAModel))),
                "ZIP"                 = c(rep(Inf, 2*nreg+2), rep(Inf, sum(ARMAModel)))
    )
    # retrieve names of marginal parameters
    MargParmsNames = switch(CountDist,
                            "Poisson"             = paste(rep("b_",nreg),0:nreg,sep=""),
                            "Negative Binomial"   = c(paste(rep("b_",nreg),0:nreg,sep=""), "k"),
                            "Mixed Poisson"       = c(paste(rep("b_1",nreg),0:nreg,sep=""),paste(rep("b_2",nreg),0:nreg,sep=""), "p"),
                            "Generalized Poisson" = c(paste(rep("b_",nreg),0:nreg,sep=""), "a"),
                            "Binomial"            = c(paste(rep("b_",nreg),0:nreg,sep=""), "n"),
                            "ZIP"                 = c(paste(rep("b_",nreg),0:nreg, sep=""), paste(rep("c_",nreg),0:nreg, sep=""))
    )
  }


  # check whether the provided initial estimates make sense
  if (!is.null(initialParam) & (prod(initialParam<=LB) | prod(initialParam>=UB)))
    stop("The specified initial parameter is outside thew feasible region.")




  # create names of the ARMA parameters
  if(nAR>0) ARNames = paste("AR_",1:ARMAModel[1], sep="")
  if(nMA>0) MANames = paste("MA_",1:ARMAModel[2], sep="")

  # put all the names together
  if(nAR>0 && nMA<1) parmnames = c(MargParmsNames, ARNames)
  if(nAR<1 && nMA>0) parmnames = c(MargParmsNames, MANames)
  if(nAR>0 && nMA>0) parmnames = c(MargParmsNames, ARNames, MANames)

  # add the parmnames on theta fix me: does this affect performance?
  if(!is.null(initialParam)) names(initialParam) = parmnames

  # value I wil set the loglik when things go bad (e.g. non invetible ARMA)
  loglik_BadValue1 = 10^8

  # value I wil set the loglik when things go bad (e.g. non invetible ARMA)
  loglik_BadValue2 = 10^9

  out = list(
    mycdf           = mycdf,
    mypdf           = mypdf,
    myinvcdf        = myinvcdf,
    MargParmIndices = MargParmIndices,
    initialParam    = initialParam,
    TrueParam       = TrueParam,
    parmnames       = parmnames,
    nMargParms      = nMargParms,
    nAR             = nAR,
    nMA             = nMA,
    n               = n,
    nreg            = nreg,
    CountDist       = CountDist,
    ARMAModel       = ARMAModel,
    ParticleNumber  = ParticleNumber,
    epsilon         = epsilon,
    nparms          = nparms,
    UB              = UB,
    LB              = LB,
    EstMethod       = EstMethod,
    DependentVar    = DependentVar,
    Regressor       = Regressor,
    Task            = Task,
    OptMethod       = OptMethod,
    OutputType      = OutputType,
    ParamScheme     = ParamScheme,
    maxdiff         = maxdiff,
    loglik_BadValue1 = loglik_BadValue1,
    loglik_BadValue2 = loglik_BadValue2
  )
  return(out)

}

# PF likelihood with resampling for MA(q)
ParticleFilter_Res_ARMA = function(theta, mod){
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

  # Retrieve parameters and save them in a list called Parms
  Parms = RetrieveParameters(theta,mod)

  # check for causality and invertibility
  if( CheckStability(Parms$AR,Parms$MA) ){
    mod$ErrorMsg = sprintf('WARNING: The ARMA polynomial must be causal and invertible.')
    warning(mod$ErrorMsg)
    return(mod$loglik_BadValue1)
    }

  # Initialize the negative log likelihood computation
  nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
                   - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1,])))

  # retrieve AR, MA orders and their max
  m = max(mod$ARMAModel)
  p = mod$ARMAModel[1]
  q = mod$ARMAModel[2]


  # Compute ARMA covariance up to lag n-1
  a        = list()
  if(!is.null(Parms$AR)){
    a$phi = Parms$AR
  }else{
    a$phi = 0
  }
  if(!is.null(Parms$MA)){
    a$theta = Parms$MA
  }else{
    a$theta = 0
  }
  a$sigma2 = 1
  gamma    = itsmr::aacvf(a,mod$n)

  # Compute coefficients of Innovations Algorithm see 5.2.16 and 5.3.9 in in Brockwell Davis book
  IA       = InnovAlg(Parms, gamma, mod)
  Theta    = IA$thetas
  Rt       = sqrt(IA$v)

  # Get the n such that |v_n-v_{n-1}|< mod$maxdiff. check me: does this guarantee convergence of Thetas?
  nTheta   = IA$n
  Theta_n  = Theta[[nTheta]]

  # allocate matrices for weights, particles and predictions of the latent series
  w        = matrix(0, mod$n, mod$ParticleNumber)
  Z        = matrix(0, mod$n, mod$ParticleNumber)
  Zhat     = matrix(0, mod$n, mod$ParticleNumber)

  # initialize particle filter weights
  w[1,]    = rep(1,mod$ParticleNumber)

  # Compute the first integral limits Limit$a and Limit$b
  Limit    = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # Initialize the particles using N(0,1) variables truncated to the limits computed above
  Z[1,]    = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))


  for (t in 2:mod$n){

    # compute Zhat_t
    Zhat[t,] = ComputeZhat_t(m,Theta,Z,Zhat,t, Parms,p,q, nTheta, Theta_n)

    # Compute integral limits
    Limit = ComputeLimits(mod, Parms, t, Zhat[t,], Rt)

    # Sample truncated normal particles
    Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b, t, Zhat[t,], Rt)

    # update weights
    w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      # print(t)
      # print(w[t,])
      message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
      return(mod$loglik_BadValue2)
    }

    # Resample the particles using common random numbers
    old_state1 = get_rand_state()
    Znew = ResampleParticles(mod, w, t, Znew)
    set_rand_state(old_state1)

    # save the current particle
    Z[t,]   = Znew

    # update likelihood
    nloglik = nloglik - log(mean(w[t,]))

  }

  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  # if (nloglik==Inf | is.na(nloglik)){
  #   nloglik = 10^8
  # }


  return(nloglik)
}

# Optimization wrapper to fit PF likelihood with resampling
FitMultiplePF_Res = function(theta, mod){
  #====================================================================================#
  # PURPOSE       Fit the Particle Filter log-likelihood with resampling.
  #               This function maximizes the PF likelihood, nfit many times for nparts
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
  n        = length(mod$DependentVar)

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

      if(mod$Task == 'Optimization'){
        # run optimization for our model --no ARMA model allowed
        optim.output <- optimx(
          par     = theta,
          fn      = ParticleFilter_Res_ARMA,
          lower   = mod$LB,
          upper   = mod$UB,
          hessian = TRUE,
          method  = mod$OptMethod,
          mod     = mod)
      }else{
        optim.output = as.data.frame(matrix(rep(NA,8+length(theta)), ncol=8+length(theta)))
        names(optim.output) = c(mapply(sprintf, rep("p%.f",length(theta)), (1:length(theta)), USE.NAMES = FALSE),
                                "value",  "fevals", "gevals", "niter", "convcode",  "kkt1", "kkt2", "xtime")

        optim.output[,1:length(theta)] = theta
        t4 = tic()
        loglikelihood = ParticleFilter_Res_ARMA(theta,mod)
        t4 = tic()-t4
        optim.output[,(length(theta)+1)] = loglikelihood
        optim.output[,(length(theta)+2)] = 1
        optim.output[,(length(theta)+3)] = 1
        optim.output[,(length(theta)+4)] = 0
        optim.output[,(length(theta)+8)] = as.numeric(t4)
      }


      # save estimates, loglik value and diagonal hessian
      ParmEst[nfit*(k-1)+j,]  = as.numeric(optim.output[1:nparms])
      loglik[nfit*(k-1) +j]   = optim.output$value
      convcode[nfit*(k-1) +j] = optim.output$convcode
      kkt1[nfit*(k-1) +j]     = optim.output$kkt1
      kkt2[nfit*(k-1) +j]     = optim.output$kkt2


      # compute Hessian
        # t5 = tic()
      if(mod$Task == "Optimization"){
        H = gHgen(par          = ParmEst[nfit*(k-1)+j,],
                fn             = ParticleFilter_Res_ARMA,
                mod            = mod)
        # t5 = tic()-t5
        # print(t5)
        # if I get all na for one row and one col of H remove it
        # H$Hn[rowSums(is.na(H$Hn)) != ncol(H$Hn), colSums(is.na(H$Hn)) != nrow(H$Hn)]

        # t6 = tic()
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
        # t6 = tic()-t6
        # print(t6)

    }
  }

  # Compute model selection criteria (assuming one fit)
  Criteria = Criteria.lgc(loglik, mod)


  if(mod$OutputType=="list"){
    #  save the results in a list
    ModelOutput = list()

    # specify output list names
    # names(ModelOutput)         = c("ParamEstimates", "StdErrors", "FitStatistics", "OptimOutput")
    ModelOutput$Model          = paste(mod$CountDist,  paste("ARMA(",mod$ARMAModel[1],",",mod$ARMAModel[2],")",sep=""), sep= "-")
    ModelOutput$ParamEstimates = ParmEst
    ModelOutput$Task           = mod$Task
    if(!is.null(mod$TrueParam)) ModelOutput$TrueParam      = mod$TrueParam
    if(!is.null(mod$initialParam)) ModelOutput$initialParam   = mod$initialParam
    if(Task=="Optimization") ModelOutput$StdErrors      = se
    ModelOutput$FitStatistics  = c(loglik, Criteria)
    if(Task=="Optimization") ModelOutput$OptimOutput    = c(convcode,kkt1,kkt2)
    ModelOutput$EstMethod      = mod$EstMethod
    ModelOutput$SampleSize     = mod$n
    if(loglikelihood==mod$loglik_BadValue1) ModelOutput$WarnMessage = "WARNING: The ARMA polynomial must be causal and invertible."
    if(loglikelihood==mod$loglik_BadValue2) ModelOutput$WarnMessage = "WARNING: Some of the weights are either too small or sum to 0."
    # assign names to all output elements
    if(!is.null(mod$TrueParam)) names(ModelOutput$TrueParam)      = mod$parmnames
    colnames(ModelOutput$ParamEstimates) = mod$parmnames

    if(Task=="Optimization")colnames(ModelOutput$StdErrors)      = paste("se(", mod$parmnames,")", sep="")
    names(ModelOutput$FitStatistics)     = c("loglik", "AIC", "BIC", "AICc")
    if(Task=="Optimization")names(ModelOutput$OptimOutput)       = c("ConvergeStatus", "kkt1", "kkt2")

  }else{
    ModelOutput  = data.frame(matrix(ncol = 4*mod$nparms+16, nrow = 1))

    # specify output names
    if(!is.null(mod$TrueParam)){
      colnames(ModelOutput) = c(
        'CountDist','ARMAModel', 'Regressor',
        paste("True_", mod$parmnames, sep=""), paste("InitialEstim_", mod$parmnames, sep=""),
        mod$parmnames, paste("se(", mod$parmnames,")", sep=""),
        'EstMethod', 'SampleSize', 'ParticleNumber', 'epsilon', 'OptMethod', 'ParamScheme',
        "loglik", "AIC", "BIC", "AICc", "ConvergeStatus", "kkt1", "kkt2")
    }
    colnames(ModelOutput) = c(
      'CountDist','ARMAModel', 'Regressor',
      paste("InitialEstim_", mod$parmnames, sep=""),
      mod$parmnames, paste("se(", mod$parmnames,")", sep=""),
      'EstMethod', 'SampleSize', 'ParticleNumber', 'epsilon', 'OptMethod',
      "loglik", "AIC", "BIC", "AICc", "ConvergeStatus", "kkt1", "kkt2")

    # Start Populating the output data frame
    ModelOutput$CountDist      = mod$CountDist
    ModelOutput$ARMAModel      = paste("ARMA(",mod$ARMAModel[1],",",mod$ARMAModel[2],")",sep="")
    ModelOutput$Regressor      = !is.null(mod$Regressor)

    offset = 4
    # true Parameters
    if(!is.null(mod$TrueParam)){
      if(mod$nMargParms>0){
        ModelOutput[, offset:(offset + mod$nMargParms -1)] = mod$TrueParam[1:mod$nMargParms]
        offset = offset + mod$nMargParms
      }

      if(mod$nAR>0){
        ModelOutput[, offset:(offset + mod$nAR -1)]        = mod$TrueParam[ (mod$nMargParms+1):(mod$nMargParms+mod$nAR) ]
        offset = offset + mod$nAR
      }

      if(mod$nMA>0){
        ModelOutput[, offset:(offset + mod$nMA -1)]        = mod$TrueParam[ (mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA)]
        offset = offset + mod$nMA
      }
    }

    # Initial Parameter Estimates
    if(mod$nMargParms>0){
      ModelOutput[, offset:(offset + mod$nMargParms -1 )] = mod$initialParam[1:mod$nMargParms ]
      offset = offset + mod$nMargParms
    }
    if(mod$nAR>0){
      ModelOutput[, offset:(offset + mod$nAR-1)]        = mod$initialParam[(mod$nMargParms+1):(mod$nMargParms+mod$nAR) ]
      offset = offset + mod$nAR
    }

    if(mod$nMA>0){
      ModelOutput[, offset:(offset + mod$nMA-1)]        = mod$initialParam[(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA) ]
      offset = offset + mod$nMA
    }

    # Parameter Estimates
    if(mod$nMargParms>0){
      ModelOutput[, offset:(offset + mod$nMargParms-1)] = ParmEst[,1:mod$nMargParms]
      offset = offset + mod$nMargParms
    }

    if(mod$nAR>0){
      ModelOutput[, offset:(offset + mod$nAR-1)]        = ParmEst[,(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]
      offset = offset + mod$nAR
    }

    if(mod$nMA>0){
      ModelOutput[, offset:(offset + mod$nMA-1)]        = ParmEst[,(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA)]
      offset = offset + mod$nMA
    }

    # Parameter Std Errors
    if(mod$nMargParms>0){
      ModelOutput[, offset:(offset + mod$nMargParms-1)] = se[,1:mod$nMargParms]
      offset = offset + mod$nMargParms
    }

    if(mod$nAR>0){
      ModelOutput[, offset:(offset + mod$nAR-1)]        = se[,(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]
      offset = offset + mod$nAR
    }

    if(mod$nMA>0){
      ModelOutput[, offset:(offset + mod$nMA-1)]        = se[,(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA)]
    }

    ModelOutput$EstMethod      = mod$EstMethod
    ModelOutput$SampleSize     = mod$n
    ModelOutput$ParticleNumber = mod$ParticleNumber
    ModelOutput$epsilon        = mod$epsilon
    ModelOutput$OptMethod      = row.names(optim.output)
    if(!is.null(mod$TrueParam)) ModelOutput$ParamScheme    = mod$ParamScheme
    ModelOutput$loglik         = loglik
    ModelOutput$AIC            = Criteria[1]
    ModelOutput$BIC            = Criteria[2]
    ModelOutput$AICc           = Criteria[3]
    ModelOutput$ConvergeStatus = convcode
    ModelOutput$kkt1           = kkt1
    ModelOutput$kkt2           = kkt2
  }

  return(ModelOutput)
}

# --what to reru +
ErrorOutput =function(mod){

  # create messages for Error Codes
  ErrorMsg = switch(mod$ErrorCode,
                    "1"  = "The specified distribution in not supported.",
                    "2"  = "The specified Task is not supported.",
                    "3"  = "The specified Estimation method is not supported.",
                    "4"  = "The specified ARMA model is not supported.",
                    "5"  = "The initial parameter length does not match the specified model.",
                    "6"  = "The specified initial parameter isn't feasible."
  )

  ModelOutput = list()
  ModelOutput$ErrorMsg       = ErrorMsg

  # create final output list when we have an error
  # if(mod$OutputType=="list"){
  #   ModelOutput = list()
  #   ModelOutput$ErrorMsg       = ErrorMsg
  #   ModelOutput$Model          = paste(mod$CountDist,  paste("ARMA(",mod$ARMAModel[1],",",mod$ARMAModel[2],")",sep=""), sep= "-")
  #   ModelOutput$ParamEstimates = NULL
  #   ModelOutput$Task           = mod$Task
  #   if(!is.null(mod$TrueParam))    ModelOutput$TrueParam      = mod$TrueParam
  #   if(!is.null(mod$initialParam)) ModelOutput$initialParam   = mod$initialParam
  #   ModelOutput$EstMethod      = mod$EstMethod
  #   ModelOutput$SampleSize     = mod$n
  #   # assign names to all output elements
  #   if(!is.null(mod$TrueParam)) names(ModelOutput$TrueParam)      = mod$parmnames
  # }

  return(ModelOutput)

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
    rc = ifelse(any(abs( polyroot(c(1, -AR))  ) < 1) || any(abs( polyroot(c(1, -MA))  ) < 1)   , 1, 0)
  }

  return(rc)
}

# compute initial estimates
InitialEstimates = function(mod){
  # require(itsmr)
  est  = rep(NA, mod$nMargParms+sum(mod$ARMAModel))
  #-----Poisson case
  if(mod$CountDist=="Poisson"){
    if(mod$nreg==0){
      est[1] = mean(mod$DependentVar)

      # Transform (1) in the JASA paper to retrieve the "observed" latent series and fit an ARMA
      if(max(mod$ARMAModel)>0) armafit = itsmr::arma(qnorm(mod$mycdf(mod$DependentVar,est[1])),mod$nAR,mod$nMA)
      if(mod$nAR>0) est[(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]                    = armafit$phi
      if(mod$nMA>0) est[(1+mod$nMargParms+mod$nAR):(mod$nMargParms+sum(mod$ARMAModel))] = armafit$theta
    }else{
      # GLM for the mean that depends on time
      # CHECK ME: If I fit a Poisson AR(3) in the the data example of the JASA paper, but the code below doesn't specify poisson family (it would pick up the default distribution that glm function has) then there will be a numerical error in the likelihood. Check it!
      glmPoisson            = glm(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)], family = "poisson")
      est[1:mod$nMargParms] = as.numeric(glmPoisson[1]$coef)

      if(mod$nAR) est[(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]                    = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$phi
      if(mod$nMA) est[(1+mod$nAR+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$theta
    }
  }

  #-----Neg Binomial case
  if(mod$CountDist=="Negative Binomial"){
    if(mod$nreg==0){
      xbar = mean(mod$DependentVar)
      sSquare = var(mod$DependentVar)

      # Method of Moments for negBin
      rEst = xbar^2/(sSquare - xbar)
      pEst = 1 - xbar/sSquare
      est[1:2] = c(rEst, pEst)
      if(mod$nAR) est[(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]                    = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$phi
      if(mod$nMA) est[(1+mod$nMargParms+mod$nAR):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$theta
    }else{
      # GLM for the mean that depends on time
      glmNegBin                 = glm.nb(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)])
      est[1:(mod$nMargParms-1)] = as.numeric(glmNegBin[1]$coef)
      # Mom on constant variance
      est[mod$nMargParms]       = NegBinMoM(mod$DependentVar,glmNegBin$fitted.values)
      if(mod$nAR) est[(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]                    = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$phi
      if(mod$nMA) est[(1+mod$nAR+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$theta
    }
  }



  if(mod$CountDist=="Mixed Poisson"){
    if(mod$nreg==0){
      # pmle for marginal parameters
      MixPois_PMLE <- pmle.pois(mod$DependentVar,2)

      pEst  = MixPois_PMLE[[1]][1]
      l1Est = MixPois_PMLE[[2]][1]
      l2Est = MixPois_PMLE[[2]][2]


      # correct estimates if they are outside the feasible region
      # if(pEst<mod$LB[1]){pEst = 1.1*mod$LB[1]}
      # if(pEst>mod$UB[1]){pEst = 0.9*mod$UB[1]}
      #
      # if(l1Est<mod$LB[2]){l1Est = 1.1*mod$LB[2]}
      # if(l2Est<mod$LB[3]){l2Est = 1.1*mod$LB[3]}

      est[1:3] = c(l1Est, l2Est, pEst)

      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$nMargParms+mod$ARMAModel[1]):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }else{
      # library(mixtools)
      mix.reg = poisregmixEM(mod$DependentVar, mod$Regressor[,2:(mod$nreg+1)])
      est[1:mod$nMargParms] = c(mix.reg$beta, mix.reg$lambda[1])
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$ARMAModel[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }
  }



  #-----Generalized Poisson case
  if(mod$CountDist=="Generalized Poisson"){
    if(mod$nreg==0){
      xbar    = mean(mod$DependentVar)
      sSquare = var(mod$DependentVar)

      # Method of Moments for negBin
      rEst = xbar^2/(sSquare - xbar)
      pEst = 1 - xbar/sSquare
      est[1:2] = c(rEst, pEst)
      if(mod$nAR) est[(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]                      = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$phi
      if(mod$nMA) est[(1+mod$nMargParms+mod$nAR):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$theta
    }else{
      est[1:mod$nMargParms] = as.numeric(glm.nb(mod$DependentVar~mod$Regressor)[1]$coef)
      if(mod$nAR) est[(mod$nMargParms+1):(1+mod$nMargParms+mod$nAR)]                    = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$phi
      if(mod$nMA) est[(1+mod$nAR+mod$nMargParms):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$nAR,mod$nMA)$theta
    }
  }

  #-----Binomial case
  if(mod$CountDist=="Binomial"){
    if(mod$nreg==0){
      xbar = mean(mod$DependentVar)
      sSquare = var(mod$DependentVar)

      # Method of Moments for Binomial E(X)=np, Var(X)=np(1-p), 1-p = Var(X)/E(X)
      pEst = 1 - sSquare/xbar
      nEst = xbar/pEst
      est[1:2] = c(pEst, nEst)
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$nMargParms+mod$ARMAModel[1]):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }else{
      failures                  = max(mod$DependentVar) - mod$DependentVar
      glmBinomial               = glm(cbind(mod$DependentVar,failures)~mod$Regressor[,2:(mod$nreg+1)], family = "binomial")
      est[1:(mod$nMargParms-1)] = as.numeric(glmBinomial$coef)
      est[mod$nMargParms]       = max(mod$DependentVar)  ##Need to use method of moment?
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$ARMAModel[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }
  }

  #-----Zero Inflation Poisson
  if(mod$CountDist=="ZIP"){
    if(mod$nreg==0){
      # pmle for marginal parameters
      ZIP_PMLE <- poisson.zihmle(mod$DependentVar, type = c("zi"),
                                 lowerbound = LB,
                                 upperbound = 10000)

      lEst = ZIP_PMLE[1]
      pEst = ZIP_PMLE[2]

      # correct estimates if they are outside the feasible region
      if(pEst<LB[1]){pEst = 1.1*mod$LB[1]}
      if(pEst>UB[1]){pEst = 0.9*mod$UB[1]}

      if(lEst<LB[2]){l1Est = 1.1*mod$LB[2]}

      est[1:2] = c(lEst, pEst)

      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$nMargParms+mod$ARMAModel[1]):(1+mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }else{
      zeroinfl_reg = zeroinfl(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)])$coefficients
      est[1:mod$nMargParms] = as.numeric(c(zeroinfl_reg$count, zeroinfl_reg$zero))
      if(mod$ARMAModel[1]) est[(mod$nMargParms+1):(mod$nMargParms+mod$ARMAModel[1])] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$phi
      if(mod$ARMAModel[2]) est[(1+mod$ARMAModel[1]+mod$nMargParms):(mod$nMargParms+sum(mod$ARMAModel))] = itsmr::arma(mod$DependentVar,mod$ARMAModel[1],mod$ARMAModel[2])$theta
    }
  }
  # add the parmnames on theta fix me: does this affect performance?
  names(est) = mod$parmnames


  return(est)
}

NegBinMoM = function(data, GLMMeanEst){
  # the GLMMeanEst is the GLM estimate of the standard log-link
  # th formula below is standard MoM for the over dispersion param in NegBin2 parametrization
  PhiMomEst = sum(GLMMeanEst^2)/(sum((data-GLMMeanEst)^2-GLMMeanEst))
  return(PhiMomEst)
}

InnovAlg = function(Parms,gamma, mod) {
  # adapt the Innovation.algorithm if ITSMR

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

  phi     = Parms$AR
  sigma2  = 1
  N       = length(gamma)
  theta_r = c(1,Parms$MA,numeric(N))

  # Innovations algorithm

  p = ifelse(is.null(Parms$AR),0,length(Parms$AR))
  q = ifelse(is.null(Parms$MA),0,length(Parms$MA))
  m = max(p,q)

  Theta   = list()
  v       = rep(NA,N+1)
  v[1]    = kappa(1,1)
  StopCondition = FALSE
  n       = 1


  while(!StopCondition && n<N ) {
    Theta[[n]] <- rep(0,n)
    Theta[[n]][n] = kappa(n+1,1)/v[1]
    if(n>q && mod$nAR==0) Theta[[n]][n]= 0
    if(n>1){
      for (k in 1:(n-1)) {
        js <- 0:(k-1)
        Theta[[n]][n-k] <- (kappa(n+1,k+1) - sum(Theta[[k]][k-js]*Theta[[n]][n-js]*v[js+1])) / v[k+1]
      }
    }
    js     = 0:(n-1)
    v[n+1] = kappa(n+1,n+1) - sum(Theta[[n]][n-js]^2*v[js+1])
    # pure MA stopping criterion
    if(mod$nAR==0 && mod$nMA>0) StopCondition = (abs(v[n+1]-v[n])< mod$maxdiff)

    # pure AR stopping criterion
    if(mod$nAR>0 && mod$nMA==0) StopCondition = (n>3*m)

    # mixed ARMA stopping criterion
    if(mod$nAR>0 && mod$nMA>0)  StopCondition = (n>4*(p+q))
    n      = n+1
  }
  v = v/v[1]

  StopCondition = (abs(v[n+1]-v[n])<mod$maxdiff)


  return(list(n=n-1,thetas=lapply(Theta[ ][1:(n-1)], function(x) {x[1:m]}),v=v[1:(n-1)]))
}

# simulate from our model
sim_lgc = function(n, CountDist, MargParm, ARParm, MAParm, Regressor=NULL){


  # Generate latent Gaussian model
  z  =arima.sim(model = list( ar = ARParm, ma=MAParm  ), n = n)
  z = z/sd(z) # standardize the data

  # number of regressors assuming there is an intercept, fix me: check the case with not intercept
  nreg = ifelse(is.null(Regressor), 0, dim(Regressor)[2]-1)

  if(nreg==0){
    # retrieve marginal inverse cdf
    myinvcdf = switch(CountDist,
                      "Poisson"             = qpois,
                      "Negative Binomial"   = function(x, theta){ qnbinom (x, theta[1], 1-theta[2]) },
                      "Generalized Poisson" = function(x, theta){ qGpois  (x, theta[1], theta[2])},
                      "Binomial"            = qbinom,
                      "Mixed Poisson"       = function(x, theta){ qmixpois(x, theta[1], theta[2], theta[3])},
                      "ZIP"                 = function(x, theta){ qzipois(x, theta[1], theta[2]) },
                      stop("The specified distribution is not supported.")
    )

    # get the final counts
    x = myinvcdf(pnorm(z), MargParm)

  }else{
    # retrieve inverse count cdf
    myinvcdf = switch(CountDist,
                      "Poisson"             = function(x, ConstMargParm, DynamMargParm){ qpois   (x, DynamMargParm)},
                      "Negative Binomial"   = function(x, ConstMargParm, DynamMargParm){ qnbinom (x, ConstMargParm, 1-DynamMargParm)},
                      "Generalized Poisson" = function(x, ConstMargParm, DynamMargParm){ qGpois  (x, ConstMargParm, DynamMargParm)},
                      "Binomial"            = function(x, ConstMargParm, DynamMargParm){ qbinom  (x, ConstMargParm, DynamMargParm)},
                      "Mixed Poisson"       = function(x, ConstMargParm, DynamMargParm){ qmixpois(x, DynamMargParm, ConstMargParm)},
                      "ZIP"                 = function(x, ConstMargParm, DynamMargParm){ qzipois (x, DynamMargParm[1], DynamMargParm[2]) },
                      stop("The specified distribution is not supported.")
    )

    # regression parameters
    beta  = MargParm[1:(nreg+1)]
    m     = exp(Regressor%*%beta)

    if(CountDist == "Poisson" && nreg>0){
      ConstMargParm  = NULL
      DynamMargParm  = m
    }

    if(CountDist == "Negative Binomial" && nreg>0){
      ConstMargParm  = 1/MargParm[nreg+2]
      DynamMargParm  = MargParm[nreg+2]*m/(1+MargParm[nreg+2]*m)
    }

    if(CountDist == "Generalized Poisson" && nreg>0){
      ConstMargParm  = MargParm[nreg+2]
      DynamMargParm  = m
    }

    if(CountDist == "Binomial" && nreg>0){
      ConstMargParm  = MargParm[nreg+2]
      DynamMargParm  = 1/(1+m)
    }

    if(CountDist == "Mixed Poisson" && nreg>0){
      ConstMargParm  = c(MargParm[nreg*2+3], 1 - MargParm[nreg*2+3])
      DynamMargParm  = cbind(exp(Regressor%*%MargParm[1:(nreg+1)]),
                             exp(Regressor%*%MargParm[(nreg+2):(nreg*2+2)]))
    }

    if(CountDist == "ZIP" && nreg>0){
      ConstMargParm  = NULL
      DynamMargParm  = cbind(exp(Regressor%*%MargParm[1:(nreg+1)]),
                             1/(1+exp(-Regressor%*%MargParm[(nreg+2):(nreg*2+2)])))
    }

    # get the final counts
    x = myinvcdf(pnorm(z), ConstMargParm, DynamMargParm)
  }
    return(x)
}

# poisson cdf using incomplete gamma and its derivative wrt to lambda
myppois = function(x, lambda){
  # compute poisson cdf as the ratio of an incomplete gamma function over the standard gamma function
  # I will also compute the derivative of the poisson cdf wrt lambda
  X  =c(lambda,x+1)
  v1 = gammainc(X)
  v2 = gamma(x+1)

  # straight forward formula from the definition of incomplete gamma integral
  v1_d = -lambda^x*exp(-lambda)
  v2_d = 0

  z  = v1/v2
  z_d = (v1_d*v2 - v2_d*v1)/v2^2
  return(c(z,z_d))
}

#---------Compute AIC, BIC, AICc
Criteria.lgc = function(loglik, mod){
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

logLik.lgc = function(object){
  return(object$FitStatistics[1])
}

AIC.lgc = function(object){
  return(object$FitStatistics[2])
}

BIC <- function(object, ...) UseMethod("BIC")

BIC.lgc = function(object){
  return(object$FitStatistics[3])
}

se <- function(object, ...) UseMethod("se")

se.lgc = function(object){
  return(object$StdErrors)
}

coefficients <- function(object, ...) UseMethod("coefficients")

coefficients.lgc = function(object){
  return(object$ParamEstimates)
}

model <- function(object, ...) UseMethod("model")

model.lgc = function(object){
  if ((object$ARMAModel[1]>0) &&  (object$ARMAModel[2]>0)){
    ARMAModel = sprintf("ARMA(%.0f, %.0f)",object$ARMAModel[1], object$ARMAModel[2])
  }
  if ((object$ARMAModel[1]>0) &&  (object$ARMAModel[2]==0)){
    ARMAModel = sprintf("AR(%.0f)",object$ARMAModel[1])
  }
  if ((object$ARMAModel[1]==0) &&  (object$ARMAModel[2]>0)){
    ARMAModel = sprintf("MA(%.0f)",object$ARMAModel[2])
  }
  if ((object$ARMAModel[1]==0) &&  (object$ARMAModel[2]==0)){
    ARMAModel = "White Noise"
  }

  a = data.frame(object$CountDist, ARMAModel)
  names(a) = c("Distribution", "Model")
  return(a)
}


ComputeLimits = function(mod, Parms, t, Zhat, Rt){
  # a and b are the arguments in the two normal cdfs in the 4th line in equation (19) in JASA paper
  Lim = list()
  # fix me: this will be ok for AR or MA models but for ARMA? is it p+q instead of max(p,q)
  index = min(t, max(mod$ARMAModel))
  if(mod$nreg==0){
    Lim$a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,t(Parms$MargParms)),0,1)) - Zhat)/(Rt[index])
    Lim$b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],t(Parms$MargParms)),0,1)) - Zhat)/Rt[index]
  }else{
    Lim$a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,Parms$ConstMargParm, Parms$DynamMargParm[t,]),0,1)) - Zhat)/Rt[index]
    Lim$b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],Parms$ConstMargParm, Parms$DynamMargParm[t,]),0,1)) - Zhat)/Rt[index]
  }

  return(Lim)
}

SampleTruncNormParticles = function(mod, a, b, t, Zhat, Rt){
  # relation (21) in JASA paper and the inverse transform method
  # check me: this can be improved?
  index = min(t, max(mod$ARMAModel))
  z = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)*Rt[index] + Zhat
  return(z)
}

ComputeZhat_t = function(m,Theta,Z,Zhat,t, Parms,p,q,nTheta, Theta_n){

  if(m>1 && t<=m) Zhat_t = Theta[[t-1]][1:(t-1)]%*%(Z[(t-1):1,]-Zhat[(t-1):1,])

  if(t>m && t<=nTheta){
    A = B= 0
    if(!is.null(Parms$AR)) A = Parms$AR%*%Z[(t-1):(t-p),]
    if(!is.null(Parms$MA)) B = Theta[[t-1]][1:q]%*%(Z[(t-1):(t-q),]-Zhat[(t-1):(t-q),])
    Zhat_t = A + B
  }
  if(t>nTheta){
    A = B = 0
    if(!is.null(Parms$AR)) A = Parms$AR%*%Z[(t-1):(t-p),]
    if(!is.null(Parms$MA)) B = Theta_n[1:q]%*%(Z[(t-1):(t-q),]-Zhat[(t-1):(t-q),])
    Zhat_t = A + B
  }

  return(Zhat_t)
}


ComputeWeights = function(mod, a, b, t, PreviousWeights){
  # equation (21) in JASA paper
  # update weights
  if(t<=max(mod$ARMAModel)){
    NewWeights = PreviousWeights*(pnorm(b,0,1) - pnorm(a,0,1))
  }else{ # fix me: if I add the wgh[t-1,] below as I should the weights become small?
    NewWeights = (pnorm(b,0,1) - pnorm(a,0,1))
  }

  return(NewWeights)
}

ResampleParticles = function(mod, wgh, t, Znew){

  # relation (26) in JASA paper and following step
  # compute normalized weights
  wghn = wgh[t,]/sum(wgh[t,])

  # effective sample size
  ESS = 1/sum(wghn^2)

  if(ESS<mod$epsilon*mod$ParticleNumber){
    ind = rmultinom(1,mod$ParticleNumber,wghn)
    Znew = rep(Znew,ind)
  }
return(Znew)
}

RetrieveParameters = function(theta,mod){

  Parms =  vector(mode = "list", length = 5)


  names(Parms) = c("MargParms", "ConstMargParm", "DynamMargParm", "AR", "MA")

  # marginal parameters
  Parms$MargParms      = theta[mod$MargParmIndices]

  # regressor parameters
  if(mod$nreg>0){
    beta  = Parms$MargParms[1:(mod$nreg+1)]
    m     = exp(mod$Regressor%*%beta)
  }

  # GLM type parameters
  if(mod$CountDist == "Negative Binomial" && mod$nreg>0){
    Parms$ConstMargParm  = 1/Parms$MargParms[mod$nreg+2]
    Parms$DynamMargParm  = Parms$MargParms[mod$nreg+2]*m/(1+Parms$MargParms[mod$nreg+2]*m)
  }

  if(mod$CountDist == "Generalized Poisson" && mod$nreg>0){
    Parms$ConstMargParm  = Parms$MargParms[mod$nreg+2]
    Parms$DynamMargParm  = m
  }

  if(mod$CountDist == "Poisson" && mod$nreg>0){
    Parms$ConstMargParm  = NULL
    Parms$DynamMargParm  = m
  }

  if(mod$CountDist == "Binomial" && mod$nreg>0){
    Parms$ConstMargParm  = Parms$MargParms[mod$nreg+2]
    Parms$DynamMargParm  = 1/(1+m)
  }

  if(mod$CountDist == "Mixed Poisson" && mod$nreg>0){
    Parms$ConstMargParm  = c(Parms$MargParms[mod$nreg*2+3], 1 - Parms$MargParms[mod$nreg*2+3])
    Parms$DynamMargParm  = cbind(exp(mod$Regressor%*%Parms$MargParms[1:(mod$nreg+1)]),
                           exp(mod$Regressor%*%Parms$MargParms[(mod$nreg+2):(mod$nreg*2+2)]))
  }

  if(mod$CountDist == "ZIP" && mod$nreg>0){
    Parms$ConstMargParm  = NULL
    Parms$DynamMargParm  = cbind(exp(mod$Regressor%*%Parms$MargParms[1:(mod$nreg+1)]),
                           1/(1+exp(-mod$Regressor%*%Parms$MargParms[(mod$nreg+2):(mod$nreg*2+2)])))
  }


  # Parms$DynamMargParm = as.matrix(Parms$DynamMargParm)

  # Parms$AR = NULL
  if(mod$ARMAModel[1]>0) Parms$AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAModel[1])]

  # Parms$MA = NULL
  if(mod$ARMAModel[2]>0) Parms$MA = theta[(mod$nMargParms+mod$ARMAModel[1]+1) :
                                      (mod$nMargParms + mod$ARMAModel[1] + mod$ARMAModel[2]) ]


  return(Parms)
}


# Mixed Poisson inverse cdf
qmixpois = function(y, lam1, lam2, p){
  yl = length(y)
  x  = rep(0,yl)
  for (n in 1:yl){
    while(pmixpois(x[n], lam1, lam2,p) <= y[n]){ # R qpois would use <y; this choice makes the function right-continuous; this does not really matter for our model
      x[n] = x[n]+1
    }
  }
  return(x)
}

# Mixed Poisson cdf
pmixpois = function(x, lam1, lam2,p){
  y = p*ppois(x,lam1) + (1-p)*ppois(x,lam2)
  return(y)
}

# mass of Mixed Poisson
dmixpois = function(x, lam1, lam2, p){
  y = p*dpois(x,lam1) + (1-p)*dpois(x,lam2)
  return(y)
}

# Generalized Poisson pdfs
dGpois = function(y,a,m){
  k = m/(1+a*m)
  return( k^y * (1+a*y)^(y-1) * exp(-k*(1+a*y)-lgamma(y+1)))
}

# Generalized Poisson cdf for one m---------#
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

# function that generates random marginal parameter values
GenModelParam = function(CountDist,BadParamProb, AROrder, MAOrder, Regressor){

  #-----------------Marginal Parameters

  if(CountDist=="Poisson"){
    # create some "bad" Poisson parameter choices
    BadLambda = c(-1,500)

    # sample with high probability from "good" choices for lambda and with low prob the "bad" choices
    prob = rbinom(1,1,BadParamProb)

    # Marginal Parameter
    if (is.null(Regressor)){
      MargParm = prob*runif(1,0,100) + (1-prob)*sample(BadLambda,1)
    }else{
      # fix me: I am hard coding 1.2 here but we should probably make this an argument
      b0 = runif(1,0,1.2)
      b1 = runif(1,0,1.2)
      MargParm = c(b0,b1)
    }
  }


  if(CountDist=="Negative Binomial"){
    # create some "bad" Poisson parameter choices
    Badr_r = c(-1,500)

    # sample a probability
    p =  runif(1,0,1)

    # sample with high probability from "good" choices for lambda and with low prob the "bad" choices
    prob = rbinom(1,1,BadParamProb)

    # Marginal Parameter
    if (is.null(Regressor)){
      r = prob*runif(1,0,100) + (1-prob)*sample(Badr_r,1)
      MargParm = c(r,p)
    }else{
      # fix me: I am hard coding 1.2 here but we should probably make this an argument
      b0 = runif(1,0,1.2)
      b1 = runif(1,0,1.2)
      k  = runif(1,0,1)
      MargParm = c(b0,b1,k)
    }
  }


  #------------------ARMA Parameters

  # create some "bad" AR parameter choices
  BadAR = c(0, -10, 0.99, 1.2)

  # create some "bad" MA parameter choices
  BadMA = c(0, 10, -0.99, -1.2)

  # set the AR Parameters
  ARParm = NULL
  MAParm = NULL
  p = rbinom(1,1,BadParamProb)

  if(AROrder) ARParm = p*runif(1,-1, 1)+ (1-p)*sample(BadAR,1)
  if(MAOrder) MAParm = p*runif(1,-1, 1)+ (1-p)*sample(BadMA,1)

  AllParms = list(MargParm, ARParm, MAParm)
  names(AllParms) = c("MargParm", "ARParm", "MAParm")
  return(AllParms)

}

GenInitVal = function(AllParms, perturbation){

  # select negative or postive perturbation based on a  coin flip
  #sign = 1-2*rbinom(1,1,0.5)
  # adding only negative sign now to ensure stability
  sign  = -1

  MargParm = AllParms$MargParm*(1+sign*perturbation)

  PerturbedValues = list(MargParm,NULL,NULL)
  names(PerturbedValues) = c("MargParm", "ARParm", "MAParm")

  if(!is.null(AllParms$ARParm))  PerturbedValues$ARParm = AllParms$ARParm*(1+sign*perturbation)
  if(!is.null(AllParms$MAParm))  PerturbedValues$MAParm = AllParms$MAParm*(1+sign*perturbation)

  return(PerturbedValues)
}

# innovations algorithm code
innovations.algorithm <- function(acvf,n.max=length(acvf)-1){
  # I found this implementation of IA online, I need to check it
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

# PF likelihood with resampling for AR(p)
ParticleFilter_Res_AR_Old = function(theta, mod){
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


  #-----------  Step 0: Retrieve values from the mod Structure --------------#

  # marginal parameters
  MargParms        = theta[mod$MargParmIndices]

  # regressor parameters
  if(mod$nreg>0){
    beta  = MargParms[1:(mod$nreg+1)]
    m     = exp(mod$Regressor%*%beta)
  }

  # GLM type parameters
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

  # ARMA parameters
  AR = NULL
  if(mod$ARMAModel[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAModel[1])]

  MA = NULL
  if(mod$ARMAModel[2]>0) MA = theta[(mod$nMargParms+mod$ARMAModel[1]+1) :
                                      (mod$nMargParms + mod$ARMAModel[1] + mod$ARMAModel[2]) ]

  # check for causality
  if( CheckStability(AR,MA) ) return(10^(8))


  # sample size and number of particles
  T1      = length(mod$DependentVar)
  N       = mod$ParticleNumber

  # Initialize the negative log likelihood computation
  nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],MargParms)),
                   - log(mod$mypdf(mod$DependentVar[1], ConstMargParm, DynamMargParm[1])))

  # Compute the theoretical covariance for the AR model for current estimate
  gt      = ARMAacf(ar = AR, ma = MA)[2:(max(mod$ARMAModel)+1)]

  # Compute the best linear predictor coefficients and errors using Durbin Levinson
  DL      = DLAcfToAR(gt)
  phit    = DL[,1]
  Rt      = sqrt(DL[,3])


  # allocate memory for particle weights and the latent Gaussian Series particles
  wgh     = matrix(0,T1,N)
  Z       = matrix(0,mod$ARMAModel[1],N)

  #======================   Start the SIS algorithm   ======================#
  # Initialize the weights and the latent Gaussian series particles
  wgh[1,] = rep(1,N)

  if(mod$nreg==0){
    a       = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1,t(MargParms)),0,1), N)
    b       = rep( qnorm(mod$mycdf(mod$DependentVar[1],t(MargParms)),0,1), N)
  }else{
    a       = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1,ConstMargParm, DynamMargParm[1,]),0,1), N)
    b       = rep( qnorm(mod$mycdf(mod$DependentVar[1],ConstMargParm, DynamMargParm[1,]),0,1), N)
  }

  Z[1,]   = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

  # =================== Loop from 2 to AR order===================== #
  if (mod$ARMAModel[1]>=2){
    for (t in 2: (mod$ARMAModel[1])){
      # STEP 1 in SIS: Compute the latent Gaussian predictions Zhat using Durbin Levinson
      if (t==2) {
        Zhat = Z[1,]*phit[1]
      } else{
        Zhat = colSums(Z[1:(t-1),]*phit[1:(t-1)])
      }

      # STEP 2 is SIS: Update the latent Gaussian series Z and the importance weights w
      if(mod$nreg==0){
        a = (qnorm(mod$mycdf(mod$DependentVar[t]-1,t(MargParms)),0,1) - Zhat)/Rt[t]
        b = (qnorm(mod$mycdf(mod$DependentVar[t],t(MargParms)),0,1) - Zhat)/Rt[t]
      }else{
        a = (qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - Zhat)/Rt[t]
        b = (qnorm(mod$mycdf(mod$DependentVar[t],ConstMargParm, DynamMargParm[t]),0,1) - Zhat)/Rt[t]
      }

      Z[1:t,] = rbind(qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)*Rt[t] + Zhat, Z[1:(t-1),])
      wgh[t,] = wgh[t-1,]*(pnorm(b,0,1) - pnorm(a,0,1))

      # update likelihood
      nloglik = nloglik - log(mean(wgh[t,]))
      # print(t)
      # print(nloglik)
    }
  }
  # =================== Loop from AR order + 1  to T ===================== #
  # From p to T1 I don't need to estimate phi anymore
  for (t in (mod$ARMAModel[1]+1):T1){

    # STEP 1 in SIS: Compute the latent Gaussian predictions Zhat using Durbin Levinson
    if(mod$ARMAModel[1]>1){# colsums doesnt work for 1-dimensional matrix
      Zhat = colSums(Z*phit)
    }else{
      Zhat =  Z*phit
    }

    # STEP 2 is SISR: Update the latent Gaussian series Z
    if(mod$nreg==0){
      a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,t(MargParms)),0,1)) - Zhat)/Rt[mod$ARMAModel[1]]
      b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],t(MargParms)),0,1)) - Zhat)/Rt[mod$ARMAModel[1]]
    }else{
      a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1)) - Zhat)/Rt[mod$ARMAModel[1]]
      b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],ConstMargParm, DynamMargParm[t]),0,1)) - Zhat)/Rt[mod$ARMAModel[1]]
    }

    Znew = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)*Rt[mod$ARMAModel[1]] + Zhat

    # compute unnormalized weights
    # wgh[t,] = wgh[t-1,]*(pnorm(b,0,1) - pnorm(a,0,1))
    wgh[t,] = (pnorm(b,0,1) - pnorm(a,0,1))
    # break if I got NA weight
    if (any(is.na(wgh[t,]))| sum(wgh[t,])==0 ){
      message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
      return(10^8)
    }

    # compute normalized weights
    wghn = wgh[t,]/sum(wgh[t,])

    # STEP 3 is SISR: Resample
    old_state1 <- get_rand_state()
    ESS = 1/sum(wghn^2)
    if(ESS<mod$epsilon*N){
      ind = rmultinom(1,N,wghn)
      # sample particles
      Znew = rep(Znew,ind)
    }
    set_rand_state(old_state1)


    # save particles
    if (mod$ARMAModel[1]>1){
      Z = rbind(Znew, Z[1:(mod$ARMAModel[1]-1),])
    }else {
      Z[1,]=Znew
    }
    # update likelihood
    nloglik = nloglik - log(mean(wgh[t,]))
    # print(t)
    # print(nloglik)
  }

  return(nloglik)
}

# PF likelihood with resampling for MA(q)
ParticleFilter_Res_MA_Old = function(theta, mod){
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
  if(mod$ARMAModel[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAModel[1])]

  MA = NULL
  if(mod$ARMAModel[2]>0) MA = theta[(mod$nMargParms+mod$ARMAModel[1]+1) : (mod$nMargParms + mod$ARMAModel[1] + mod$ARMAModel[2]) ]

  # check for causality
  if( CheckStability(AR,MA) ) return(10^(-6))


  T1 = length(mod$DependentVar)
  N = mod$ParticleNumber          # number of particles


  # allocate matrix to collect all particle weights
  wgh = matrix(0,length(mod$DependentVar),N)

  # Compute integral limits
  if(mod$nreg==0){
    a = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1,t(MargParms)),0,1), N)
    b = rep( qnorm(mod$mycdf(mod$DependentVar[1],t(MargParms)),0,1), N)
  }else{
    a = rep( qnorm(mod$mycdf(mod$DependentVar[1]-1, ConstMargParm, DynamMargParm[1]) ,0,1), N)
    b = rep( qnorm(mod$mycdf(mod$DependentVar[1], ConstMargParm, DynamMargParm[1]) ,0,1), N)
  }

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)


  # run innovations Algorithm for MA models that are not WN
  if(mod$ARMAModel[2]>0) Inn = matrix(0,N,mod$ARMAModel[2])   # I will save here the q many innovations (Z - Zhat) --see (5.3.9) BD book
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
    # zhat = -Theta[[1]][1]*zprev
    zhat = Theta[[1]][1]*zprev
    Inn[,ARMAModel[2]] = zprev
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
    if(mod$ARMAModel[2]>1) Inn[,1:(mod$ARMAModel[2]-1)] = Inn[,2:(mod$ARMAModel[2])]

    # compute limits of truncated normal distribution
    if(mod$nreg==0){
      a = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t]-1,MargParms),0,1) - zhat)/vt
      b = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t],MargParms),0,1) -   zhat)/vt
    }else{
      a = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t]-1,ConstMargParm, DynamMargParm[t]),0,1) - zhat)/vt
      b = as.numeric(qnorm(mod$mycdf(mod$DependentVar[t],ConstMargParm, DynamMargParm[t]),0,1) -   zhat)/vt
    }

    # draw errors from truncated normal
    err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)

    # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
    znew = zhat + vt*err

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

    # compute new innovation
    Inn[,mod$ARMAModel[2]] = (znew-zhat)

    # update zhat--fix me can probably be vectorized
    if (is.null(MA) && is.null(AR)){
      zhat = 0
    }else{
      S = 0
      for(j in 1:min(t,mod$ARMAModel[2])){
        S = S+Theta[[t]][j]*Inn[,mod$ARMAModel[2]-j+1]
      }
      zhat = S
    }

    set_rand_state(old_state1)

    # update likelihood
    nloglik = nloglik - log(mean(wgh[t,]))
  }

  # likelihood approximation
  if(mod$nreg<1){
    nloglik = nloglik - log(mod$mypdf(mod$DependentVar[1],MargParms))
  }else{
    nloglik = nloglik - log(mod$mypdf(mod$DependentVar[1], ConstMargParm, DynamMargParm[1]))
  }

  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))


  if (nloglik==Inf | is.na(nloglik)){
    nloglik = 10^8
  }


  return(nloglik)
}


# PF likelihood with resampling for AR(p) - written more concicely
ParticleFilter_Res_AR = function(theta, mod){
  #--------------------------------------------------------------------------#
  # PURPOSE:  Use particle filtering with resampling
  #           to approximate the likelihood of the
  #           a specified count time series model with an underlying AR(p)
  #           dependence structure.
  #
  #
  # INPUTS:
  #    theta: parameter vector
  #      mod: a list containing all t he information for the model, such as
  #           count distribution. ARMA model, etc
  # OUTPUT:
  #    loglik: approximate log-likelihood
  #
  # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias
  # DATE:    July  2020
  #--------------------------------------------------------------------------#

  # keep track of the random seed to use common random numbers
  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # Retrieve parameters ans save them in a li
  Parms = RetrieveParameters(theta,mod)

  # check for causality
  if( CheckStability(Parms$AR,Parms$MA) ){
    message(sprintf('WARNING: The ARMA polynomial must be causal and invertible.'))
    return(mod$loglik_BadValue1)
  }

  # Initialize the negative log likelihood computation
  nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
                   - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1,])))

  # Compute the theoretical covariance for the AR model for current estimate
  gt    = ARMAacf(ar = Parms$AR, ma = Parms$MA,lag.max = mod$n)

  # Compute the best linear predictor coefficients and errors using Durbin Levinson
  Phi = list()
  for (t in 2:mod$n){
    CurrentDL     = DLAcfToAR(gt[2:t])
    Phi[[t-1]] = CurrentDL[,1]
  }
  Rt           = c(1,sqrt(as.numeric(CurrentDL[,3])))

  # allocate memory for particle weights and the latent Gaussian Series particles, check me: do I weights for 1:T or only 2?
  w     = matrix(0, mod$n, mod$ParticleNumber)
  Z     = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)

  #======================   Start the SIS algorithm   ======================#
  # Initialize the weights
  w[1,] = rep(1,mod$ParticleNumber)

  # Compute the first integral limits Limit$ a and Limit$b
  Limit = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # Initialize particles from truncated normal distribution
  Z[1,] = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # =================== Loop over t ===================== #
  for (t in 2:mod$n){
    # Compute the latent Gaussian predictions Zhat_t using Innovations Algorithm
    Zhat  =         Phi[[t-1]]%*%Z[1:(t-1),]

    # Compute integral limits
    Limit = ComputeLimits(mod, Parms, t, Zhat, Rt)

    # Sample truncated normal particles
    Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b,t, Zhat, Rt)

    # update weights
    w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
      return(mod$loglik_BadValue2)
    }

    # Resample the particles using common random numbers
    old_state1 = get_rand_state()
    Znew = ResampleParticles(mod, w, t, Znew)
    set_rand_state(old_state1)

    # Combine current particles, with particles from previous iterations
    Z = rbind(Znew, matrix(Z[1:(t-1),],ncol = mod$ParticleNumber))

    # update log-likelihood
    nloglik = nloglik - log(mean(w[t,]))
  }

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

  # Retrieve parameters and save them in a list called Parms
  Parms = RetrieveParameters(theta,mod)

  # check for causality
  if( CheckStability(Parms$AR,Parms$MA) ){
    message(sprintf('WARNING: The ARMA polynomial must be causal and invertible.'))
    return(mod$loglik_BadValue1)
  }


  # Initialize the negative log likelihood computation
  nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
                   - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1,])))

  # Compute covariance up to lag n-1
  gt    = as.vector(ARMAacf(ar = Parms$AR, ma = Parms$MA, lag.max = mod$n))

  # Compute coefficients of Innovations Algorithm see 5.2.16 and 5.3.9 in in Brockwell Davis book
  IA    = innovations.algorithm(gt)
  Theta = IA$thetas
  Rt    = sqrt(IA$vs)

  # allocate matrices for weights, particles and innovations which are equal to Z-Zhat
  w     = matrix(0, mod$n, mod$ParticleNumber)
  Z     = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)
  Inn   = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)

  # particle filter weights
  w[1,]   = rep(1,mod$ParticleNumber)

  # Compute the first integral limits Limit$ a and Limit$b
  Limit = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
  Z[1,]   = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # Compute the first innovation (Zhat_1=0)
  Inn[1,] = Z[1,]


  for (t in 2:mod$n){

    # Compute the latent Gaussian predictions Zhat_t using Innovations Algorithm - see 5.3.9 in Brockwell Davis book
    if(mod$ParticleNumber==1){
      if(t==2){
        Zhat  =         Inn[1:(min(t-1,mod$nMA)),] %*% Theta[[t-1]][1:(min(t-1,mod$nMA))]
      }else{
        Zhat  = colSums(Inn[1:(min(t-1,mod$nMA)),] %*% Theta[[t-1]][1:(min(t-1,mod$nMA))])
      }
    }else{
      if(t==2){
        Zhat  =         Inn[1:(min(t-1,mod$nMA)),] * Theta[[t-1]][1:(min(t-1,mod$nMA))]
      }else{
        Zhat  = colSums(Inn[1:(min(t-1,mod$nMA)),] * Theta[[t-1]][1:(min(t-1,mod$nMA))])
      }
    }

    # Compute integral limits
    Limit = ComputeLimits(mod, Parms, t, Zhat, Rt)

    # Sample truncated normal particles
    Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b, t, Zhat, Rt)

    # update weights
    w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
      return(mod$loglik_BadValue2)
    }

    # Resample the particles using common random numbers
    old_state1 = get_rand_state()
    Znew = ResampleParticles(mod, w, t, Znew)
    set_rand_state(old_state1)

    # Compute the new Innovation
    InnNew = Znew - Zhat

    # Combine current particles, with particles from previous iterations
    Inn[1:min(t,mod$nMA),] = rbind(matrix(InnNew,ncol=mod$ParticleNumber), matrix(Inn[1:min(t-1,mod$nMA-1),],ncol = mod$ParticleNumber))

    # update likelihood
    nloglik = nloglik - log(mean(w[t,]))

  }

  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  # if (nloglik==Inf | is.na(nloglik)){
  #   nloglik = 10^8
  # }


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

  loglik = ParticleFilter_Res_ARMA(theta, mod)

  # # Pure AR model - uses DL
  # if(mod$ARMAModel[1]>0 && mod$ARMAModel[2]==0 )  loglik = ParticleFilter_Res_AR(theta, mod)
  #
  # # Pure MA model - uses a non optimized INALg
  # if(mod$ARMAModel[1]==0 && mod$ARMAModel[2]>=0 ) loglik = ParticleFilter_Res_MA(theta, mod)
  #
  # # ARMA model - uses optimized InAlg
  # if(mod$ARMAModel[1]>0 && mod$ARMAModel[2]>0 )   loglik = ParticleFilter_Res_ARMA(theta, mod)

  return(loglik)
}



