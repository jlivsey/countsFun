#' #' @keywords internal
#' @importFrom stats ARMAacf arima.sim coef dbinom dnbinom dpois filter glm arima residuals
#' @importFrom stats pbinom plogis pnbinom pnorm ppois qbinom qnbinom qnorm nobs
#' @importFrom stats qpois rbinom rmultinom runif sd terms var ts
#' @importFrom MASS glm.nb
#' @importFrom VGAM vglm genpoisson2 loglink
#' @importFrom VGAM dgenpois2 qgenpois2 pgenpois2 rgenpois2 zipoisson
#' @importFrom iZID poisson.zihmle
#' @importFrom mixtools poisregmixEM
#' @importFrom extraDistr dmixpois pmixpois rmixpois dzip pzip qzip
#' @importFrom optimx optimx
#' @importFrom optextras gHgen
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom truncnorm rtruncnorm
"_PACKAGE"


#' Create and Validate Model Specification for Particle Filtering
#'
#' Parses and validates the full model specification and specified options.
#' This function gathers model components (data, distribution,
#' dependence structure, optimization parameters, etc.) and returns a structured list used as
#' input to other core functions in the package.
#'
#' @param DependentVar Numeric vector or data frame. The response (count) time series.
#' @param Regressor Optional matrix or data frame of covariates.
#' @param Intercept Logical. Whether to include an intercept in the model.
#' @param EstMethod Character. Estimation method. Currently only \code{"PFR"} is supported.
#' @param ARMAModel Integer vector of length 2. Specifies the order of the AR and MA components.
#' @param CountDist Character. The name of the count distribution (e.g., "Poisson", "Negative Binomial", "ZIP").
#' @param ParticleNumber Integer. Number of particles to be used in the filter.
#' @param epsilon Numeric. Resampling occurs when the effective sample size (ESS) drops below \code{epsilon * N}.
#' @param initialParam Numeric vector. Initial parameter values for optimization.
#' @param TrueParam Numeric vector. True parameter values (used in simulation/synthesis, ignored otherwise).
#' @param Task Character. One of \code{"Evaluation"}, \code{"Optimization"}, \code{"Simulation"}, or \code{"Synthesis"}.
#' @param SampleSize Integer. Required when \code{Task == "Synthesis"}.
#' @param OptMethod Character. Optimization method, e.g., \code{"L-BFGS-B"}.
#' @param OutputType Character. Type of output, typically \code{"list"}.
#' @param ParamScheme Integer. Indicates the parameter configuration scheme.
#' @param maxdiff Numeric. Convergence criterion for optimization.
#' @param ntrials Integer. Number of trials (used in Binomial model).
#' @param verbose Logical. If \code{TRUE} (default), informative messages are printed during execution.
#' @param nsim Integer. Number of replications. Required when \code{Task == "Simulation"}.
#' Set to \code{FALSE} to suppress messages.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named list containing model specifications, validated inputs, marginal distribution
#' functions (PDF, CDF, inverse CDF), initial and true parameters, parameter bounds,
#' parameter names, and utility values for likelihood evaluation.
#'
#' @details
#' This function validates the inputs (e.g., distribution support, ARMA structure, initial
#' parameter lengths), computes the number of model parameters, constructs marginal distribution
#' functions based on the count distribution, and configures parameter bounds. It is a required
#' pre-processing step before fitting or simulating count time series models.
#'
#'
#' @export
ModelScheme = function(DependentVar = NULL, Regressor=NULL, Intercept = NULL, EstMethod="PFR", ARMAModel=c(0,0), CountDist="NULL",
                       ParticleNumber = 5, epsilon = 0.5, initialParam=NULL, TrueParam=NULL, Task="Evaluation", SampleSize=NULL,
                       OptMethod="L-BFGS-B", OutputType="list", ParamScheme=1, maxdiff=10^(-8),ntrials= NULL,verbose=TRUE,nsim=NULL,...){

  # Distribution list
  if( !(CountDist %in% c("Poisson", "Negative Binomial", "Generalized Poisson", "Generalized Poisson 2", "Mixed Poisson",
                         "ZIP", "Binomial")))
    stop("The specified distribution in not supported.")

  # Specify Task
  #if( !(Task %in% c("Evaluation", "Optimization", "Synthesis"))){
  if( !(Task %in% c("Evaluation", "Optimization", "Simulation", "Synthesis"))){
    Task = "Evaluation"
    message("The selected Task is not supported. The Task has been set to 'Evaluation'")
  }

  # Estimation Method for Optimization and
  if( !(EstMethod %in% c("PFR")) & Task %in% c("Optimization", "Simulation") ){
    EstMethod = "PFR"
    message("The selected EstMethod is not supported. EstMethod has been set to'PFR'")
  }

  # If user specified EstMethod but asks for evaluation Task then no estimation takes place
  # check me: I will reconsider this, it seems a not very high priority comunication to the user,
  # it  will spame me in test.
    # if (Task %in% c("Evaluation") & EstMethod %in% c("PFR")){
  #   EstMethod = "None"
  #   message("For the Evaluation Task no estimation takes place. EstMethod has been set to None.\n",
  #   "To fit the model specify the Optimzation Task.")
  # }


  # check that the ARMA order has dimension 2, the ARMA orders are integers
  if(length(ARMAModel)!=2 | !prod(ARMAModel %% 1 == c(0,0)) )
    stop("The specified ARMA model is not supported.")

  # if evaluation or Optimization is requested dependent variable cannot be NULL
  if (Task %in% c("Evaluation", "Optimization") && is.null(DependentVar)){
    stop("Dependent variable is NULL")
  }

  # retrieve sample size and make the dependent variable numeric if Evaluation or Optimization has been requested
  if(Task %in% c("Synthesis", "Simulation")){
    n = SampleSize
  }else{
    n = length(DependentVar)
    if(is.data.frame(DependentVar)){
      DependentVar = DependentVar[,1]
      n = length(DependentVar)
    }
  }

  # fix me: do I need a warning here is there is no samplesize in simulation
  if(!is.null(Regressor)){
    if (Intercept==TRUE && sum(as.matrix(Regressor)[,1])!=n ){
      Regressor = as.data.frame(cbind(rep(1,n),Regressor))
      names(Regressor)[1] = "Intercept"
    }
    else{
      Regressor = as.data.frame(Regressor)
    }
  }

  # number of regressors
  nreg = ifelse(is.null(Regressor), 0, DIM(Regressor)[2]-as.numeric(Intercept))

  if(is.null(Intercept) || is.null(Regressor)){
    nint = 1
  }else{
    nint = as.numeric(Intercept)
  }

  # get the indices of marginal parameters
  MargParmIndices = switch(CountDist,
                           "Poisson"               = 1:(1+nreg+nint-1),
                           "Negative Binomial"     = 1:(2+nreg+nint-1),
                           "Generalized Poisson"   = 1:(2+nreg+nint-1),
                           "Generalized Poisson 2" = 1:(2+nreg+nint-1),
                           "Binomial"              = 1:(1+nreg+nint-1),
                           "Mixed Poisson"         = 1:(3+(nreg+nint-1)*2),
                           "ZIP"                   = 1:(2+nreg+nint-1)
  )

  # retrieve the number of  parameters
  nMargParms = length(MargParmIndices)
  nparms     = nMargParms + sum(ARMAModel)
  nAR        = ARMAModel[1]
  nMA        = ARMAModel[2]

  if(Task %in% c("Synthesis", "Simulation")){
    if (nparms!=length(TrueParam)){
      stop("The length of the specified true parameter does not comply with the
           model or the number of regressors.")
    }
  }

  # cannot allow for sample size to be smaller than model parameters
  if (n<=(nparms+2))
    stop(sprintf("The provided sample size is equal to %.0f while the
                 the specified model has %.0f parameters",n,nparms))

  # check if initial param length is wrong
  if(!is.null(initialParam) & length(initialParam)!=nparms)
    stop("The specified initial parameter has wrong length.")

  # parse all information in the case without Regressors or in the case with Regressors
  if(nreg<1){
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"               = ppois,
                   "Negative Binomial"     = function(x, theta,...){   pnbinom(x, theta[1], 1-theta[2],...)},
                   "Generalized Poisson"   = function(x, theta,...){    pGpois(x, theta[1], theta[2],...)},
                   "Generalized Poisson 2" = function(x, theta,...,
                                                 lower.tail = TRUE,
                                                     log.p = FALSE){ pval = pgenpois2(x, theta[2], theta[1],lower.tail = lower.tail,...)
                                                                     if (log.p){
                                                                       pval = log(pval)
                                                                     }
                                                                       pval
                                                                       },
                   "Binomial"              = function(x, theta,...){    pbinom(x, ntrials, theta[1],...)},
                   "Mixed Poisson"         = function(x, theta,...){ pmixpois1(x, theta[1], theta[2], theta[3],...)},
                   "ZIP"                   = function(x, theta,...){          pzip(x, theta[1], theta[2],...)}
    )

    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"               = dpois,
                   "Negative Binomial"     = function(x, theta,...){   dnbinom(x, theta[1], 1-theta[2],...)},
                   "Generalized Poisson"   = function(x, theta,...){    dGpois(x, theta[1], theta[2],...)},
                   "Generalized Poisson 2" = function(x, theta,...){ dgenpois2(x, theta[2], theta[1],...)},
                   "Binomial"              = function(x, theta,...){    dbinom(x, ntrials, theta[1],...)},
                   "Mixed Poisson"         = function(x, theta,...){ dmixpois1(x, theta[1], theta[2], theta[3],...)},
                   "ZIP"                   = function(x, theta,...){          dzip(x, theta[1], theta[2],...)}
    )



    # retrieve marginal inverse cdf
    myinvcdf = switch(CountDist,
                      "Poisson"               = qpois,
                      "Negative Binomial"     = function(x, theta){   qnbinom(x, theta[1], 1-theta[2])},
                      "Generalized Poisson"   = function(x, theta){    qGpois(x, theta[1], theta[2])},
                      "Generalized Poisson 2" = function(x, theta){ qgenpois2(x, theta[2], theta[1])},
                      "Binomial"              = function(x, theta){    qbinom(x, ntrials, theta[1])},
                      "Mixed Poisson"         = function(x, theta){ qmixpois1(x, theta[1], theta[2], theta[3])},
                      "ZIP"                   = function(x, theta){      qzip(x, theta[1], theta[2])}
    )

    # lower bound constraints
    LB = switch(CountDist,
                "Poisson"               = c(0.01,              rep(-Inf, sum(ARMAModel))),
                "Negative Binomial"     = c(0.01, 0.01,        rep(-Inf, sum(ARMAModel))),
                "Generalized Poisson"   = c(0.01, 0.01,        rep(-Inf, sum(ARMAModel))),
                "Generalized Poisson 2" = c(0.01, 0.01,        rep(-Inf, sum(ARMAModel))),
                "Binomial"              = c(0.01,              rep(-Inf, sum(ARMAModel))),
                "Mixed Poisson"         = c(0.01, 0.01, 0.01,  rep(-Inf, sum(ARMAModel))),
                "ZIP"                   = c(0.01, 0.01,        rep(-Inf, sum(ARMAModel)))
    )
    # upper bound constraints
    UB = switch(CountDist,
                "Poisson"               = c(Inf,            rep( Inf, sum(ARMAModel))),
                "Negative Binomial"     = c(Inf, 0.99,      rep( Inf, sum(ARMAModel))),
                "Generalized Poisson"   = c(Inf, Inf,       rep( Inf, sum(ARMAModel))),
                "Generalized Poisson 2" = c(Inf, Inf,       rep( Inf, sum(ARMAModel))),
                "Binomial"              = c(Inf,            rep( Inf, sum(ARMAModel))),
                "Mixed Poisson"         = c(Inf, Inf, 0.99, rep( Inf, sum(ARMAModel))),
                "ZIP"                   = c(Inf, 0.99,      rep( Inf, sum(ARMAModel)))
    )
    # names of marginal parameters
    MargParmsNames = switch(CountDist,
                            "Poisson"               = c("lambda"),
                            "Negative Binomial"     = c("r","p"),
                            "Mixed Poisson"         = c("lambda_1", "lambda_2", "p"),
                            "Generalized Poisson"   = c("a", "mu"),
                            "Generalized Poisson 2" = c("a", "mu"),
                            "Binomial"              = c("p"),
                            "ZIP"                   = c("lambda", "p")
    )
  }else{
    # retrieve marginal cdf
    mycdf = switch(CountDist,
                   "Poisson"               = function(x, ConstMargParm, DynamMargParm,...){     ppois(x, DynamMargParm,...)},
                   "Negative Binomial"     = function(x, ConstMargParm, DynamMargParm,...){   pnbinom(x, ConstMargParm, 1-DynamMargParm,...)},
                   "Generalized Poisson"   = function(x, ConstMargParm, DynamMargParm,...){    pGpois(x, ConstMargParm, DynamMargParm,...)},
                   "Generalized Poisson 2" = function(x, ConstMargParm, DynamMargParm,...){ pgenpois2(x, DynamMargParm, ConstMargParm,...)},
                   "Binomial"              = function(x, ConstMargParm, DynamMargParm,...){    pbinom(x, ntrials, DynamMargParm,...)},
                   "Mixed Poisson"         = function(x, ConstMargParm, DynamMargParm,...){ pmixpois1(x, DynamMargParm[1], DynamMargParm[2], ConstMargParm,...)},
                   "ZIP"                   = function(x, ConstMargParm, DynamMargParm,...){          pzip(x, DynamMargParm, ConstMargParm,...)}
    )
    # retrieve marginal pdf
    mypdf = switch(CountDist,
                   "Poisson"               = function(x, ConstMargParm, DynamMargParm,...){     dpois(x, DynamMargParm,...)},
                   "Negative Binomial"     = function(x, ConstMargParm, DynamMargParm,...){   dnbinom(x, ConstMargParm, 1-DynamMargParm,...)},
                   "Generalized Poisson"   = function(x, ConstMargParm, DynamMargParm,...){    dGpois(x, ConstMargParm, DynamMargParm,...)},
                   "Generalized Poisson 2" = function(x, ConstMargParm, DynamMargParm,...){ dgenpois2(x, DynamMargParm, ConstMargParm,...)},
                   "Binomial"              = function(x, ConstMargParm, DynamMargParm,...){    dbinom(x, ntrials, DynamMargParm,...)},
                   "Mixed Poisson"         = function(x, ConstMargParm, DynamMargParm,...){ dmixpois1(x, DynamMargParm[1], DynamMargParm[2], ConstMargParm,...)},
                   "ZIP"                   = function(x, ConstMargParm, DynamMargParm,...){          dzip(x, DynamMargParm, ConstMargParm,...)}
    )
    # retrieve marginal inverse cdf
    myinvcdf = switch(CountDist,
                   "Poisson"               = function(x, ConstMargParm, DynamMargParm){     qpois(x, DynamMargParm)},
                   "Negative Binomial"     = function(x, ConstMargParm, DynamMargParm){   qnbinom(x, ConstMargParm, 1-DynamMargParm)},
                   "Generalized Poisson"   = function(x, ConstMargParm, DynamMargParm){    qGpois(x, ConstMargParm, DynamMargParm)},
                   "Generalized Poisson 2" = function(x, ConstMargParm, DynamMargParm){ qgenpois2(x, DynamMargParm, ConstMargParm)},
                   "Binomial"              = function(x, ConstMargParm, DynamMargParm){    qbinom(x, ntrials, DynamMargParm)},
                   "Mixed Poisson"         = function(x, ConstMargParm, DynamMargParm){ qmixpois1(x, DynamMargParm[,1], DynamMargParm[,2], ConstMargParm)},
                   "ZIP"                   = function(x, ConstMargParm, DynamMargParm){      qzip(x, DynamMargParm, ConstMargParm)}
    )
    # lower bound contraints
    LB = switch(CountDist,
                "Poisson"               = rep(-Inf, sum(ARMAModel)+nreg+nint),
                "Negative Binomial"     = c(rep(-Inf, nreg+nint), 0.001, rep(-Inf, sum(ARMAModel))),
                "Generalized Poisson"   = c(rep(-Inf, nreg+nint), 0.001, rep(-Inf, sum(ARMAModel))),
                "Generalized Poisson 2" = c(rep(-Inf, nreg+nint), 0.001, rep(-Inf, sum(ARMAModel))),
                "Binomial"              = rep(-Inf, sum(ARMAModel)+nreg+nint),
                "Mixed Poisson"         = c(rep(-Inf, 2*(nreg+nint)), 0.001, rep(-Inf, sum(ARMAModel))),
                "ZIP"                   = c(rep(-Inf, nreg+nint), 0.001, rep(-Inf, sum(ARMAModel)))
    )
    # upper bound constraints
    UB = switch(CountDist,
                "Poisson"               = rep(Inf, sum(ARMAModel)+nreg+nint),
                "Negative Binomial"     = c(rep(Inf, nreg+nint), Inf, rep(Inf, sum(ARMAModel))),
                "Generalized Poisson"   = c(rep(Inf, nreg+nint), Inf, rep(Inf, sum(ARMAModel))),
                "Generalized Poisson 2" = c(rep(Inf, nreg+nint), Inf, rep(Inf, sum(ARMAModel))),
                "Binomial"              = rep(Inf, sum(ARMAModel)+nreg+nint),
                "Mixed Poisson"         = c(rep(Inf, 2*(nreg+nint)), 0.49, rep(Inf, sum(ARMAModel))),
                "ZIP"                   = c(rep(Inf, nreg+nint), Inf, rep(Inf, sum(ARMAModel))),
    )
    # retrieve names of marginal parameters
    MargParmsNames = switch(CountDist,
                            "Poisson"               = paste(rep("b_",nreg),(1-nint):nreg,sep=""),
                            "Negative Binomial"     = c(paste(rep("b_",nreg),(1-nint):nreg,sep=""), "k"),
                            "Mixed Poisson"         = c(paste(rep("b_1",nreg),(1-nint):nreg,sep=""),paste(rep("b_2",nreg),(1-nint):nreg,sep=""), "p"),
                            "Generalized Poisson"   = c(paste(rep("b_",nreg),(1-nint):nreg,sep=""), "a"),
                            "Generalized Poisson 2" = c(paste(rep("b_",nreg),(1-nint):nreg,sep=""), "a"),
                            "Binomial"              = paste(rep("b_",nreg),(1-nint):nreg,sep=""),
                            "ZIP"                   = c(paste(rep("b_",nreg),(1-nint):nreg,sep=""), "p"),
    )
  }

  # check whether the provided initial estimates make sense - this is for Evaluation, Simulation, Optimization Task
  if (!is.null(initialParam)) {
    if (max(initialParam<=LB) | max(initialParam>=UB)){
      stop("The specified initial parameter is outside the feasible region.")
    }
  }

  # check whether the provided initial estimates make sense - This is for Synthesis Tasks
  if (Task=="Synthesis") {
    if (max(TrueParam<=LB) | max(TrueParam>=UB)){
      stop("The specified parameter is outside the feasible region.")
    }
  }

  # create names of the ARMA parameters
  if(nAR>0) ARNames = paste("AR_",1:ARMAModel[1], sep="")
  if(nMA>0) MANames = paste("MA_",1:ARMAModel[2], sep="")

  # put all the names together
  if(nAR>0 && nMA<1) parmnames = c(MargParmsNames, ARNames)
  if(nAR<1 && nMA>0) parmnames = c(MargParmsNames, MANames)
  if(nAR>0 && nMA>0) parmnames = c(MargParmsNames, ARNames, MANames)
  if(nAR==0 && nMA==0) parmnames = c(MargParmsNames)

  # add the parmnames on theta fix me: does this affect performance?
  if(!is.null(initialParam)) names(initialParam) = parmnames

  # value I will set the loglik when things go bad (e.g. non invertible ARMA)
  loglik_BadValue1 = 10^8

  # value I will set the loglik when things go bad (e.g. non invertible ARMA)
  loglik_BadValue2 = 10^9

  out = list(
    mycdf            = mycdf,
    mypdf            = mypdf,
    myinvcdf         = myinvcdf,
    MargParmIndices  = MargParmIndices,
    initialParam     = initialParam,
    TrueParam        = TrueParam,
    nsim             = nsim,
    parmnames        = parmnames,
    nMargParms       = nMargParms,
    nAR              = nAR,
    nMA              = nMA,
    n                = n,
    nreg             = nreg,
    nint             = nint,
    CountDist        = CountDist,
    ARMAModel        = ARMAModel,
    ParticleNumber   = ParticleNumber,
    epsilon          = epsilon,
    nparms           = nparms,
    UB               = UB,
    LB               = LB,
    EstMethod        = EstMethod,
    DependentVar     = DependentVar,
    Regressor        = Regressor,
    Task             = Task,
    OptMethod        = OptMethod,
    OutputType       = OutputType,
    ParamScheme      = ParamScheme,
    maxdiff          = maxdiff,
    loglik_BadValue1 = loglik_BadValue1,
    loglik_BadValue2 = loglik_BadValue2,
    ntrials          = ntrials,
    Intercept        = Intercept,
    verbose          = verbose
    )
  return(out)

}


#' Particle Filter likelihood function with Resampling for ARMA Models
#'
#' Uses particle filtering with resampling to approximate the likelihood of
#' a specified count time series model with an underlying ARMA dependence structure.
#'
#' @param theta Numeric vector. Parameter vector for the model.
#' @param mod List. A list containing all model specifications, including:
#'   \describe{
#'     \item{DependentVar}{Numeric. The dependent time series variable to be modeled.}
#'     \item{Task}{Character. The task requested by the user (e.g., Evaluation, Optimization, Synthesis, etc.).}
#'     \item{ParticleNumber}{Integer. Number of particles to use.}
#'     \item{Regressor}{Optional independent variable(s).}
#'     \item{CountDist}{Character string. Specifies the count marginal distribution.}
#'   }
#'
#' @return Numeric. Approximate log-likelihood of the model.
#'
#' @examples
#' # Example: Compute the log-likelihood of a Poisson-ARMA(1,0) model
#' CountDist      <- "Poisson"
#' n              <- 100
#' ARMAModel      <- c(1,0)
#' MargParm       <- 3
#' ARParm         <- 0.5
#' MAParm         <- NULL
#' initialParam   <- c(3, ARParm, MAParm)
#'
#' # Simulate data
#' set.seed(2)
#' DependentVar <- sim_lgc(n, CountDist, MargParm, ARParm, MAParm)
#'
#' # Build model specification
#' mod <- ModelScheme(
#'   DependentVar   = DependentVar,
#'   CountDist      = CountDist,
#'   ARMAModel      = ARMAModel,
#'   initialParam   = initialParam
#' )
#'
#' # Compute log-likelihood at a given parameter point
#' ll <- ParticleFilter_Res_ARMA(initialParam, mod)
#' ll
#' @export
ParticleFilter_Res_ARMA = function(theta, mod){

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # Retrieve parameters and save them in a list called Parms
  Parms = RetrieveParameters(theta,mod)

  # check for causality
  if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))

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
  #Z[1,]    = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
  Z[1,]    = SampleTruncNormParticles(mod, Limit, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))


  for (t in 2:mod$n){

    # compute Zhat_t
    #Zhat[t,] = ComputeZhat_t(m,Theta,Z,Zhat,t, Parms,p,q, nTheta, Theta_n)
    Zhat[t,] = ComputeZhat_t(mod, IA, Z, Zhat,t, Parms)

    # Compute integral limits
    Limit = ComputeLimits(mod, Parms, t, Zhat[t,], Rt)

    # Sample truncated normal particles
    #Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b, t, Zhat[t,], Rt)
    Znew  = SampleTruncNormParticles(mod, Limit, t, Zhat[t,], Rt)

    # update weights
    #w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])
    w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      #print(t)
      #print(w[t,])
      message(sprintf('WARNING: At t=%.0f some of the weights are either too small or sum to 0',t))
      return(10^8)
    }
    # Resample the particles using common random numbers
    old_state1 = get_rand_state()
    Znew = ResampleParticles(mod, w, t, Znew)
    set_rand_state(old_state1)

    # save the current particle
    Z[t,]   = Znew

    #print(t)
    #print(nloglik)
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


#' Optimization wrapper to fit PF likelihood with resampling (new version)
#'
#' Fits the particle filter log-likelihood using resampling. This version expects
#' a complete model object (created by \code{\link{ModelScheme}}) and performs
#' one optimization per value of \code{mod$ParticleNumber}.
#'
#' @param theta Numeric vector. Initial parameter values.
#' @param mod A list containing model-related metadata and control settings, typically
#'   returned by \code{\link{ModelScheme}}.
#'
#' @return A list containing parameter estimates, standard errors, log-likelihood,
#'   AIC, BIC, AICc, and optimization diagnostics.
#'
#' @examples
#' # Specify model characteristics
#' CountDist      <- "Poisson"
#' n              <- 100
#' ARMAModel      <- c(1, 0)
#' MargParm       <- 3
#' ARParm         <- 0.5
#' MAParm         <- NULL
#' initialParam   <- c(3.2, 0.9 * ARParm, 0.8 * MAParm)
#' Task           <- "Optimization"
#'
#' # Simulate data
#' set.seed(2)
#' DependentVar <- sim_lgc(n, CountDist, MargParm, ARParm, MAParm)
#'
#' # Prepare the model object
#' mod <- ModelScheme(
#'   DependentVar   = DependentVar,
#'   CountDist      = CountDist,
#'   ARMAModel      = ARMAModel,
#'   initialParam   = initialParam,
#'   Task           = Task
#' )
#'
#' # Fit the model
#' fit <- FitMultiplePF_Res(initialParam, mod)
#' fit
#'
#' @export
FitMultiplePF_Res = function(theta, mod){

  # retrieve parameter, sample size etc
  nparts    = length(mod$ParticleNumber)
  nparms    = length(theta)
  nfit      = 1
  n         = length(mod$DependentVar)

  # allocate memory to save parameter estimates, hessian values, and loglik values
  ParmEst   = matrix(0,nrow=nfit*nparts,ncol=nparms)
  se        = matrix(NA,nrow=nfit*nparts,ncol=nparms)
  loglik    = rep(0,nfit*nparts)
  convcode  = rep(0,nfit*nparts)
  kkt1      = rep(0,nfit*nparts)
  kkt2      = rep(0,nfit*nparts)

  Criteria  = matrix(0,nrow=nfit*nparts,ncol=3)

  # Each realization will be fitted nfit*nparts many times
  for (j in 1:nfit){
    set.seed(j)
    # for each fit repeat for different number of particles
    for (k in 1:nparts){
      # FIX ME: I need to somehow update this in mod. (number of particles to be used). I t now works only for 1 choice of ParticleNumber
      # ParticleNumber = mod$ParticleNumber[k]

      # Update the ParticleNumber inside mod for current iteration
      mod_current                = mod
      mod_current$ParticleNumber = mod$ParticleNumber[k]



      if(mod_current$Task == 'Optimization'){
        # run optimization for our model --no ARMA model allowed
        optim.output <- optimx(
          par     = theta,
          fn      = ParticleFilter,
          lower   = mod_current$LB,
          upper   = mod_current$UB,
          hessian = TRUE,
          method  = mod_current$OptMethod,
          mod     = mod_current)
        loglikelihood = optim.output["value"]
      }else{
        optim.output = as.data.frame(matrix(rep(NA,8+length(theta)), ncol=8+length(theta)))
        names(optim.output) = c(mapply(sprintf, rep("p%.f",length(theta)), (1:length(theta)), USE.NAMES = FALSE),
                                "value",  "fevals", "gevals", "niter", "convcode",  "kkt1", "kkt2", "xtime")

        optim.output[,1:length(theta)] = theta
        startTime = Sys.time()
        loglikelihood = ParticleFilter(theta,mod_current)
        endTime = Sys.time()
        runTime = difftime(endTime, startTime, units = 'secs')
        optim.output[,(length(theta)+1)] = loglikelihood
        optim.output[,(length(theta)+2)] = 1
        optim.output[,(length(theta)+3)] = 1
        optim.output[,(length(theta)+4)] = 0
        optim.output[,(length(theta)+8)] = as.numeric(runTime)
      }


      # save estimates, loglik value and diagonal hessian
      ParmEst[nfit*(k-1)+j,]  = as.numeric(optim.output[1:nparms])
      loglik[nfit*(k-1) +j]   = optim.output$value
      convcode[nfit*(k-1) +j] = optim.output$convcode
      kkt1[nfit*(k-1) +j]     = optim.output$kkt1
      kkt2[nfit*(k-1) +j]     = optim.output$kkt2


      # compute Hessian
      # t5 = tic()
      if(mod_current$Task == "Optimization"){
        H = gHgen(par          = ParmEst[nfit*(k-1)+j,],
                  fn           = ParticleFilter,
                  mod          = mod_current)
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
            se[nfit*(k-1)+j,NAIndex] = rep(NA, length(NAIndex))
          }
        }
      }
      # t6 = tic()-t6
      # print(t6)
      # compute model selection criteria
      Criteria[nfit*(k-1)+j,]  =  Criteria.lgc(loglik[nfit*(k-1) +j], mod_current)
    }
  }

  # commenting this out to add it in the loop and allow for multiple fits
  # Compute model selection criteria (assuming one fit)
  # Criteria = Criteria.lgc(loglik, mod_current)

  FitResults             = list()
  FitResults$ParmEst     = ParmEst
  FitResults$se          = se
  FitResults$FitStatistics  = cbind(loglik, Criteria)
  FitResults$OptimOutput = cbind(convcode,kkt1,kkt2)

  colnames(FitResults$FitStatistics)     = c("loglik", "AIC", "BIC", "AICc")
  if(mod_current$Task=="Optimization")colnames(FitResults$OptimOutput)       = c("ConvergeStatus", "kkt1", "kkt2")

  return(FitResults)
}


#' Prepare Output from Model Estimation Wrapper
#'
#' Formats the output of a model fitting or evaluation task. Depending on the
#' value of \code{mod$OutputType}, the function returns either a structured
#' list (default) or a wide-format data frame containing parameter estimates,
#' standard errors, model configuration, and fit statistics.
#'
#' @param mod A list of model specifications and settings, typically created using
#'   \code{\link{ModelScheme}}. Contains all inputs used in estimation or simulation,
#'   such as the count distribution, ARMA structure, optimization method, and task type.
#' @param FitResults A list containing outputs from the model fitting process.
#'   Expected components include:
#'   \itemize{
#'     \item \code{ParmEst}: A matrix of parameter estimates.
#'     \item \code{se}: A matrix of standard errors (if applicable).
#'     \item \code{FitStatistics}: A named vector of fit measures (e.g., log-likelihood, AIC).
#'     \item \code{OptimOutput}: A list or vector of optimizer diagnostics (e.g., convergence code, KKT conditions).
#'   }
#'
#' @return If \code{mod$OutputType == "list"}, returns a named list with:
#' \itemize{
#'   \item \code{ParamEstimates}: Parameter estimates matrix.
#'   \item \code{StdErrors}: Standard error matrix (if \code{mod$Task == "Optimization"}).
#'   \item \code{FitStatistics}: Log-likelihood, AIC, BIC, AICc, etc.
#'   \item \code{OptimOutput}: Optimizer output (if applicable).
#'   \item \code{Model}: Description of the model (e.g., distribution + ARMA order).
#'   \item \code{Task}, \code{EstMethod}, \code{SampleSize}, etc.
#'   \item \code{WarnMessage}: A message indicating specific issues during fitting (if any).
#' }
#'
#' If \code{mod$OutputType != "list"}, returns a single-row data frame where each column
#' represents a parameter, fit statistic, or model metadata field.
#'
#' @details
#' This function is designed to be called within a wrapper or high-level function that
#' runs particle filtering, optimization, or simulation. It ensures consistent formatting
#' and labeling of outputs for downstream use (e.g., summaries, diagnostics, plots).
#'
#' @note This function is intended for internal use by the package and is not typically called directly by users.
#' @keywords internal
#' @export
PrepareOutput = function(mod, FitResults){


  # if(mod$OutputType=="list"){
  if(mod$Task %in% c("Evaluation","Optimization")){
    #  save the results in a list
    ModelOutput = list()

    # populate information from fitting
    ModelOutput$ParamEstimates = FitResults$ParmEst
    ModelOutput$FitStatistics     = FitResults$FitStatistics

    if(mod$Task=="Optimization"){
      ModelOutput$StdErrors    = FitResults$se
      ModelOutput$OptimOutput = FitResults$OptimOutput
    }

    # populate information from parsing the initial inputs
    ModelOutput$Model          = paste(mod$CountDist,  paste("ARMA(",mod$ARMAModel[1],",",mod$ARMAModel[2],")",sep=""), sep= "-")
    ModelOutput$Task           = mod$Task
    if(!is.null(mod$TrueParam)) ModelOutput$TrueParam       = mod$TrueParam
    if(!is.null(mod$initialParam)) ModelOutput$initialParam = mod$initialParam
    ModelOutput$EstMethod      = mod$EstMethod
    ModelOutput$SampleSize     = mod$n


    if(FitResults$FitStatistics[,"loglik"]==mod$loglik_BadValue1) ModelOutput$WarnMessage = "WARNING: The ARMA polynomial must be causal and invertible."
    if(FitResults$FitStatistics[,"loglik"]==mod$loglik_BadValue2) ModelOutput$WarnMessage = "WARNING: Some of the weights are either too small or sum to 0."
    # assign names to all output elements
    if(!is.null(mod$TrueParam)) names(ModelOutput$TrueParam)      = mod$parmnames
    colnames(ModelOutput$ParamEstimates) = mod$parmnames

    if(mod$Task=="Optimization")colnames(ModelOutput$StdErrors)      = paste("se(", mod$parmnames,")", sep="")
    #names(ModelOutput$FitStatistics)     = c("loglik", "AIC", "BIC", "AICc")
    #if(mod$Task=="Optimization")names(ModelOutput$OptimOutput)       = c("ConvergeStatus", "kkt1", "kkt2")

    ModelOutput$mod            = mod
  }else if(mod$Task %in% c("Simulation")){
    ModelOutput  = data.frame(matrix(ncol = 4*mod$nparms+16, nrow = mod$nsim))

    # specify output names
    if(!is.null(mod$TrueParam)){
      colnames(ModelOutput) = c(
        'CountDist','ARMAModel', 'Regressor',
        paste("True_", mod$parmnames, sep=""), paste("InitEst_", mod$parmnames, sep=""),
        mod$parmnames, paste("se(", mod$parmnames,")", sep=""),
        'EstMethod', 'SampleSize', 'ParticleNumber', 'epsilon', 'OptMethod', 'ParamScheme',
        "loglik", "AIC", "BIC", "AICc", "ConvergeStatus", "kkt1", "kkt2")
    }
    # if mod$Task = Simulation then I must have TrueParam
    # colnames(ModelOutput) = c(
    #   'CountDist','ARMAModel', 'Regressor',
    #   paste("InitialEstim_", mod$parmnames, sep=""),
    #   mod$parmnames, paste("se(", mod$parmnames,")", sep=""),
    #   'EstMethod', 'SampleSize', 'ParticleNumber', 'epsilon', 'OptMethod',
    #   "loglik", "AIC", "BIC", "AICc", "ConvergeStatus", "kkt1", "kkt2")

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

    LastOffset = offset
    for (i in seq_along(FitResults)) {

      fit_i <- FitResults[[i]]
      offset = LastOffset

      # Initial Parameter Estimates
      if(mod$nMargParms>0){
        ModelOutput[i, offset:(offset + mod$nMargParms -1 )] = fit_i$initialParam[1:mod$nMargParms ]
        offset = offset + mod$nMargParms
      }
      if(mod$nAR>0){
        ModelOutput[i, offset:(offset + mod$nAR-1)]        = fit_i$initialParam[(mod$nMargParms+1):(mod$nMargParms+mod$nAR) ]
        offset = offset + mod$nAR
      }

      if(mod$nMA>0){
        ModelOutput[i, offset:(offset + mod$nMA-1)]        = fit_i$initialParam[(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA) ]
        offset = offset + mod$nMA
      }

      # Parameter Estimates
      if(mod$nMargParms>0){
        ModelOutput[i, offset:(offset + mod$nMargParms-1)] = fit_i$ParmEst[,1:mod$nMargParms]
        offset = offset + mod$nMargParms
      }

      if(mod$nAR>0){
        ModelOutput[i, offset:(offset + mod$nAR-1)]        = fit_i$ParmEst[,(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]
        offset = offset + mod$nAR
      }

      if(mod$nMA>0){
        ModelOutput[i, offset:(offset + mod$nMA-1)]        = fit_i$ParmEst[,(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA)]
        offset = offset + mod$nMA
      }

      # Parameter Std Errors
      if(mod$nMargParms>0){
        ModelOutput[i, offset:(offset + mod$nMargParms-1)] = fit_i$se[,1:mod$nMargParms]
        offset = offset + mod$nMargParms
      }

      if(mod$nAR>0){
        ModelOutput[i, offset:(offset + mod$nAR-1)]        = fit_i$se[,(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]
        offset = offset + mod$nAR
      }

      if(mod$nMA>0){
        ModelOutput[i, offset:(offset + mod$nMA-1)]        = fit_i$se[,(mod$nMargParms+mod$nAR+1):(mod$nMargParms+mod$nAR+mod$nMA)]
      }

      ModelOutput[i,]$loglik         = fit_i$FitStatistics[,"loglik"]
      ModelOutput[i,]$AIC            = fit_i$FitStatistics[,"AIC"]
      ModelOutput[i,]$BIC            = fit_i$FitStatistics[,"BIC"]
      ModelOutput[i,]$AICc           = fit_i$FitStatistics[,"AICc"]
      ModelOutput[i,]$ConvergeStatus = fit_i$OptimOutput[,"ConvergeStatus"]
      ModelOutput[i,]$kkt1           = fit_i$OptimOutput[,"kkt1"]
      ModelOutput[i,]$kkt2           = fit_i$OptimOutput[,"kkt2"]

    }
    ModelOutput$EstMethod      = mod$EstMethod
    ModelOutput$SampleSize     = mod$n
    ModelOutput$ParticleNumber = mod$ParticleNumber
    ModelOutput$epsilon        = mod$epsilon
    ModelOutput$OptMethod      = mod$OptMethod
    if(!is.null(mod$TrueParam)) ModelOutput$ParamScheme    = mod$ParamScheme
  }

  return(ModelOutput)

}

#' Get current random seed state
#' Add some functions that I will need for common random numbers
#'
#' @keywords internal
get_rand_state <- function() {
  # Using `get0()` here to have `NULL` output in case object doesn't exist.
  # Also using `inherits = FALSE` to get value exactly from global environment
  # and not from one of its parent.
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

#' Set random seed state
#'
#' @keywords internal
set_rand_state <- function(state) {
  # Assigning `NULL` state might lead to unwanted consequences
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
  }
}


#' Check Causality and Invertibility of ARMA Model
#'
#' Checks whether the specified ARMA model is stable — i.e., whether the AR polynomial
#' is causal and the MA polynomial is invertible. This is done by examining the
#' location of the roots of the characteristic polynomials.
#'
#' @param AR Numeric vector. Autoregressive (AR) coefficients of the model. Can be \code{NULL} if the model has no AR part.
#' @param MA Numeric vector. Moving average (MA) coefficients of the model. Can be \code{NULL} if the model has no MA part.
#'
#' @return Integer. Returns \code{1} if the model is **not stable** (i.e., at least one root of the AR or MA polynomial
#'          lies inside the unit circle), and \code{0} if the model is stable.
#'
#' @details
#' The function uses \code{\link[base]{polyroot}} to compute the roots of the AR and MA characteristic polynomials:
#' \itemize{
#'   \item The AR polynomial is defined as \eqn{1 - AR_1 z - AR_2 z^2 - \dots}.
#'   \item The MA polynomial is defined similarly.
#' }
#' If any root lies **within** the unit circle (i.e., has modulus < 1), the model is considered unstable.
#'
#' @examples
#' # Stable ARMA(1,1)
#' CheckStability(AR = 0.5, MA = 0.3)  # returns 0
#'
#' # Unstable ARMA(1,1)
#' CheckStability(AR = 1.2, MA = 0.3)  # returns 1
#'
#' # AR-only model
#' CheckStability(AR = c(0.9), MA = NULL)
#'
#' # MA-only model
#' CheckStability(AR = NULL, MA = c(1.5))
#'
#' @export
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


#' Compute Initial Parameter Estimates
#'
#' Generates initial parameter estimates based on the model specifications provided in \code{mod}.
#' These starting values are typically used as input to optimization routines or simulations in
#' count time series models with ARMA dependencies and various marginal count distributions.
#'
#' @param mod List. A list containing all model specifications, including:
#'   \describe{
#'     \item{DependentVar}{Numeric. The dependent time series variable to be modeled.}
#'     \item{Task}{Character. The task requested by the user (e.g., Evaluation, Optimization, Synthesis, etc.).}
#'     \item{ParticleNumber}{Integer. Number of particles to use.}
#'     \item{Regressor}{Optional independent variable(s).}
#'     \item{CountDist}{Character string. Specifies the count marginal distribution.}
#'   }
#'
#' @return A numeric vector of initial parameter estimates consistent with the model structure
#'   described in \code{mod}. The vector includes:
#'   \itemize{
#'     \item Marginal distribution parameters
#'     \item Regression coefficients (if any)
#'     \item AR and MA coefficients
#'   }
#'   The names of the vector elements (if assigned) match those in \code{mod$parmnames}.
#'
#' @details
#' The initial estimates are heuristically determined and may depend on the distribution
#' type, sample statistics (e.g., mean, variance), or defaults chosen for stability.
#' These estimates are used as starting points in numerical optimization procedures.
#' GLM and MoM estimates are used for marginal parameters and Yulew-Walker for ARMA.
#'
#' @examples
#' mod <- ModelScheme(
#'   DependentVar = rpois(100, lambda = 3),
#'   CountDist = "Poisson",
#'   ARMAModel = c(1, 1),
#'   Intercept = TRUE
#' )
#' init <- InitialEstimates(mod)
#' print(init)
#'
#' @export
InitialEstimates = function(mod){
  # require(itsmr)

  # allocate memory
  est  = rep(NA, mod$nMargParms+sum(mod$ARMAModel))

  #------------First compute estimates of marginal parameters
  #-----Poisson case
  if(mod$CountDist=="Poisson"){
    if(mod$nreg==0){
      # Method of moments estimate for Poisson parameter (via the mean)
      est[1] = mean(mod$DependentVar)
    }else{
      # GLM for the mean that depends on time

      # CHECK ME: If I fit a Poisson AR(3) in the the data example of the JASA paper, but the code below doesn't
      # specify poisson family (it would pick up the default distribution that glm function has) then there will
      # be a numerical error in the likelihood. Check it!

      # to fix #36 I will introduce an ifelse below. Another way to do things would be to pass the dataframe
      # through the mod structure and use the data frame and the variable names in classic lm fashion. However,
      # with my design of passing the DependentVar and the Regrerssor as separate arguments in mod, I need the
      # ifelse statement below.

      # glmPoisson            = glm(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)], family = "poisson")
      if(mod$nint){
        glmPoisson            = glm(mod$DependentVar~ as.matrix(mod$Regressor[,2:(mod$nreg+1)]), family = "poisson")
      }else{
        glmPoisson            = glm(mod$DependentVar~0 + as.matrix(mod$Regressor[,1:mod$nreg]), family = "poisson")
      }

      est[1:mod$nMargParms] = as.numeric(glmPoisson[1]$coef)

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

    }else{
      # GLM for the mean that depends on time
      #glmNB                     = glm.nb(mod$DependentVar~mod$Regressor[,2:(mod$nreg+1)])

      if(mod$nint){
        glmNB            = glm.nb(mod$DependentVar ~ as.matrix(mod$Regressor[,2:(mod$nreg+1)]))
      }else{
        glmNB            = glm.nb(mod$DependentVar ~ 0 + as.matrix(mod$Regressor[,1:mod$nreg]))
      }

      # check me: I think that glmNegBin is a typo below and I should be using glmNB
      # est[1:(mod$nMargParms-1)] = as.numeric(glmNegBin[1]$coef)
      est[1:(mod$nMargParms-1)] = as.numeric(glmNB[1]$coef)

      # MoM for the over dispersion param in NegBin2 parametrization
      est[mod$nMargParms]       = sum(glmNB$fitted.values^2)/(sum((mod$DependentVar-glmNB$fitted.values)^2-glmNB$fitted.values))

    }
  }

  if(mod$CountDist=="Mixed Poisson"){
    if(mod$nreg==0){
      # pmle for marginal parameters
      MixPois_PMLE <- pmle.pois(as.numeric(mod$DependentVar),2)

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
    }else{
      #library(mixtools)
      # MP_fit = poisregmixEM(mod$DependentVar, mod$Regressor[,2:(mod$nreg+1)])

      if(mod$nint){
        MP_fit            = poisregmixEM(mod$DependentVar, as.matrix(mod$Regressor[,2:(mod$nreg+1)]))
      }else{
        MP_fit            = poisregmixEM(mod$DependentVar, as.matrix(mod$Regressor[,1:mod$nreg]),addintercept = FALSE)
      }

      if(MP_fit$lambda[1]<0.5){
        est[1:mod$nMargParms] = as.numeric(c(MP_fit$beta, MP_fit$lambda[1]))
      }
      else{
        # check me: should I give an error here?
        # check me: does this work for more regressors? I think there will be an error with the indices below
        if(mod$nint){
          est[1:mod$nMargParms] = as.numeric(c(MP_fit$beta[3:4],MP_fit$beta[1:2], 1-MP_fit$lambda[1]))
        }else{
          est[1:mod$nMargParms] = as.numeric(c(MP_fit$beta[2],MP_fit$beta[1], 1-MP_fit$lambda[1]))
        }
        }
      #could also use the flexmix package
      #check me: for now we allow mixture of only two components
      #MP_fit <- flexmix(mod$DependentVar ~ mod$Regressor[,2:(mod$nreg+1)], k = 2, model = FLXMRglm(family = "poisson"))
      #refit1 = refit(MP_fit)

      #beta1Hat = as.numeric(refit1@coef[1:2])
      #beta2Hat = as.numeric(refit1@coef[3:4])
      #ProbHat  = MP_fit@prior[1]

      #est[1:mod$nMargParms] = c(beta1Hat, beta2Hat,ProbHat)

    }
  }

  #-----Generalized Poisson case
  if(mod$CountDist=="Generalized Poisson" || mod$CountDist=="Generalized Poisson 2"){
    if(mod$nreg==0){
      xbar    = mean(mod$DependentVar)
      sSquare = var(mod$DependentVar)


      # the GenPois density is parametrized through the mean with the pair (alpha,mu) - see (2.4) in Famoye 1994
      # solve (2.5) in Famoye 1994 wrt alpha and replace sample estimates
      alpha_1 = (sqrt(sSquare/xbar) - 1)/xbar
      alpha_2 = (-sqrt(sSquare/xbar) - 1)/xbar

      # to choose between the two solutions I ll check "overdispersion"
      if (sSquare>=xbar) alpha=alpha_1
      if (sSquare<xbar) alpha =alpha_2

      # FIX ME: the GenPois densities in VGAM do not allow for negative alpha yet
      if (alpha<0) alpha = 10^(-6)
      est[1:2] = c(alpha, xbar)

    }else{
      # run GenPois GLM using VGAM package - CHECK ME: shouls surface maxit as an option to the user?
      if(mod$nint){
        fit = VGAM::vglm( mod$DependentVar ~ as.matrix(mod$Regressor[,2:(1+mod$nreg)]), genpoisson2, maxit=60)
      }else{
        fit = VGAM::vglm( mod$DependentVar ~  0 + as.matrix(mod$Regressor[,1:mod$nreg]), genpoisson2(zero=0), maxit=60,)
      }
      # save linear predictor coefficients
      est[1:(mod$nMargParms-1)]  = as.numeric(coef(fit, matrix = TRUE)[,1])

      # save dispersion coefficient
      est[mod$nMargParms] = loglink(as.numeric(coef(fit, matrix = TRUE)[1,2]),inverse = TRUE)
    }
  }

  #-----Binomial case
  if(mod$CountDist=="Binomial"){
    if(mod$nreg==0){
      xbar = mean(mod$DependentVar)
      # Method of Moments for Binomial E(X)=rp, where  r = ntrials
      pEst = xbar/mod$ntrials
      est[1] = pEst
    }else{
      #glmBinomial               = glm(cbind(mod$DependentVar,mod$ntrials-mod$DependentVar) ~ mod$Regressor[,2:(mod$nreg+1)] , family = 'binomial')
      if(mod$nint){
        glmBinomial = glm(cbind(mod$DependentVar,mod$ntrials-mod$DependentVar) ~ as.matrix(mod$Regressor[,2:(mod$nreg+1)]) , family = 'binomial')
      }else{
        glmBinomial = glm(cbind(mod$DependentVar,mod$ntrials-mod$DependentVar) ~  0 + as.matrix(mod$Regressor[,1:mod$nreg]) , family = 'binomial')

      }
      est[1:mod$nMargParms] = as.numeric(glmBinomial$coef)
    }
  }

  #-----Zero Inflation Poisson
  if(mod$CountDist=="ZIP"){
    if(mod$nreg==0){
      # pmle for marginal parameters
      ZIP_PMLE <- poisson.zihmle(mod$DependentVar, type = c("zi"), lowerbound = 0.01, upperbound = 10000)

      lEst = ZIP_PMLE[1]
      pEst = ZIP_PMLE[2]
      est[1:2] = c(lEst, pEst)

    }else{
      #zeroinfl_reg <- zeroinfl( mod$DependentVar~ mod$Regressor[,2:(mod$nreg+1)] | 1,)
      #zeroinfl_reg = vglm( mod$DependentVar~ mod$Regressor[,2:(mod$nreg+1)], zipoisson(zero=1))

      if(mod$nint){
        zeroinfl_reg = vglm( mod$DependentVar ~ as.matrix(mod$Regressor[,2:(mod$nreg+1)]), zipoisson(zero=1))
      }else{
        zeroinfl_reg = vglm( mod$DependentVar ~ 0 + as.matrix(mod$Regressor[,1:mod$nreg]), zipoisson(zero=0))
      }

      est[1:(mod$nMargParms-1)] = as.numeric(coef(zeroinfl_reg))[2:(mod$nMargParms)]
      est[mod$nMargParms] = plogis(as.numeric(coef(zeroinfl_reg))[1])
    }
  }

  #------------ARMA Initial Estimates
  # Transform (1) in the JASA paper to retrieve the "observed" latent series and fit an ARMA
  # check me: Oct 2025 - the following needs to be refactored
    if(mod$nreg==0){

    Cxt = mod$mycdf(mod$DependentVar,est[1:mod$nMargParms])
    # if Cxt = 1, I will need to deal with this somehow. In the Binomial case it seems very likely to get 1
    # but this can also happen for other distributions. So far, in the initial estimation I have only seen the issue
    # taking place in the Binomial case, however, I have come across this issue on the Particle filter estimation
    # for misspecified models. For now I will simply set
    if (mod$CountDist=="Binomial"){
      Cxt[Cxt==1] = 1-10^(-16)
      Cxt[Cxt==0] = 0+10^(-16)
    }
    # I am doing it for the binomial only so that I will get an error if it happens for another case.
    if(max(mod$ARMAModel)>0) armafit = itsmr::arma(qnorm(Cxt),mod$nAR,mod$nMA)
    if(mod$nAR>0) est[(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]                    = armafit$phi
    if(mod$nMA>0) est[(1+mod$nMargParms+mod$nAR):(mod$nMargParms+sum(mod$ARMAModel))] = armafit$theta
  }else{
    # put the parameters in appropriate format
    Params = RetrieveParameters(est,mod)

    # see comment above regarding Cxt=1 and Cxt = 0
    # while working on #27 and testing misspecified models I came across the same issue for Poisson.
    Cxt = mod$mycdf(mod$DependentVar,Params$ConstMargParm, Params$DynamMargParm)
    if (mod$CountDist=="Binomial" || mod$CountDist=="Mixed Poisson" || mod$CountDist=="Poisson"){
      Cxt[Cxt==1] = 1-10^(-16)
      Cxt[Cxt==0] = 0+10^(-16)
    }
    # Transform (1) in the JASA paper to retrieve the "observed" latent series and fit an ARMA
    if(mod$nAR || mod$nMA) ARMAFit = itsmr::arma(qnorm(Cxt),mod$nAR,mod$nMA)
    if(mod$nAR) est[(mod$nMargParms+1):(mod$nMargParms+mod$nAR)]                    = ARMAFit$phi
    if(mod$nMA) est[(1+mod$nMargParms+mod$nAR):(mod$nMargParms+sum(mod$ARMAModel))] = ARMAFit$theta
  }

  # add the parmnames on theta fix me: does this affect performance?
  names(est) = mod$parmnames


  return(est)
}

#' Innovations Algorithm for Time Series Prediction
#'
#' Applies the Innovations Algorithm to compute the optimal one-step-ahead predictors
#' and associated prediction error variances based on the autocovariance sequence
#' of a stationary time series. The algorithm uses the ARMA parameters provided in
#' \code{Parms} and the autocovariance vector \code{gamma} to generate recursive
#' estimates of innovation coefficients and variances.
#'
#' @param Parms A list containing the ARMA model parameters:
#'   \itemize{
#'     \item \code{AR}: Vector of autoregressive (AR) coefficients.
#'     \item \code{MA}: Vector of moving average (MA) coefficients.
#'   }
#' @param gamma Numeric vector. The autocovariance sequence of the observed time series.
#'   Typically length \eqn{N}, where \eqn{N} is the time series length.
#'
#' @param mod List. A list containing all model specifications, including:
#'   \describe{
#'     \item{DependentVar}{Numeric. The dependent time series variable to be modeled.}
#'     \item{Task}{Character. The task requested by the user (e.g., Evaluation, Optimization, Synthesis, etc.).}
#'     \item{ParticleNumber}{Integer. Number of particles to use.}
#'     \item{Regressor}{Optional independent variable(s).}
#'     \item{CountDist}{Character string. Specifies the count marginal distribution.}
#'   }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{n}}{Number of recursion steps performed.}
#'   \item{\code{thetas}}{List of innovation coefficients (one per step). Each element is a vector of length \code{q}.}
#'   \item{\code{v}}{Vector of innovation variances at each step (normalized by \code{v[1]}).}
#' }
#'
#' @details
#' The recursion is terminated based on a fixed tolerance (\code{mod$maxdiff}) for variance convergence in MA models,
#' or after a fixed number of steps in AR-only models. The \code{kappa(i,j)} function computes the autocovariance
#' terms needed for recursion using model coefficients and the autocovariance vector \code{gamma}.
#'
#' @examples
#' # ARMA(1,1) example
#' gamma <- ARMAacf(ar = 0.5, ma = 0.4, lag.max = 20)
#' mod <- list(nAR = 1, nMA = 1, maxdiff = 1e-8)
#' Parms <- list(AR = 0.5, MA = 0.4)
#' result <- InnovAlg(Parms, gamma, mod)
#' str(result)
#'
#' @export
InnovAlg <- function(Parms, gamma, mod) {
  # Extract AR and MA coefficients safely
  phi <- if (!is.null(Parms$AR)) Parms$AR else numeric(0)
  theta <- if (!is.null(Parms$MA)) Parms$MA else numeric(0)

  p <- length(phi)
  q <- length(theta)
  m <- max(p, q)

  sigma2 <- 1
  N <- length(gamma)
  theta_r <- c(1, theta, rep(0, N))  # MA coefficients with padding

  # Define the kappa function (autocovariance)
  kappa <- function(i, j) {
    if (j > m) {
      idx1 <- 1:(q + 1)
      idx2 <- (i - j + 1):(i - j + q + 1)
      if (min(idx2) < 1 || max(idx2) > length(theta_r)) return(0)
      return(sum(theta_r[idx1] * theta_r[idx2]))
    } else if (i > 2 * m) {
      return(0)
    } else if (i > m) {
      idx <- abs((1 - i + j):(p - i + j)) + 1
      idx <- idx[idx >= 1 & idx <= length(gamma)]
      return((gamma[i - j + 1] - sum(phi * gamma[idx])) / sigma2)
    } else {
      return(gamma[i - j + 1] / sigma2)
    }
  }

  # Initialize
  Theta <- list()
  v <- rep(NA_real_, N + 1)
  v[1] <- kappa(1, 1)
  StopCondition <- FALSE
  n <- 1

  while (!StopCondition && n < N) {
    Theta[[n]] <- rep(0, q)

    denom <- v[1]
    if (!is.finite(denom) || denom == 0) break

    Theta[[n]][n] <- kappa(n + 1, 1) / denom
    if (n > q && mod$nAR == 0) Theta[[n]][n] <- 0

    if (n > 1) {
      for (k in 1:(n - 1)) {
        js <- 0:(k - 1)
        v_k1 <- v[k + 1]
        if (!is.finite(v_k1) || v_k1 == 0) next

        Theta[[n]][n - k] <- (kappa(n + 1, k + 1) -
                                sum(Theta[[k]][k - js] * Theta[[n]][n - js] * v[js + 1])) / v_k1
      }
    }

    js <- 0:(n - 1)
    v[n + 1] <- kappa(n + 1, n + 1) - sum((Theta[[n]][n - js])^2 * v[js + 1])

    # Stopping rules depending on model structure
    if (mod$nAR == 0 && mod$nMA > 0) {
      StopCondition <- abs(v[n + 1] - v[n]) < mod$maxdiff
    } else if (mod$nAR > 0 && mod$nMA == 0) {
      StopCondition <- n > 3 * m
    } else {
      StopCondition <- n >= N || abs(v[n + 1] - v[n]) < mod$maxdiff
    }

    n <- n + 1
  }

  v <- v / v[1]
  return(list(
    n = n - 1,
    thetas = lapply(Theta[1:(n - 1)], function(x) x[1:q]),
    v = v[1:(n - 1)]
  ))
}


#' Simulate Count Time Series from Specified Model
#'
#' Simulates a univariate count time series of length \code{n} from a specified model
#' defined by a marginal distribution, ARMA dependence structure, and optional regression
#' covariates.
#'
#' @param n Integer. Length of the time series to simulate.
#' @param CountDist Character. Name of the marginal count distribution to simulate from.
#'   Supported values include \code{"Poisson"}, \code{"Negative Binomial"}, \code{"ZIP"},
#'   \code{"Generalized Poisson"}, \code{"Mixed Poisson"}, \code{"Binomial"}, etc.
#' @param MargParm Numeric vector. Parameters of the marginal distribution (e.g., mean, dispersion).
#' @param ARParm Numeric vector. Autoregressive (AR) parameters. Can be empty or \code{NULL} for no AR component.
#' @param MAParm Numeric vector. Moving average (MA) parameters. Can be empty or \code{NULL} for no MA component.
#' @param Regressor Optional matrix or data frame of covariates. Used for dynamic regression models.
#' @param Intercept Logical or numeric. Whether to include an intercept term. If numeric, should be 0 or 1.
#' @param ntrials Integer. Number of trials for the Binomial model (required if \code{CountDist == "Binomial"}).
#' @param ... Additional arguments (currently unused).
#'
#' @return A ts object containing the simulated count time series.
#'
#' @details
#' The simulation is based on a latent Gaussian copula framework, where the marginal distribution
#' is applied to transformed latent variables with dependence induced by an ARMA process.
#'
#' @examples
#' # Simulate Poisson-ARMA(1,1) series
#' sim <- sim_lgc(
#'   n = 100,
#'   CountDist = "Poisson",
#'   MargParm = 3,
#'   ARParm = 0.5,
#'   MAParm = -0.3
#' )
#' ts.plot(sim)
#'
#' @export
sim_lgc = function(n, CountDist, MargParm, ARParm, MAParm, Regressor=NULL, Intercept=NULL, ntrials=NULL,...){

  # combine all parameters in on vector
  TrueParam = c(MargParm,ARParm,MAParm)

  # set ARMAModel orders
  ARMAModel = c(length(ARParm), length(MAParm))

  # add a column of ones in the Regressors if Intercept is present
  if (!is.null(Intercept) && Intercept==TRUE && sum(as.matrix(Regressor)[,1])!=n){
    Regressor = as.data.frame(cbind(rep(1,n),Regressor))
    names(Regressor)[1] = "Intercept"
  }

  # parse the model information - needed for distribution functions
  mod = ModelScheme(SampleSize = n,
                     CountDist = CountDist,
                     ARMAModel = ARMAModel,
                     TrueParam = TrueParam,
                     Regressor = Regressor,
                     Intercept = Intercept,
                          Task = "Synthesis",
                       ntrials = ntrials)

  # reorganize the parameters in expected format
  Parms = RetrieveParameters(TrueParam,mod)

  # Generate latent Gaussian model
  z  =arima.sim(model = list( ar = ARParm, ma=MAParm), n = n)
  z = z/sd(z) # standardize the data

  # apply the inverse count cdf - check me: in Binomial case ntrials comes as global
  if(mod$nreg==0){
    x = mod$myinvcdf(pnorm(z), Parms$MargParms)
  }else{
    x = mod$myinvcdf(pnorm(z), Parms$ConstMargParm, Parms$DynamMargParm)
  }
  # FIX ME: Should I be returning ts objects?
  # return(as.numeric(x))
  return(ts(x))
}


#' Compute AIC, BIC, and AICc for Model Evaluation
#'
#' Computes the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC),
#' and corrected AIC (AICc) for a fitted model, based on the log-likelihood, sample size,
#' and number of parameters.
#'
#' @param loglik Numeric. The log-likelihood value of the fitted model.
#' @param mod A list containing model-related metadata, typically created by
#'   \code{\link{ModelScheme}}. It must include:
#'   \itemize{
#'     \item \code{nparms}: Number of estimated parameters.
#'     \item \code{n}: Sample size.
#'     \item \code{EstMethod}: Estimation method used ("GL" or "PFR").
#'   }
#'
#' @return A numeric vector of length 3 containing the values:
#' \describe{
#'   \item{\code{AIC}}{Akaike Information Criterion}
#'   \item{\code{BIC}}{Bayesian Information Criterion}
#'   \item{\code{AICc}}{Corrected Akaike Information Criterion}
#' }
#'
#' @examples
#' mod <- list(nparms = 5, n = 100, EstMethod = "PFR")
#' loglik <- -123.45
#' Criteria.lgc(loglik, mod)
#'
#' @export
Criteria.lgc = function(loglik, mod){

  # initial versions of the function included Gaussial Likelihood as estimation method
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

#' @title Log-Likelihood for lgc Model
#' @description Returns the log-likelihood value stored in a fitted `lgc` model object.
#' @param object An object of class `lgc`.
#' @param ... Additional arguments (currently unused).
#' @return Numeric. Log-likelihood value.
#' @exportS3Method logLik lgc
logLik.lgc = function(object,...){
  return(object$FitStatistics[1])
}

#' @title AIC for lgc Model
#' @description Returns the AIC value for a fitted `lgc` model.
#' @param object An object of class `lgc`.
#' @param ... Additional arguments (currently unused).
#' @return Numeric. AIC value.
#' @exportS3Method AIC lgc
AIC.lgc = function(object,...){
  return(object$FitStatistics[2])
}

BIC <- function(object, ...) UseMethod("BIC")

#' @title BIC for lgc Model
#' @description Returns the BIC value for a fitted `lgc` model.
#' @param object An object of class `lgc`.
#' @param ... Additional arguments (currently unused).
#' @return Numeric. BIC value.
#' @exportS3Method BIC lgc
BIC.lgc = function(object,...){
  return(object$FitStatistics[3])
}

se <- function(object, ...) UseMethod("se")

#' @title Standard Errors for lgc Model
#' @description Returns the standard errors from a fitted `lgc` model.
#' @param object An object of class `lgc`.
#' @param ... Additional arguments (currently unused).
#' @return Matrix of standard errors.
#' @exportS3Method se lgc
se.lgc = function(object,...){
  return(object$StdErrors)
}

coefficients <- function(object, ...) UseMethod("coefficients")

#' @title Coefficients for lgc Model
#' @description Returns the estimated coefficients from a fitted `lgc` model.
#' @param object An object of class `lgc`.
#' @param ... Additional arguments (currently unused).
#' @return Matrix of parameter estimates.
#' @exportS3Method coefficients lgc
coefficients.lgc = function(object,...){
  return(object$ParamEstimates)
}

model <- function(object, ...) UseMethod("model")

#' @title Model Summary for lgc Object
#' @description Returns a data frame summarizing the distribution and ARMA model structure.
#' @param object An object of class `lgc`.
#' @param ... Additional arguments (currently unused).
#' @return A data frame with two columns: distribution and model type.
#' @exportS3Method model lgc
model.lgc = function(object,...){
  # if ((object$ARMAModel[1]>0) &&  (object$ARMAModel[2]>0)){
  #   ARMAModel = sprintf("ARMA(%.0f, %.0f)",object$ARMAModel[1], object$ARMAModel[2])
  # }
  # if ((object$ARMAModel[1]>0) &&  (object$ARMAModel[2]==0)){
  #   ARMAModel = sprintf("AR(%.0f)",object$ARMAModel[1])
  # }
  # if ((object$ARMAModel[1]==0) &&  (object$ARMAModel[2]>0)){
  #   ARMAModel = sprintf("MA(%.0f)",object$ARMAModel[2])
  # }
  # if ((object$ARMAModel[1]==0) &&  (object$ARMAModel[2]==0)){
  #   ARMAModel = "White Noise"
  # }
  #
  # a = data.frame(object$CountDist, ARMAModel)
  # names(a) = c("Distribution", "Model")
  # return(a)
  return(object$Model)
}


#' @title Extract Residuals from an lgc Model
#' @description Computes residuals from a fitted \code{lgc} model using the
#'   ARMA structure specified in the model.
#' @param object An object of class \code{lgc}.
#' @param ... Additional arguments (currently unused).
#' @return A numeric vector of residuals of length \code{nobs(object)}.
#' @examples
#' \dontrun{
#'   residuals(mylgc)
#' }
#' @exportS3Method residuals lgc
residuals.lgc <- function(object, ...) {
  if (!inherits(object, "lgc")) stop("Object must be of class 'lgc'")
  res <- ComputeResiduals(object)
  res
}



#' @title Summarize an lgc Model
#' @description Provides a summary of a fitted \code{lgc} model, including
#'   parameter estimates and fit statistics.
#' @param object An object of class \code{lgc}.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns the original \code{lgc} object.
#' @examples
#' \dontrun{
#'   summary(mylgc)
#' }
#' @exportS3Method summary lgc
summary.lgc <- function(object, ...) {
  cat("Model:\n")
  print(object$Model)

  cat("\nParameter Estimates:\n")
  print(object$ParamEstimates)

  if (!is.null(object$SE)) {
    cat("\nStandard Errors:\n")
    print(object$SE)
  }

  cat("\nFit Statistics:\n")
  print(object$FitStatistics)

  invisible(object)
}

# print <- function(object, ...) UseMethod("print")

#' @title Print a Latent Gaussian Count (LGC) Model
#' @description
#' Nicely formats and prints a fitted \code{lgc} model object, including model type,
#' estimation details, fit statistics, and parameter estimates. Designed to provide
#' a concise summary suitable for quick inspection by statisticians and
#' time-series practitioners.
#'
#' @param x An object of class \code{"lgc"}.
#' @param digits Integer; number of significant digits to display (default = 4).
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This method provides a human-readable printout of an LGC model object.
#' For more comprehensive model information, including diagnostic metrics
#' and statistical inference, use \code{summary.lgc()}.
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @examples
#' \dontrun{
#'   # Fit an example LGC model
#'   mylgc <- lgc(formula = y ~ x, data = simdata, EstMethod = "PFR")
#'
#'   # Print the model summary
#'   print(mylgc)
#' }
#'
#' @seealso \code{\link{summary.lgc}} for a detailed summary method.
#' @export
#' @method print lgc
print.lgc <- function(x, digits = 4, ...) {
  # cat("========================================\n")
  # cat("     Latent Gaussian Count Model (LGC)  \n")
  # cat("========================================\n")
  cat("Model:", x$Model, "\n")
  if (!is.null(x$CountDist)) cat("Count Dist.:    ", x$CountDist, "\n")
  if (!is.null(x$ARMAModel)) cat("ARMA Structure: ", x$ARMAModel, "\n")
  #if (!is.null(x$SampleSize)) cat("Sample Size:    ", x$SampleSize, "\n")
  #if (!is.null(x$EstMethod)) cat("Estimation:     ", x$EstMethod, "\n")

  # Add optimization method if available
  # if (!is.null(x$mod) && !is.null(x$mod$OptMethod)) {
  #   cat("Optimizer:      ", x$mod$OptMethod, "\n")
  # } else if (!is.null(x$OptMethod)) {
  #   # fallback if optimizer stored directly in object
  #   cat("Optimizer:      ", x$OptMethod, "\n")
  # }
  #
  # cat("\n")
  cat("\n")
  ## --- Parameter estimates (+ optional SEs) ---
  if (!is.null(x$ParamEstimates)) {
    cat("Parameter Estimates:\n")

    # pull estimates
    est <- as.numeric(x$ParamEstimates)
    names(est) <- colnames(x$ParamEstimates)

    # pull std errors if present
    if (!is.null(x$StdErrors)) {
      se <- as.numeric(x$StdErrors)
      names(se) <- names(est)
    } else {
      se <- rep(NA_real_, length(est))
    }

    # how wide should each column be?
    # we take the max of: parameter name length vs formatted number length
    digits <- 4
    est_fmt <- formatC(est, digits = digits, format = "f")
    se_fmt  <- ifelse(is.na(se), "", formatC(se, digits = digits, format = "f"))

    col_widths <- pmax(
      nchar(names(est)),
      nchar(est_fmt),
      nchar(se_fmt),
      6L # minimum width so it's not too cramped
    ) + 2L  # add some padding between columns

    # helper: print a row given a character vector (one cell per param)
    print_row <- function(cells, left_label = "", left_width = 6) {
      cat(format(left_label, width = left_width, justify = "right"))
      for (j in seq_along(cells)) {
        cat(format(cells[j], width = col_widths[j], justify = "right"))
      }
      cat("\n")
    }

    # 1) header row: parameter names
    print_row(names(est), left_label = "", left_width = 6)

    # 2) estimates row
    print_row(est_fmt, left_label = "", left_width = 6)

    # 3) std. errors row (only if we have any)
    if (any(!is.na(se))) {
      print_row(se_fmt, left_label = "s.e.", left_width = 6)
    }

    cat("\n")
  } else {
    cat("  (No parameter estimates available)\n\n")
  }



  ## --- Fit statistics ---
  if (!is.null(x$FitStatistics)) {
    cat("Fit Statistics:\n")
    fs <- as.data.frame(x$FitStatistics)
    print(round(fs, digits), row.names = FALSE)
    cat("\n")
  }

  ## --- Optimization diagnostics ---
  # if (!is.null(x$OptimOutput)) {
  #   cat("Optimization Diagnostics:\n")
  #   opt <- as.data.frame(x$OptimOutput)
  #   print(opt, row.names = FALSE)
  # }

  ## --- Other info ---
  # if (!is.null(x$Converged))  cat("Converged:      ", ifelse(x$Converged, "Yes", "No"), "\n")
  # if (!is.null(x$NumParams))  cat("Num Parameters: ", x$NumParams, "\n")



  invisible(x)
}

#' @title Number of Observations in lgc Model
#' @description Returns the number of observations used in fitting the model.
#' @param object An object of class \code{lgc}.
#' @param ... Additional arguments (currently unused).
#' @return An integer, the number of observations.
#' @exportS3Method stats::nobs lgc
nobs.lgc <- function(object, ...) {
  return(as.integer(object$SampleSize))
}


#' Compute Truncation Limits for Latent Gaussian Count Models
#'
#' Computes the lower and upper truncation limits used in the sampling of the latent variable
#' in the particle filter algorithm. Corresponds to equation (19) in the JASA paper.
#'
#' @param mod A list containing model specification and data, typically created by \code{\link{ModelScheme}}.
#' @param Parms A list of model parameters retrieved via \code{\link{RetrieveParameters}}.
#' @param t Integer. Current time point.
#' @param Zhat Numeric vector of predicted latent values up to time \code{t}.
#' @param Rt Numeric vector of innovation standard deviations from the innovations algorithm.
#'
#' @return A list with elements \code{a} and \code{b} representing the lower and upper limits
#' for truncated normal sampling.
#' @keywords internal
#' @seealso \code{\link{SampleTruncNormParticles}}
#' @export
ComputeLimits = function(mod, Parms, t, Zhat, Rt){

  # a and b are the arguments in the two normal cdfs in the 4th line in equation (19) in JASA paper
  Lim = list()
  # fix me: this will be ok for AR or MA models but for ARMA? is it p+q instead of max(p,q)
  # fix me: why do I have t(Parms$MargParms)?
  index = min(t, max(mod$ARMAModel))

  # add the following for White Noise models
  if(max(mod$ARMAModel)==0){index=1}

  if(mod$nreg==0){
    Lim$a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,t(Parms$MargParms)),0,1)) - Zhat)/(Rt[index])
    Lim$b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],t(Parms$MargParms)),0,1)) - Zhat)/Rt[index]
  }else{
    Lim$a = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t]-1,Parms$ConstMargParm, Parms$DynamMargParm[t,]),0,1)) - Zhat)/Rt[index]
    Lim$b = as.numeric((qnorm(mod$mycdf(mod$DependentVar[t],Parms$ConstMargParm, Parms$DynamMargParm[t,]),0,1)) - Zhat)/Rt[index]
  }

  return(Lim)
}

#' Sample Truncated Normal Particles
#'
#' Generates samples from a truncated normal distribution using inverse transform sampling,
#' as required by the particle filter for latent Gaussian models.
#'
#' @param mod Model list as returned by \code{\link{ModelScheme}}.
#' @param Limit A list of truncation limits with elements \code{a} and \code{b}, typically
#' returned by \code{\link{ComputeLimits}}.
#' @param t Integer. Current time index.
#' @param Zhat Numeric vector of predicted latent values up to time \code{t}.
#' @param Rt Numeric vector of innovation standard deviations from the innovations algorithm.
#'
#' @return A numeric vector of sampled latent variables.
#' @keywords internal
#' @export
SampleTruncNormParticles = function(mod, Limit, t, Zhat, Rt){
  # relation (21) in JASA paper and the inverse transform method
  # check me: this can be improved?
  index = min(t, max(mod$ARMAModel))
  if(max(mod$ARMAModel)==0){index=1}
  z = qnorm(runif(length(Limit$a),0,1)*(pnorm(Limit$b,0,1)-pnorm(Limit$a,0,1))+pnorm(Limit$a,0,1),0,1)*Rt[index] + Zhat
  return(z)
}

#' Compute Filtered Latent Value at Time t
#'
#' Computes the conditional mean of the latent variable at time \code{t} using the innovations
#' algorithm output.
#'
#' @param mod Model list as returned by \code{\link{ModelScheme}}.
#' @param IA A list returned by \code{\link{InnovAlg}}, containing innovation variances and
#' theta coefficients.
#' @param Z Matrix of latent particle values.
#' @param Zhat Matrix of predicted latent means.
#' @param t Integer. Time point at which to compute the latent mean.
#' @param Parms List of model parameters from \code{\link{RetrieveParameters}}.
#'
#' @return A numeric vector with the filtered latent value at time \code{t}.
#' @keywords internal
#' @export
ComputeZhat_t = function(mod, IA, Z, Zhat,t, Parms){

  Theta    = IA$thetas
  nTheta   = length(IA$thetas)
  Theta_n  = Theta[[nTheta]]

  m = max(mod$ARMAModel)
  p = mod$ARMAModel[1]
  q = mod$ARMAModel[2]

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

#' Compute Particle Weights
#'
#' Computes the importance weights of the particles at time \code{t}, based on the
#' truncation limits.
#'
#' @param mod Model list as returned by \code{\link{ModelScheme}}.
#' @param Limit A list of truncation limits (from \code{\link{ComputeLimits}}).
#' @param t Integer. Current time index.
#' @param PreviousWeights Numeric vector of weights at time \code{t-1}.
#'
#' @return A numeric vector of updated particle weights at time \code{t}.
#' @keywords internal
#' @export
ComputeWeights = function(mod, Limit, t, PreviousWeights){
  # equation (21) in JASA paper
  # update weights
  if(t<=max(mod$ARMAModel)){
    NewWeights = PreviousWeights*(pnorm(Limit$b,0,1) - pnorm(Limit$a,0,1))
  }else{ # fix me: if I add the wgh[t-1,] below as I should the weights become small?
    NewWeights = (pnorm(Limit$b,0,1) - pnorm(Limit$a,0,1))
  }

  return(NewWeights)
}

#' Resample Particles Using Effective Sample Size (ESS) Criterion
#'
#' Performs resampling of particles when the effective sample size falls below a
#' threshold defined by the model. Follows relation (26) in the JASA paper.
#'
#' @param mod Model list as returned by \code{\link{ModelScheme}}.
#' @param wgh Matrix of particle weights for all time points.
#' @param t Integer. Current time index.
#' @param Znew Vector of new particle samples to be resampled.
#'
#' @return A resampled vector of particles of length equal to \code{mod$ParticleNumber}.
#' @keywords internal
#' @export
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

#' Extract model Parameters from Flat Vector
#'
#' Decomposes a flat parameter vector \code{theta} into structured components
#' according to the model specification provided in \code{mod}. This includes
#' marginal distribution parameters, autoregressive (AR) and moving average (MA)
#' components, and dynamic GLM-type parameters when covariates are present.
#'
#' @param theta Numeric vector. The complete parameter vector, typically as passed
#'   to or from an optimizer.
#' @param mod A list containing the model specification, as produced by
#'   \code{\link{ModelScheme}}. It must contain fields such as:
#'   \itemize{
#'     \item \code{MargParmIndices}: Indices of marginal parameters in \code{theta}
#'     \item \code{CountDist}: Name of the count distribution
#'     \item \code{Regressor}: Covariate matrix (if applicable)
#'     \item \code{nreg}, \code{nint}, \code{nMargParms}: Model structure details
#'     \item \code{ARMAModel}: Integer vector specifying AR and MA orders
#'   }
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{\code{MargParms}}{Marginal distribution parameters (from \code{theta}).}
#'   \item{\code{ConstMargParm}}{Constant component(s) of the marginal distribution (if applicable).}
#'   \item{\code{DynamMargParm}}{Dynamic component(s), such as covariate effects (if applicable).}
#'   \item{\code{AR}}{Autoregressive coefficients (if present).}
#'   \item{\code{MA}}{Moving average coefficients (if present).}
#' }
#'
#' @keywords internal
#' @export
RetrieveParameters = function(theta,mod){

  Parms =  vector(mode = "list", length = 5)


  names(Parms) = c("MargParms", "ConstMargParm", "DynamMargParm", "AR", "MA")

  # marginal parameters
  Parms$MargParms      = theta[mod$MargParmIndices]

  # regressor parameters
  if(mod$nreg>0){
    beta  = Parms$MargParms[1:(mod$nreg+mod$nint)]
    m     = exp(as.matrix(mod$Regressor)%*%beta)
  }

  # GLM type parameters
  if(mod$CountDist == "Negative Binomial" && mod$nreg>0){
    Parms$ConstMargParm  = 1/Parms$MargParms[mod$nreg+mod$nint+1]
    Parms$DynamMargParm  = Parms$MargParms[mod$nreg+mod$nint+1]*m/(1+Parms$MargParms[mod$nreg+mod$nint+1]*m)
  }

  if(mod$CountDist == "Generalized Poisson" && mod$nreg>0){
    Parms$ConstMargParm  = Parms$MargParms[mod$nreg+mod$nint+1]
    Parms$DynamMargParm  = m
  }

  if(mod$CountDist == "Generalized Poisson 2" && mod$nreg>0){
    Parms$ConstMargParm  = Parms$MargParms[mod$nreg+mod$nint+1]
    Parms$DynamMargParm  = m
  }

  if(mod$CountDist == "Poisson" && mod$nreg>0){
    Parms$ConstMargParm  = NULL
    Parms$DynamMargParm  = m
  }

  if(mod$CountDist == "Binomial" && mod$nreg>0){
    Parms$ConstMargParm  = NULL
    Parms$DynamMargParm  = m/(1+m)
  }

  if(mod$CountDist == "Mixed Poisson" && mod$nreg>0){
    Parms$ConstMargParm  = c(Parms$MargParms[2*(mod$nreg+mod$nint)+1])
    Parms$DynamMargParm  = cbind(m, exp(as.matrix(mod$Regressor)%*%Parms$MargParms[(mod$nreg+mod$nint+1):((mod$nreg+mod$nint)*2)]))
  }

  if(mod$CountDist == "ZIP" && mod$nreg>0){
    Parms$ConstMargParm  = Parms$MargParms[mod$nreg+mod$nint+1]
    Parms$DynamMargParm  = m
  }


  # Parms$DynamMargParm = as.matrix(Parms$DynamMargParm)

  # Parms$AR = NULL
  if(mod$ARMAModel[1]>0) Parms$AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAModel[1])]

  # Parms$MA = NULL
  if(mod$ARMAModel[2]>0) Parms$MA = theta[(mod$nMargParms+mod$ARMAModel[1]+1) :
                                            (mod$nMargParms + mod$ARMAModel[1] + mod$ARMAModel[2]) ]


  return(Parms)
}


#------------------------------ MIXED POISSON----------------------------------------------#
#' Vectorized Quantile Function for Mixed Poisson Distribution
#'
#' Computes quantiles of the mixed Poisson distribution using a vectorized approach
#' with \code{mapply()}.
#'
#' @param p Numeric vector of probabilities (between 0 and 1).
#' @param lam1 Numeric or vector. Mean(s) of the first Poisson component.
#' @param lam2 Numeric or vector. Mean(s) of the second Poisson component.
#' @param prob Numeric between 0 and 1. Mixing probability for the first component.
#'
#' @return A numeric vector of quantiles corresponding to \code{p}.
#'
#' @details
#' This function is a vectorized version of \code{\link{qmixpois1}} using \code{mapply()}
#' to compute quantiles efficiently.
#'
#' @examples
#' qmixpois1(c(0.1, 0.5, 0.9), lam1 = 2, lam2 = 5, prob = 0.6)
#'
#' @export
qmixpois1 = function(p, lam1, lam2, prob) {
  # Ensure p, lam1, and lam2 are vectors
  p = as.vector(p)

  # Replicate lam1 and lam2 if they are constants
  if (length(lam1) == 1) {
    lam1 = rep(lam1, length(p))
  }
  if (length(lam2) == 1) {
    lam2 = rep(lam2, length(p))
  }

  # Function to calculate quantile for a single set of parameters
  find_quantile = function(pi, lambda1, lambda2, prob) {
    x = 0
    while (pmixpois1(x, lambda1, lambda2, prob) <= pi) {
      x = x + 1
    }
    return(x)
  }

  # Use mapply to apply the function over vectors p, lam1, and lam2
  quantiles = mapply(find_quantile, p, lam1, lam2, MoreArgs = list(prob = prob))

  return(quantiles)
}

#' Mixture of Two Poisson Distributions (CDF)
#'
#' Computes the cumulative distribution function of a two-component mixture of Poisson distributions.
#'
#' @param x Non-negative integer(s). Vector of quantiles.
#' @param lam1 Mean of the first Poisson component (must be > 0).
#' @param lam2 Mean of the second Poisson component (must be > 0).
#' @param p Mixing probability for the first component (0 <= p <= 1).
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)}.
#'   Otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities p are given as \code{log(p)}.
#'
#' @return Numeric vector of mixture CDF values (or log-probabilities).
#' @examples
#' # CDF at 3 under the mixture
#' pmixpois1(3, lam1 = 2, lam2 = 5, p = 0.4)
#'
#' # Upper tail probability in log scale
#' pmixpois1(3, lam1 = 2, lam2 = 5, p = 0.4, lower.tail = FALSE, log.p = TRUE)
#'
#' @export
pmixpois1 <- function(x, lam1, lam2, p, lower.tail = TRUE, log.p = FALSE) {
  cdf <- p * ppois(x, lam1) + (1 - p) * ppois(x, lam2)
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  cdf
}


#' Mixture of Two Poisson Distributions (PMF)
#'
#' Computes the probability mass function of a two-component mixture of Poisson distributions.
#'
#' @param x Non-negative integer(s). Vector of quantiles.
#' @param lam1 Mean of the first Poisson component (must be > 0).
#' @param lam2 Mean of the second Poisson component (must be > 0).
#' @param p Mixing probability for the first component (0 <= p <= 1).
#' @param log Logical; if \code{TRUE}, probabilities p are given as \code{log(p)}.
#'
#' @return Numeric vector of mixture probabilities (or log-probabilities).
#' @examples
#' # Probability of observing exactly 3 under the mixture
#' dmixpois1(3, lam1 = 2, lam2 = 5, p = 0.4)
#'
#' # Log-probability
#' dmixpois1(3, lam1 = 2, lam2 = 5, p = 0.4, log = TRUE)
#'
#' # Vectorized input
#' dmixpois1(0:5, lam1 = 2, lam2 = 5, p = 0.4)
#'
#' @export
dmixpois1 <- function(x, lam1, lam2, p, log = FALSE) {
  dens <- p * dpois(x, lam1) + (1 - p) * dpois(x, lam2)
  if (log) dens <- log(dens)
  dens
}


#-----------------------------------------------------------------------------------------------#


#------------------------------------    Generalized Poisson        -----------------------------#
# dGpois, pGpois, qGpois and rGpois are used when CountDist = "Generalize Poisson". I have now
# implemented the "Generalized Poisson 2" and in the future the following may not be needed.

#' Generalized Poisson PMF
#'
#' Computes the probability mass function (PMF) of the Generalized Poisson distribution
#' using the mean-parametrization (Famoye 1994, relation 2.4).
#'
#' @param y Non-negative integer vector of quantiles.
#' @param a Dispersion parameter.
#' @param m Mean parameter.
#' @param log Logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return A numeric vector of (log-)probabilities.
#'
#' @export
dGpois <- function(y, a, m, log = FALSE) {
  k <- m / (1 + a * m)
  dens <- k^y * (1 + a * y)^(y - 1) * exp(-k * (1 + a * y) - lgamma(y + 1))
  if (log) dens <- log(dens)
  dens
}


#' Generalized Poisson CDF
#'
#' Computes the cumulative distribution function (CDF) of the Generalized Poisson distribution.
#'
#' @param x Integer or numeric vector. Values at which to evaluate the CDF.
#' @param a Numeric. Dispersion parameter.
#' @param m Numeric scalar or vector. Mean parameter(s).
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)}.
#'   Otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @return A numeric vector or matrix of cumulative probabilities.
#'
#' @seealso \code{\link{dGpois}}
#'
#' @examples
#' pGpois(3, 0.4, 2)
#' pGpois(3, 0.4, 2, lower.tail = FALSE, log.p = TRUE)
#'
#' @export
pGpois <- function(x, a, m, lower.tail = TRUE, log.p = FALSE) {

  compute_cdf <- function(xi, ai, mi) {
    M <- max(xi, 0)
    probs <- dGpois(0:M, ai, mi)
    cum_probs <- cumsum(probs)

    out <- rep(0, length(xi))
    out[xi >= 0] <- cum_probs[xi[xi >= 0] + 1]

    # handle lower.tail
    if (!lower.tail) {
      out <- 1 - out
    }
    # handle log.p
    if (log.p) {
      out <- log(out)
    }
    out
  }

  if (length(m) > 1 || length(a) > 1) {
    mapply(compute_cdf, x = x, ai = a, mi = m, SIMPLIFY = TRUE)
  } else {
    compute_cdf(x, a, m)
  }
}


#' Generalized Poisson Quantile Function (vectorized and scalar)
#'
#' Computes the quantile function (inverse CDF) of the Generalized Poisson distribution
#' for one or more values of the mean parameter \code{m}.
#'
#' @param p Numeric vector of probabilities (between 0 and 1).
#' @param a Numeric. Dispersion parameter.
#' @param m Numeric scalar or vector. Mean parameter(s).
#'
#' @return A numeric vector or matrix of quantiles.
#'
#' @seealso \code{\link{pGpois}}
#'
#' @examples
#' qGpois(0.8,0.4,3)
#'
#' @export
qGpois <- function(p, a, m) {

  compute_quantile <- function(pi, ai, mi) {
    cdf_vals <- pGpois(0:100, ai, mi)
    quantiles <- rep(NA_real_, length(pi))

    for (i in seq_along(pi)) {
      ix <- which(cdf_vals >= pi[i])[1]
      quantiles[i] <- if (!is.na(ix)) ix - 1 else NA_real_
    }

    return(quantiles)
  }

  if (length(m) > 1 || length(a) > 1) {
    mapply(compute_quantile, pi = p, ai = a, mi = m, SIMPLIFY = TRUE)
  } else {
    compute_quantile(p, a, m)
  }
}
#-----------------------------------------------------------------------------------------------#


#' Generate Random Model Parameters for Testing
#'
#' Randomly generates marginal distribution and ARMA parameters for simulation and
#' testing purposes. Includes a mechanism to occasionally generate "bad" parameters
#' (e.g., invalid or edge-case values) with controlled probability, useful for
#' stress-testing estimation and validation routines. This function was uised in initial stages
#' of testing and may need to be reworked for more sophisticated testing in the future
#'
#' @param CountDist Character string. The name of the count distribution. Supported
#'   values include \code{"Poisson"}, \code{"Negative Binomial"}, \code{"Mixed Poisson"}, and \code{"ZIP"}.
#' @param BadParamProb Numeric between 0 and 1. Probability of sampling "good" vs. "bad"
#'   parameter values. A value close to 1 favors well-behaved parameters.
#' @param AROrder Integer. Order of the autoregressive (AR) component. If 0, no AR terms are generated.
#' @param MAOrder Integer. Order of the moving average (MA) component. If 0, no MA terms are generated.
#' @param Regressor Optional matrix or data frame of covariates. If supplied, dynamic
#'   models with covariate-based parameterization are used.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{\code{MargParm}}{Vector of marginal distribution parameters. Structure depends on \code{CountDist} and whether \code{Regressor} is provided.}
#'   \item{\code{ARParm}}{Autoregressive parameter(s), or \code{NULL} if \code{AROrder = 0}.}
#'   \item{\code{MAParm}}{Moving average parameter(s), or \code{NULL} if \code{MAOrder = 0}.}
#' }
#'
#' @details
#' For testing robustness of model fitting procedures, this function sometimes samples
#' marginal or ARMA parameters from values known to be problematic (e.g., negative rates,
#' non-stationary AR roots). This is controlled by the \code{BadParamProb} parameter.
#'
#' When \code{Regressor} is not \code{NULL}, the function returns parameters assuming a log-linear
#' model with regression effects (e.g., intercept and slope coefficients).
#'
#' @keywords internal
#' @export
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


  if(CountDist=="Mixed Poisson"){
    # create some "bad" Poisson parameter choices
    BadLambda = c(-1,500)

    # sample a probability
    p =  runif(1,0,1)

    # sample with high probability from "good" choices for lambda and with low prob the "bad" choices
    prob = rbinom(2,1,BadParamProb)

    # Marginal Parameter
    if (is.null(Regressor)){
      lambda1 = prob[1]*runif(1,0,100) + (1-prob[1])*sample(BadLambda,1)
      lambda2 = prob[2]*runif(1,0,100) + (1-prob[2])*sample(BadLambda,1)
      MargParm = c(lambda1, lambda2, p)
    }else{
      # fix me: I am hard coding 1.2 here but we should probably make this an argument
      b0 = runif(1,0,1.2)
      b1 = runif(1,0,1.2)
      c0 = runif(1,0,1.2)
      c1 = runif(1,0,1.2)
      k  = runif(1,0,1)
      MargParm = c(b0,b1,c0,c1,k)
    }
  }


  if(CountDist=="ZIP"){
    # create some "bad" ZIP parameter choices
    BadLambda = c(-1,500)

    # sample a probability
    p =  runif(1,0,1)

    # sample with high probability from "good" choices for lambda and with low prob the "bad" choices
    prob = rbinom(1,1,BadParamProb)

    # Marginal Parameter
    if (is.null(Regressor)){
      lambda = prob*runif(1,0,100) + (1-prob)*sample(BadLambda,1)
      MargParm = c(lambda, p)
    }else{
      # fix me: I am hard coding 1.2 here but we should probably make this an argument
      b0 = runif(1,0,1.2)
      b1 = runif(1,0,1.2)
      c0 = runif(1,0,1.2)
      c1 = runif(1,0,1.2)
      MargParm = c(b0,b1,c0, c1)
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

#' Generate Perturbed Initial Parameter Values for Testing
#'
#' Applies a small perturbation to a given set of model parameters to generate
#' initial values for optimization or testing purposes.
#'
#' @param AllParms A named list containing model parameters with components:
#'   \itemize{
#'     \item \code{MargParm}: vector of marginal distribution parameters
#'     \item \code{ARParm}: vector of autoregressive (AR) parameters, or \code{NULL}
#'     \item \code{MAParm}: vector of moving average (MA) parameters, or \code{NULL}
#'   }
#' @param perturbation Numeric scalar. Proportion by which to perturb the parameters
#'   (e.g., \code{0.1} applies a 10% change).
#'
#' @return A named list with the same structure as \code{AllParms}, containing
#' perturbed parameter values. If any of the original parameters are \code{NULL},
#' they remain \code{NULL} in the output.
#'
#' @details
#' The perturbation is applied multiplicatively and negatively:
#' \deqn{\theta_{\text{init}} = \theta \cdot (1 - \text{perturbation})}
#' This helps avoid starting optimization from true parameters, which could
#' mask numerical or convergence issues.
#'
#' @examples
#' true_parms <- list(MargParm = c(1, 0.5), ARParm = c(0.3), MAParm = NULL)
#' GenInitVal(true_parms, perturbation = 0.1)
#'
#' @keywords internal
#' @export
GenInitVal = function(AllParms, perturbation){

  # select negative or positive perturbation based on a  coin flip
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


#' Computes residuals from a fitted latent Gaussian count (LGC) model using the approach
#' described in relation (41) of Jia et al. (2021).
#'
#' @param lgc A list object containing the outcome of a model fit, typically returned by
#' the package's main wrapper function \code{\link{lgc}}.
#'
#' @return A numeric vector of residuals, computed via the inverse-normal transformation of the
#' cumulative distribution function (CDF) of the observed counts.
#'
#' @details
#' This function computes residuals as defined in relation (41) of Jia et al. (2021)
#'
#' @references
#' Jia, Y., Kechagias, S., Livsey, J., Lund, R., & Pipiras, V. (2021).
#' Latent Gaussian Count Time Series.
#' \emph{Journal of the American Statistical Association}, 118(541), 596–606.
#' \doi{10.1080/01621459.2021.1944874}
#'
#' @export
ComputeResiduals= function(lgc){

  # retrieve marginal distribution parameters
  theta      = lgc$ParamEstimates
  mod        = lgc$mod
  Parms      = RetrieveParameters(theta,mod)
  nMargParms = length(Parms$MargParms)
  nparms     = length(theta)

  # check if the number of parameters matches the model setting
  if(nMargParms + sum(mod$ARMAModel)!=nparms) stop('The length of the provided estimated parameters does not
                                                   match the model specification.')

  # allocate memory
  Zhat <- rep(NA,mod$n)

  # pGpois function doesn't allow for vector lambda so a for-loop is needed, however the change shouldn't
  # be too hard. note the formula subtracts 1 from the counts mod$DependentVar so it may lead to negative numbers
  # that's why there is and if else below.
  for (i in 1:mod$n){
    k <- mod$DependentVar[i]
    if (k != 0) {
      # Compute limits
      if(mod$nreg==0){
        C_xt_minus1 = mod$mycdf(k-1,t(Parms$MargParms))
        C_xt        = mod$mycdf(k,t(Parms$MargParms))
        }else{
          C_xt_minus1 = mod$mycdf(k-1, Parms$ConstMargParm, Parms$DynamMargParm[i])
          C_xt        = mod$mycdf(k, Parms$ConstMargParm, Parms$DynamMargParm[i])
      }
      a           = qnorm(C_xt_minus1,0,1)
      b           = qnorm(C_xt ,0,1)
      ResidNum    = exp(-a^2/2)-exp(-b^2/2)
      ResidDenom  = sqrt(2*pi)*(C_xt-C_xt_minus1)
    }else{
      if(mod$nreg==0){
        # Compute integral limits
        C_xt        = mod$mycdf(0,t(Parms$MargParms))
        b           = qnorm(C_xt,0,1)
      }else{
        C_xt        = mod$mycdf(0, Parms$ConstMargParm, Parms$DynamMargParm[i])
        b           = qnorm(C_xt ,0,1)
      }
      ResidNum    = -exp(-b^2/2)
      ResidDenom  = sqrt(2*pi)*C_xt
    }
    Zhat[i]     =  ResidNum/ResidDenom
  }

  # apply AR filter--Fix me allow for MA as well -commenting out the code below to use arima
  # function. I am also adding centering step.
  # residualOld = data.frame(filter(Zhat,c(1,-Parms$AR))[1:(mod$n-mod$ARMAModel[1])])
  # names(residualOld) = "residual"

  # demean the data - initially the centering step was skipped
  DemeanRes =  Zhat - mean(Zhat)

  # Fit ARMA with fixed parameters
  fit = arima(DemeanRes,
               order = c(mod$nAR,0,mod$nMA),
               fixed = c(Parms$AR,Parms$MA, 0), # AR1, AR2, AR3, intercept
               transform.pars = FALSE)

  residual = as.numeric(residuals(fit))
  names(residual) <- NULL

  return(residual)
}


#' Constructs a  string that describes the marginal count distribution.
#' This function is used to help with printing warnings.
#'
#' @param theta Numeric vector. The full parameter vector, from which the marginal
#'   parameters will be extracted using model structure from \code{mod}.
#' @param mod A list object containing model specifications. Must include:
#'   \itemize{
#'     \item \code{CountDist}: character name of the marginal distribution (e.g., "Poisson")
#'     \item \code{parmnames}: character vector of parameter names
#'     \item \code{nMargParms}: number of marginal parameters
#'   }
#'
#' @return A character string of the form \code{"DistName(param1=..., param2=..., ...)"}.
#'
#' @examples
#' mod <- list(
#'   CountDist = "Poisson",
#'   parmnames = c("lambda"),
#'   nMargParms = 1
#' )
#' theta <- c(2.35)
#' CurrentDist(theta, mod)
#' # Returns: "Poisson(lambda=2.35)"
#'
#' @keywords internal
#' @export
CurrentDist = function(theta,mod){
  # check me: does it work for models with regressors?
  name = sprintf("%s=%.2f",mod$parmnames[1],theta[1])

  if (mod$nMargParms>1){
    for(i in 2:mod$nMargParms){
      #name = paste(name,initialParam[i],sep=",")
      a = sprintf("%s=%.2f",mod$parmnames[i],theta[i])
      name = paste(name,a,sep=",")
    }
  }
  name = paste(mod$CountDist, "(", name,")", sep="")
  return(name)
}


#' Modifies the standard particle filtering procedure to return the predictive distribution
#' at each time point of a latent Gaussian count time series model, instead of just the log-likelihood.
#'
#' @param theta Numeric vector. Parameter vector containing marginal and ARMA parameters.
#' @inheritParams ParticleFilter_Res_ARMA
#'
#' @return A numeric matrix of dimension \code{2 x n}, where each column corresponds to a time point:
#'   \itemize{
#'     \item First row: lower tail probability for the observed count
#'     \item Second row: predictive probability mass at \code{Y_t + 1}
#'   }
#'
#' @details
#' This function runs a particle filter with resampling to estimate the predictive distribution
#' of the latent state at each time point. It leverages the Innovations Algorithm to calculate
#' ARMA-based forecasts of the latent process and approximates conditional distributions using
#' truncated normal draws. See Section 3.4 in Jia et al. (2021)
#'
#' @references
#' Jia, Y., Kechagias, S., Livsey, J., Lund, R., & Pipiras, V. (2021).
#' Latent Gaussian Count Time Series.
#' \emph{Journal of the American Statistical Association}, 118(541), 596–606.
#' \doi{10.1080/01621459.2021.1944874}
#'
#' @seealso \code{\link{ComputeWeights}}, \code{\link{SampleTruncNormParticles}}, \code{\link{InnovAlg}}
#'
#' @examples
# # specify model
#' n              = 10
#' Regressor      = data.frame(runif(n),runif(n))
#' Intercept      = TRUE
#' ARMAModel      = c(2,0)
#' ARParm         = c(0.5, 0.2)
#' MAParm         = NULL
#' CountDist      = "Poisson"
#' b0             = 1
#' b1             = 4
#' b2             = 2
#' MargParm       = c(b0,b1,b2)
#'
#' # simulate data
#' set.seed(2)
#' DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, Regressor, Intercept)
#'
#' mod = ModelScheme(DependentVar   = DependentVar,
#'                   Regressor      = Regressor,
#'                   Intercept      = Intercept,
#'                   CountDist      = CountDist,
#'                  ARMAModel      = ARMAModel)
#'
#' # select a parameter point
#' theta <- c(MargParm, ARParm, MAParm)
#'
#' # compute the Predictive distribution
#' PDvalues(theta, mod)
#'
#' @export
PDvalues = function(theta, mod){
  #------------------------------------------------------------------------------------#
  # PURPOSE:  Compute predictive distribution
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
  if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))

  # Initialize the negative log likelihood computation
  nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
                   - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1])))

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

  # to collect the values of predictive distribution of interest
  preddist = matrix(0,2,mod$n)

  # initialize particle filter weights
  w[1,]    = rep(1,mod$ParticleNumber)

  # Compute the first integral limits Limit$a and Limit$b
  Limit    = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # Initialize the particles using N(0,1) variables truncated to the limits computed above
  #Z[1,]    = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
  Z[1,]    = SampleTruncNormParticles(mod, Limit, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))


  for (t in 2:mod$n){

    # compute Zhat_t
    Zhat[t,] = ComputeZhat_t(mod, IA, Z, Zhat,t, Parms)

    # compute a,b and temp
    temp  = rep(0,(mod$DependentVar[t]+1))
    index = min(t, max(mod$ARMAModel))
    for (x in 0:mod$DependentVar[t]){
      # Compute integral limits
      if(mod$nreg==0){
        a = as.numeric(qnorm(mod$mycdf(x-1,Parms$MargParms),0,1) -  Zhat[t,])/Rt[index]
        b = as.numeric(qnorm(mod$mycdf(x,  Parms$MargParms),0,1) -  Zhat[t,])/Rt[index]
      }else{
        a = as.numeric(qnorm(mod$mycdf(x-1,Parms$ConstMargParm, Parms$DynamMargParm[t]),0,1) -  Zhat[t,])/Rt[index]
        b = as.numeric(qnorm(mod$mycdf(x,  Parms$ConstMargParm, Parms$DynamMargParm[t]),0,1) -  Zhat[t,])/Rt[index]
      }
      temp[x+1] = mean(pnorm(b,0,1) - pnorm(a,0,1))
    }

    # Compute integral limits
    Limit = ComputeLimits(mod, Parms, t, Zhat[t,], Rt)

    # Sample truncated normal particles
    Znew  = SampleTruncNormParticles(mod, Limit, t, Zhat[t,], Rt)

    # update weights
    w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      message(sprintf('WARNING: At t=%.0f some of the weights are either too small or sum to 0',t))
      return(10^8)
    }
    # Resample the particles using common random numbers
    old_state1 = get_rand_state()
    Znew = ResampleParticles(mod, w, t, Znew)
    set_rand_state(old_state1)

    # save the current particle
    Z[t,]   = Znew


    if (mod$DependentVar[t]==0){
      preddist[,t] = c(0,temp[1])
    }else{
      preddist[,t] = cumsum(temp)[mod$DependentVar[t]:(mod$DependentVar[t]+1)]
    }

  }

  return(preddist)
}

#' Compute Probability Integral Transform (PIT) Histogram Values
#'
#' Computes histogram values for the Probability Integral Transform (PIT) as described
#' in relations (39)–(40) of the JASA paper on latent Gaussian count time series modeling.
#' This function processes predictive distributions and outputs the PIT histogram bins
#' for calibration assessment.
#'
#' @param H Integer. Number of histogram bins to divide the \[0,1\] PIT range.
#' @param predDist A numeric matrix of dimension \code{2 x n}, typically returned by
#'   \code{\link{PDvalues}}. Each column corresponds to a time point, with:
#'   \itemize{
#'     \item First row: lower tail probability of the observed count
#'     \item Second row: cumulative probability including the observed count
#'   }
#'
#' @return A numeric vector of length \code{H} representing the PIT histogram values (bin heights).
#' These can be plotted to assess calibration of predictive distributions.
#'
#' @details
#' The PIT is used to assess the calibration of probabilistic forecasts. For discrete data,
#' the PIT is computed using randomized methods or as a piecewise function, as detailed in
#' the latent Gaussian count models literature.
#'
#' This implementation applies the formulas:
#' \deqn{PIT(Y_t) \sim U(0,1)} if predictive distributions are correctly specified.
#'
#' @references
#' Jia, Y., Kechagias, S., Livsey, J., Lund, R., & Pipiras, V. (2021).
#' Latent Gaussian Count Time Series.
#' \emph{Journal of the American Statistical Association}, 118(541), 596–606.
#' \doi{10.1080/01621459.2021.1944874}
#'
#' @seealso \code{\link{PDvalues}}, \code{\link{ComputeResiduals}}
#'
#' @examples
#' \dontrun{
#' predDist <- PDvalues(theta, mod)
#' hist(PITvalues(10, predDist), breaks = 10, main = "PIT Histogram")
#' }
#'
#' @export
PITvalues = function(H, predDist){
  PITvalues = rep(0,H)

  predd1 = predDist[1,]
  predd2 = predDist[2,]
  Tr = length(predd1)

  for (h in 1:H){
    id1 = (predd1 < h/H)*(h/H < predd2)
    id2 = (h/H >= predd2)
    tmp1 = (h/H-predd1)/(predd2-predd1)
    tmp1[!id1] = 0
    tmp2 = rep(0,Tr)
    tmp2[id2] = 1
    PITvalues[h] = mean(tmp1+tmp2)
  }
PITvalues = c(0,PITvalues)
return(diff(PITvalues))
}


#' Parse Formula into Response and Regressors
#'
#' Parses a standard R formula and extracts the response variable, predictor variables,
#' and intercept status.
#'
#' @param formula An object of class \code{formula}, typically in the form \code{y ~ x1 + x2}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{DependentVar}}{Character string of the response variable name.}
#'   \item{\code{Regressor}}{Character vector of predictor variable names, or \code{NULL} if none.}
#'   \item{\code{intercept}}{Logical, \code{TRUE} if an intercept is included, \code{FALSE} otherwise.}
#' }
#'
#' @details
#' This function is useful for preprocessing model formulas in custom modeling functions.
#' It uses \code{\link{terms}} to extract structure from the formula object.
#'
#'
#' @seealso \code{\link{terms}}, \code{\link{model.frame}}
#'
#' @export
parse_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("Input must be a formula.")
  }

  # Extract the terms object from the formula
  terms_obj <- terms(formula)

  # Extract the response variable (left-hand side)
  DependentVar <- as.character(formula[[2]])

  # Extract the predictor terms (right-hand side)
  Regressor <- attr(terms_obj, "term.labels")

  # If there are no predictors, return NULL
  if (length(Regressor) == 0) {
    Regressor <- NULL
  }

  # Check if the intercept is included (1 if included, 0 if excluded)
  has_intercept <- attr(terms_obj, "intercept") == 1

  # Return a list with response, predictors, and intercept status
  return(list(
    DependentVar = DependentVar,
    Regressor = Regressor,
    intercept = has_intercept
  ))
}

#' Dimension function for vector/matrix inputs
#' Treats vectors as length x 1 column vectors
#'
#' @param x vector (column matrix) or matrix
#'
#' @return the vector's dimension
#' @keywords internal
#' @export
DIM <- function(x){
  if (is.null(dim(x)))
    c(length(x), 1)
  else
    dim(x)
}

#===========================================================================================#
#--------------------------------------------------------------------------------#
# On August 2024 for Issue #27, I created a modified likelihood function called
# ParticleFilter with three changes:
#
# 1. The value 1 for the count CDF calculations is changed to 1-10^(-16), so that
#    when inverse normal is applied in the copula, we do not receive an inf value.
# 2. When the limits of the truncated normal distribution become too large (say>7),
#    I set them equal to 7-epsilon, so when I apply the normal cdf and subsequnetly the
#    inverse cdf i avoid the inf values.
# 3. If the weights in the particle filter likelihood become 0 I set them equal to
#    10^(-64).
#
# I have added the changes in the ParticleFilter, ComputeLimits_MisSpec
# and SampleTruncNormParticles_MisSpec files and will continue  working with the
# standard versions below until I perform adequate testing.


#' @title Particle Filter Latent Gaussian Count Models
#'
#' @description
#' Implements a particle filtering algorithm for latent Gaussian count time series models
#' with potentially misspecified ARMA dependence structures.
#' The algorithm sequentially updates latent Gaussian states, particle weights, and the log-likelihood,
#' using the Innovations Algorithm to obtain one-step-ahead predictors and innovation variances.
#' Corrections are applied to truncated normal sampling limits and particle weights
#' when numerical instability or boundary issues occur.
#'
#' @param theta Numeric vector. Current parameter estimates including marginal, AR, and MA parameters.
#'
#' @param mod List. Model specification object containing all necessary elements for filtering, such as:
#'   \describe{
#'     \item{\code{DependentVar}}{Numeric vector of observed counts.}
#'     \item{\code{CountDist}}{Character string specifying the count distribution (e.g., \code{"Poisson"}, \code{"Binomial"}).}
#'     \item{\code{ARMAModel}}{Numeric vector \code{c(p, q)} specifying AR and MA orders.}
#'     \item{\code{ParticleNumber}}{Integer. Number of particles used in the filter.}
#'     \item{\code{n}}{Integer. Length of the observed time series.}
#'     \item{\code{mycdf}}{User-supplied function returning the marginal CDF of the count distribution.}
#'     \item{\code{mypdf}}{User-supplied function returning the marginal PDF of the count distribution.}
#'     \item{\code{Regressor}}{Optional numeric vector or matrix of regressors.}
#'     \item{\code{nreg}}{Integer. Number of regressors.}
#'     \item{\code{verbose}}{Logical. If \code{TRUE}, prints diagnostic messages during filtering.}
#'     \item{\code{maxdiff}}{Numeric tolerance used in the Innovations Algorithm convergence criterion.}
#'   }
#'
#' @details
#' The function evaluates the (negative) log-likelihood for an LGC model under ARMA dynamics.
#'
#' @return
#' Returns the numeric value of the negative log-likelihood (\eqn{-\ell(\theta)}).
#' If numerical issues occur (e.g., unstable ARMA coefficients, degenerate weights, or non-finite log-likelihood),
#' the function returns a large penalty value (\eqn{10^8}) to guide optimization away from invalid parameter regions.
#
#' @export
ParticleFilter = function(theta, mod){

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # Retrieve parameters and save them in a list called Parms
  Parms = RetrieveParameters(theta,mod)

  # check for causality
  if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))

  # Initialize the negative log likelihood computation
  if(mod$nreg==0){
    nloglik = - log(max(mod$mypdf(mod$DependentVar[1],Parms$MargParms),.Machine$double.xmin))
  }else{
    nloglik = - log(max(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1]),.Machine$double.xmin))
  }

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

  # initialize a list to keep track of manual corrections
  Corrections = list(
            weights = NULL,
           weights2 = NULL,
             C_xt_1 = NULL,
               C_xt = NULL,
                a_t = NULL,
                b_t = NULL,
    SampleParticles = NULL
  )

  # Compute the first integral limits Limit$a and Limit$b
  Limit = ComputeLimits_MisSpec(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # if correction for C_xt_1 was necessary for time point t=1, then save it
  if (Limit$Correction_C_xt_1){
    Corrections$C_xt_1 = 1
  }

  # if correction for C_xt was necessary for time point t=1, then save it
  if (Limit$Correction_C_xt){
    Corrections$C_xt = 1
  }

  # Initialize the particles using N(0,1) variables truncated to the limits computed above
  #Z[1,]    = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
  SampleParticles = SampleTruncNormParticles_MisSpec(mod, Limit, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
  Z[1,]    = SampleParticles$z

  # if correction for lower limit was necessary for time point t=1, then save it
  if (SampleParticles$Correction_a){
    Corrections$a_t = 1
  }

  # if correction for upper limit was necessary for time point t=1, then save it
  if (SampleParticles$Correction_b){
    Corrections$b_t = 1
  }

  for (t in 2:mod$n){

    # compute Zhat_t
    #Zhat[t,] = ComputeZhat_t(m,Theta,Z,Zhat,t, Parms,p,q, nTheta, Theta_n)
    Zhat[t,] = ComputeZhat_t(mod, IA, Z, Zhat,t, Parms)


    # Compute integral limits
    Limit = ComputeLimits_MisSpec(mod, Parms, t, Zhat[t,], Rt)

    # if correction for C_xt_1 was necessary for time point t=1, then save it
    if (Limit$Correction_C_xt_1){
      Corrections$C_xt_1 = c(Corrections$C_xt_1,t)
    }

    # if correction for C_xt was necessary for time point t=1, then save it
    if (Limit$Correction_C_xt){
      Corrections$C_xt = c(Corrections$C_xt,t)
    }

    # Sample truncated normal particles
    #Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b, t, Zhat[t,], Rt)
    SampleParticles = SampleTruncNormParticles_MisSpec(mod, Limit, Parms, t, Zhat[t,], Rt)
    Znew  = SampleParticles$z

    # if correction for lower limit was necessary for time point t, then save it
    if (SampleParticles$Correction_a){
      Corrections$a_t = c(Corrections$a_t,t)
    }

    # if correction for upper limit was necessary for time point t, then save it
    if (SampleParticles$Correction_b){
      Corrections$b_t = c(Corrections$b_t,t)
    }


    # update weights
    #w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])
    w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])

    # check me: In misspecified models, the weights may get equal to 0. Is it ok
    # for me to do the following? how is this different from allowing zero weights and
    # returning a large likelihood?

    if (sum(w[t,])==0){
      w[t,] = rep(10^(-64),mod$ParticleNumber)
      Corrections$weights = c(Corrections$weights,t)
    }

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      Corrections$weights2 = c(Corrections$weights2,t)
      if(mod$verbose){
        message(sprintf('WARNING: At t=%.0f some of the weights are either too small or sum to 0.\n',t))
      }
      return(10^8)
    }

    # Resample the particles using common random numbers
    old_state1 = get_rand_state()
    Znew = ResampleParticles(mod, w, t, Znew)
    set_rand_state(old_state1)


    # save the current particle
    Z[t,]   = Znew

    #print(t)
    #print(Z[t,])
    # print(nloglik)
    # update likelihood
    nloglik = nloglik - log(mean(w[t,]))
  }

  # report messages
  if(mod$verbose){
    ReportDiagnostics(Corrections, Parms, mod)
  }
  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  # if (nloglik==Inf | is.na(nloglik)){
  #   nloglik = 10^8
  # }

  return(nloglik)
}

# gather the messages that should be reported in a likelihood computation
ReportDiagnostics = function(Corrections, Parms, mod){

  # if there were corrections performed for C_xt_1, print an appropriate message
  if(!is.null(Corrections$weights)){
    message(
      sprintf(
        paste0(
          "At t = %s,\n",
          "weights sum to 0. They are all resetted to 10^(-64).\n\n"
        ),
        paste(  Corrections$weights, collapse = ", ")
      )
    )
  }

  # if there were corrections performed for C_xt_1, print an appropriate message
  if(!is.null(Corrections$C_xt_1)){
    message(
      sprintf(
        paste0(
          "The Count CDF takes the value 1 for some time series values. When plugged in the inverse Normal\n",
          "CDF this produces an Inf value. To avoid this, the Count CDF is set to 1-10^(-16).\n\n",

          "This behavior can be caused from a misspecified marginal distribution,\n",
          "missing covariates, or the optimizer exploring unreasonable parameter values.\n\n",

          "A %s model was fitted.\n",
          "Corrrections occured at time point(s) t = %s,\n",
          "where the dependent variable takes value(s): %s."
        ),
        CurrentDist(Parms$MargParms, mod),
        paste(Corrections$C_xt_1, collapse = ", "),
        paste(mod$DependentVar[Corrections$C_xt_1], collapse = ", ")
      )
    )
  }

  # if there were corrections performed for C_xt, print an appropriate message
  if(is.null(Corrections$C_xt_1) & !is.null(Corrections$C_xt)){
    message(
      sprintf(
        paste0(
          "The Count CDF takes the value 1 for some time series values. When plugged in the inverse Normal\n",
          "CDF this produces an Inf value. To avoid this, the Count CDF is set to 1-10^(-16).\n\n",

          "This behavior can be caused from a misspecified marginal distribution,\n",
          "missing covariates, or the optimizer exploring unreasonable parameter values.\n\n",

          "A %s model was fitted.\n",
          "Corrrections occured at time point(s) t = %s,\n",
          "where the dependent variable takes value(s): %s."
        ),
        CurrentDist(Parms$MargParms, mod),
        paste(Corrections$C_xt, collapse = ", "),
        paste(mod$DependentVar[Corrections$C_xt], collapse = ", ")
      )
    )
  }

  # if there were corrections performed for lower limit, print an appropriate message
  if(is.null(Corrections$C_xt_1) & is.null(Corrections$C_xt) & !is.null(Corrections$a_t)){
    message(
      sprintf(
        paste0(
          "The lower limit of the truncated Normal distribution used to generate particles exceeded\n",
          "the value 7, which could to an Inf value from the inverse CDF\n",
          "(e.g., qpois(pnorm(20),1)). To avoid this, the limit was set to 7 - 1e-11.\n\n",

          "Large limit values may result from a misspecified marginal distribution,\n",
          "missing covariates, or the optimizer exploring unreasonable parameter values.\n\n",

          "A %s model was fitted.\n",
          "Corrrections occured at time point(s) t = %s,\n",
          "the dependent variable takes value(s): %s."
        ),
        CurrentDist(Parms$MargParms, mod),
        paste(Corrections$a_t, collapse = ", "),
        paste(mod$DependentVar[Corrections$a_t], collapse = ", ")
      )
    )
  }

  # if there were corrections performed for upper limit, but no corrections for the lower limit print an
  # appropriate message.
  if(is.null(Corrections$C_xt_1) & is.null(Corrections$C_xt) & is.null(Corrections$a_t) & !is.null(Corrections$b_t)){
    message(
      sprintf(
        paste0(
          "The upper limit of the truncated Normal distribution used to generate particles exceeded\n",
          "the value 7, which could to an Inf value from the inverse CDF.\n",
          "(e.g., qpois(pnorm(20),1)). To avoid this, the limit was set to 7 - 1e-11.\n\n",

          "Large limit values may result from a misspecified marginal distribution,\n",
          "missing covariates, or the optimizer exploring unreasonable parameter values.\n\n",

          "A %s model was fitted.\n",
          "Corrrections occured at time point(s) t = %s,\n",
          "the dependent variable takes value(s): %s."
        ),
        CurrentDist(Parms$MargParms, mod),
        paste(Corrections$b_t, collapse = ", "),
        paste(mod$DependentVar[Corrections$b_t], collapse = ", ")
      )
    )
  }




}

ComputeLimits_MisSpec <- function(mod, Parms, t, Zhat, Rt) {
  Lim <- list()

  # index for AR/MA models
  index <- min(t, max(mod$ARMAModel))
  if (max(mod$ARMAModel) == 0) index <- 1   # white noise

  Lim$Correction_C_xt_1 <- FALSE
  Lim$Correction_C_xt   <- FALSE

  # --- Step 1: compute log CDF and log PMF ---
  y <- mod$DependentVar[t]

  if (mod$nreg == 0) {
    logC_1 <- mod$mycdf(y - 1, t(Parms$MargParms), log.p = TRUE)   # log P(X ≤ y-1)
    logpmf <- mod$mypdf(y,     t(Parms$MargParms), log = TRUE)     # log P(X = y)
    logC1U <- mod$mycdf(y - 1, t(Parms$MargParms), lower.tail = FALSE, log.p = TRUE) # log P(X > y-1)
    logCU  <- mod$mycdf(y,     t(Parms$MargParms), lower.tail = FALSE, log.p = TRUE) # log P(X > y)
  } else {
    logC_1 <- mod$mycdf(y - 1, Parms$ConstMargParm, Parms$DynamMargParm[t, ], log.p = TRUE)
    logpmf <- mod$mypdf(y,     Parms$ConstMargParm, Parms$DynamMargParm[t, ], log = TRUE)
    logC1U <- mod$mycdf(y - 1, Parms$ConstMargParm, Parms$DynamMargParm[t, ],
                        lower.tail = FALSE, log.p = TRUE)
    logCU  <- mod$mycdf(y,     Parms$ConstMargParm, Parms$DynamMargParm[t, ],
                        lower.tail = FALSE, log.p = TRUE)
  }

  # --- Step 2: log-space addition for C = log(P(X ≤ y)) ---
  logC <- logspace_add(logC_1, logpmf)

  # --- Step 3: safe qnorm using log-probs directly ---
  safe_q_from_logs <- function(log_lower, log_upper) {
    # log_lower = log P(X ≤ x)
    # log_upper = log P(X > x)
    if (is.finite(log_lower) && log_lower < 0) {
      # small lower tail (probability close to 0)
      return(qnorm(log_lower, log.p = TRUE))
    } else if (is.finite(log_upper) && log_upper < 0) {
      # small upper tail (probability close to 0 → lower tail ~ 1)
      return(-qnorm(log_upper, log.p = TRUE))
    } else {
      # middle case: convert once safely
      return(qnorm(exp(log_lower)))
    }
  }

  # --- Step 4: compute limits ---
  qa <- safe_q_from_logs(logC_1, logC1U)
  qb <- safe_q_from_logs(logC,   logCU)

  Lim$a <- (qa - Zhat) / Rt[index]
  Lim$b <- (qb - Zhat) / Rt[index]

  # --- Step 5: sanity check ---
  if (any(Lim$a > Lim$b) || any(is.na(Lim$a)) || any(is.na(Lim$b))) {
    warning("Problematic limits at time ", t)
  }

  return(Lim)
}


logspace_add <- function(a, b) {
  # Computes log(exp(a) + exp(b)) safely

  if (a > b) {
    # Factor out exp(a) so that exp(b - a) <= 1
    # This avoids overflow if 'a' is very large
    # Use log1p(x) instead of log(1 + x) for accuracy when x is tiny
    return(a + log1p(exp(b - a)))
  } else {
    # Factor out exp(b) so that exp(a - b) <= 1
    # Again use log1p for numerical stability
    return(b + log1p(exp(a - b)))
  }
}


SampleTruncNormParticles_MisSpec = function(mod, Limit, Parms, t, Zhat, Rt){
  # relation (21) in JASA paper and the inverse transform method
  # check me: this can be improved?
  index = min(t, max(mod$ARMAModel))
  if(max(mod$ARMAModel)==0){index=1}

  # Initialize flags that track whether the limits are set to 7 - 10^(-11)
  Correction_a = FALSE
  Correction_b = FALSE

  # check me: this may need to be surfaced in the wrapper. It is small positive constant I will
  # use to get away from zero or one. I need to think this more rigorously.
  #epsilon = 10^(-11)

  #-----------------------------------------------------------------#
  # If the upper limit is less than a threshold (here minus 7) then so will the lower limit
  # if (all(Limit$b<=-38.46) ) {
  #   Limit$b[Limit$b<=-38.46] = -38.46 + epsilon
  #   Limit$a[Limit$a<=-38.46] = -38.46 + epsilon/2
  #   Correction_a  = TRUE
  #   Correction_b  = TRUE
  # }

  # If only the lower limit is less than the threshold we dont need to do anything
  # if (min(Limit$b)>-37 &  min(Limit$a)<=-37 ) {
  #   diff = min(Limit$b - Limit$a)
  #   Limit$a[Limit$a<=-7] = -7 + min(epsilon,diff/2)
  #   Correction_a  = TRUE
  # }
  #-----------------------------------------------------------------#


  #-----------------------------------------------------------------#
  # Check me: in the case of missspecifed models, I may get really large values for the limits
  # wchich means that the pnorm will equal 1 and then the qnorm will yield Inf. Is this ok?
  # If the lower limit is over the threshold, then so will the upper limit
  # if (all(Limit$a>=8.29)) {
  #   Limit$a[Limit$a>=8.29] = 8.29 - epsilon
  #   Limit$b[Limit$b>=8.29] = 8.29 - epsilon/2
  #   Correction_a  = TRUE
  #   Correction_b  = TRUE
  # }

  # If only the upper limit is over the threshold we dont need to do anything
  # if (max(Limit$a)<7 & max(Limit$b)>=7 ) {
  #   diff = min(Limit$b - Limit$a)
  #   Limit$b[Limit$b>=7] = 7 - min(epsilon,diff/2)
  #   Correction_b  = TRUE
  # }
  #-----------------------------------------------------------------#

  z = sample_truncnorm(Limit)*Rt[index] + Zhat
  # z = qnorm(runif(length(Limit$a),0,1)*(pnorm(Limit$b,0,1)-pnorm(Limit$a,0,1))+pnorm(Limit$a,0,1),0,1)*Rt[index] + Zhat
  #z = rtruncnorm(length(Limit$a), Limit$a, Limit$b, mean = 0, sd = 1)*Rt[index] + Zhat

  # create output
  return = list(
    z = z,
    Correction_a = Correction_a,
    Correction_b = Correction_b
  )

  return(return)
}

sample_truncnorm = function(Limit) {
  n <- length(Limit$a)

  # --- manual inverse CDF ---
  u <- runif(n, 0, 1)
  u_shifted <- u * (pnorm(Limit$b, 0, 1) - pnorm(Limit$a, 0, 1)) + pnorm(Limit$a, 0, 1)
  z_manual <- qnorm(u_shifted, mean = 0, sd = 1)

  # --- find problematic draws (Inf / -Inf / NA) ---
  bad_idx <- which(!is.finite(z_manual))

  if (length(bad_idx) > 0) {
    # Replace bad values with rtruncnorm equivalents
    z_manual[bad_idx] <- rtruncnorm(
      length(bad_idx),
      a = Limit$a[bad_idx],
      b = Limit$b[bad_idx],
      mean = 0, sd = 1
    )
  }

  # apply scaling and shift
  z <- z_manual
  return(z)
}

# keeping an older version of this function in case issues show up - I update it in Sep 19.
InnovAlgOld = function(Parms,gamma, mod) {
  # adapt the Innovation.algorithm from ITSMR
  # check me: explain in the doc the exact modifications
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
    Theta[[n]] <- rep(0,q)
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
    if(mod$nAR==0 && mod$nMA>0) StopCondition = (abs(v[n+1]-v[n])< mod$maxdiff)
    if(mod$nAR>0 && mod$nMA==0) StopCondition = (n>3*m)
    n      = n+1
  }
  v = v/v[1]

  StopCondition = (abs(v[n+1]-v[n])< mod$maxdiff)


  return(list(n=n-1,thetas=lapply(Theta[ ][1:(n-1)], function(x) {x[1:q]}),v=v[1:(n-1)]))
}

# writing a log lik function with the old Innov Alg so I can test
ParticleFilter_Res_ARMAOld = function(theta, mod){

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # Retrieve parameters and save them in a list called Parms
  Parms = RetrieveParameters(theta,mod)

  # check for causality
  if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))

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
  IA       = InnovAlgOld(Parms, gamma, mod)
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
  #Z[1,]    = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
  Z[1,]    = SampleTruncNormParticles(mod, Limit, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))


  for (t in 2:mod$n){

    # compute Zhat_t
    #Zhat[t,] = ComputeZhat_t(m,Theta,Z,Zhat,t, Parms,p,q, nTheta, Theta_n)
    Zhat[t,] = ComputeZhat_t(mod, IA, Z, Zhat,t, Parms)

    # Compute integral limits
    Limit = ComputeLimits(mod, Parms, t, Zhat[t,], Rt)

    # Sample truncated normal particles
    #Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b, t, Zhat[t,], Rt)
    Znew  = SampleTruncNormParticles(mod, Limit, t, Zhat[t,], Rt)

    # update weights
    #w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])
    w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      #print(t)
      #print(w[t,])
      message(sprintf('WARNING: At t=%.0f some of the weights are either too small or sum to 0',t))
      return(10^8)
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

# # simulate from our model - the old function
# sim_lgc_old = function(n, CountDist, MargParm, ARParm, MAParm, Regressor=NULL, Intercept=NULL){
#
#   # Generate latent Gaussian model
#   z  =arima.sim(model = list( ar = ARParm, ma=MAParm  ), n = n)
#   z = z/sd(z) # standardize the data
#
#   # add a column of ones in the Regressors if Intercept is present
#   if (!is.null(Intercept) && Intercept==TRUE && sum(as.matrix(Regressor)[,1])!=n){
#     Regressor = as.data.frame(cbind(rep(1,n),Regressor))
#     names(Regressor)[1] = "Intercept"
#   }
#
#
#   # number of regressors
#   nreg = ifelse(is.null(Regressor), 0, dim(Regressor)[2]-as.numeric(Intercept))
#
#   if(nreg==0){
#     # retrieve marginal inverse cdf
#     myinvcdf = switch(CountDist,
#                       "Poisson"               = qpois,
#                       "Negative Binomial"     = function(x, theta){ qnbinom (x, theta[1], 1-theta[2]) },
#                       "Generalized Poisson"   = function(x, theta){ qGpois  (x, theta[1], theta[2])},
#                       "Generalized Poisson 2" = function(x, theta){ qGpois  (x, theta[2], theta[1])},
#                       "Binomial"              = qbinom,
#                       "Mixed Poisson"         = function(x, theta){qmixpois1(x, theta[1], theta[2], theta[3])},
#                       "ZIP"                   = function(x, theta){ qzipois(x, theta[1], theta[2]) },
#                       stop("The specified distribution is not supported.")
#     )
#
#     # get the final counts
#     x = myinvcdf(pnorm(z), MargParm)
#
#   }else{
#     # retrieve inverse count cdf
#     myinvcdf = switch(CountDist,
#                       "Poisson"               = function(x, ConstMargParm, DynamMargParm){ qpois   (x, DynamMargParm)},
#                       "Negative Binomial"     = function(x, ConstMargParm, DynamMargParm){ qnbinom (x, ConstMargParm, 1-DynamMargParm)},
#                       "Generalized Poisson"   = function(x, ConstMargParm, DynamMargParm){ qGpois  (x, ConstMargParm, DynamMargParm)},
#                       "Generalized Poisson 2" = function(x, ConstMargParm, DynamMargParm){ qgenpois2  (x, DynamMargParm, ConstMargParm)},
#                       "Binomial"              = function(x, ConstMargParm, DynamMargParm){ qbinom  (x, ConstMargParm, DynamMargParm)},
#                       "Mixed Poisson"         = function(x, ConstMargParm, DynamMargParm){qmixpois1(x, DynamMargParm[,1], DynamMargParm[,2],  ConstMargParm)},
#                       "ZIP"                   = function(x, ConstMargParm, DynamMargParm){ qzipois (x, DynamMargParm[,1], DynamMargParm[,2]) },
#                       stop("The specified distribution is not supported.")
#     )
#
#     # regression parameters
#     beta  = MargParm[1:(nreg+1)]
#     m     = exp(as.matrix(Regressor)%*%beta)
#
#     if(CountDist == "Poisson" && nreg>0){
#       ConstMargParm  = NULL
#       DynamMargParm  = m
#     }
#
#     if(CountDist == "Negative Binomial" && nreg>0){
#       ConstMargParm  = 1/MargParm[nreg+2]
#       DynamMargParm  = MargParm[nreg+2]*m/(1+MargParm[nreg+2]*m)
#     }
#
#     if(CountDist == "Generalized Poisson" && nreg>0){
#       ConstMargParm  = MargParm[nreg+2]
#       DynamMargParm  = m
#     }
#
#
#     if(CountDist == "Generalized Poisson 2" && nreg>0){
#       ConstMargParm  = MargParm[nreg+2]
#       DynamMargParm  = m
#     }
#
#     if(CountDist == "Binomial" && nreg>0){
#       # ConstMargParm  = MargParm[nreg+2]
#       DynamMargParm  = m/(1+m)
#     }
#
#     if(CountDist == "Mixed Poisson" && nreg>0){
#       ConstMargParm  = c(MargParm[nreg*2+3], 1 - MargParm[nreg*2+3])
#       DynamMargParm  = cbind(exp(as.matrix(Regressor)%*%MargParm[1:(nreg+1)]),
#                              exp(as.matrix(Regressor)%*%MargParm[(nreg+2):(nreg*2+2)]))
#     }
#
#     if(CountDist == "ZIP" && nreg>0){
#       ConstMargParm  = NULL
#       DynamMargParm  = cbind(exp(as.matrix(Regressor)%*%MargParm[1:(nreg+1)]),
#                              1/(1+exp(-as.matrix(Regressor)%*%MargParm[(nreg+2):(nreg*2+2)])))
#     }
#
#     # get the final counts
#     x = myinvcdf(pnorm(z), ConstMargParm, DynamMargParm)
#   }
#   return(as.numeric(x))
# }


# PF likelihood with resampling for ARMA(p,q) - with adhoc truncations, handling case when C==1 not in log space
ParticleFilter_Res_ARMA_MisSpecOld = function(theta, mod){

  old_state <- get_rand_state()
  on.exit(set_rand_state(old_state))

  # Retrieve parameters and save them in a list called Parms
  Parms = RetrieveParameters(theta,mod)

  # check for causality
  if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))

  # Initialize the negative log likelihood computation
  if(mod$nreg==0){
    nloglik = - log(max(mod$mypdf(mod$DependentVar[1],Parms$MargParms),.Machine$double.xmin))
  }else{
    nloglik = - log(max(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1]),.Machine$double.xmin))
  }

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

  # initialize a list to keep track of manual corrections
  Corrections = list(
    weights = NULL,
    weights2 = NULL,
    C_xt_1 = NULL,
    C_xt = NULL,
    a_t = NULL,
    b_t = NULL,
    SampleParticles = NULL
  )

  # Compute the first integral limits Limit$a and Limit$b
  Limit = ComputeLimits_MisSpecOld(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))

  # if correction for C_xt_1 was necessary for time point t=1, then save it
  if (Limit$Correction_C_xt_1){
    Corrections$C_xt_1 = 1
  }

  # if correction for C_xt was necessary for time point t=1, then save it
  if (Limit$Correction_C_xt){
    Corrections$C_xt = 1
  }

  # Initialize the particles using N(0,1) variables truncated to the limits computed above
  #Z[1,]    = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
  SampleParticles = SampleTruncNormParticles_MisSpec(mod, Limit, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
  Z[1,]    = SampleParticles$z

  # if correction for lower limit was necessary for time point t=1, then save it
  if (SampleParticles$Correction_a){
    Corrections$a_t = 1
  }

  # if correction for upper limit was necessary for time point t=1, then save it
  if (SampleParticles$Correction_b){
    Corrections$b_t = 1
  }

  for (t in 2:mod$n){

    # compute Zhat_t
    #Zhat[t,] = ComputeZhat_t(m,Theta,Z,Zhat,t, Parms,p,q, nTheta, Theta_n)
    Zhat[t,] = ComputeZhat_t(mod, IA, Z, Zhat,t, Parms)


    # Compute integral limits
    Limit = ComputeLimits_MisSpecOld(mod, Parms, t, Zhat[t,], Rt)

    # if correction for C_xt_1 was necessary for time point t=1, then save it
    if (Limit$Correction_C_xt_1){
      Corrections$C_xt_1 = c(Corrections$C_xt_1,t)
    }

    # if correction for C_xt was necessary for time point t=1, then save it
    if (Limit$Correction_C_xt){
      Corrections$C_xt = c(Corrections$C_xt,t)
    }

    # Sample truncated normal particles
    #Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b, t, Zhat[t,], Rt)
    SampleParticles = SampleTruncNormParticles_MisSpec(mod, Limit, Parms, t, Zhat[t,], Rt)
    Znew  = SampleParticles$z

    # if correction for lower limit was necessary for time point t, then save it
    if (SampleParticles$Correction_a){
      Corrections$a_t = c(Corrections$a_t,t)
    }

    # if correction for upper limit was necessary for time point t, then save it
    if (SampleParticles$Correction_b){
      Corrections$b_t = c(Corrections$b_t,t)
    }

    # update weights
    #w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])
    w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])

    # check me: In misspecified models, the weights may get equal to 0. Is it ok
    # for me to do the following? how is this different from allowing zero weights and
    # returning a large likelihood?

    if (sum(w[t,])==0){
      w[t,] = rep(10^(-64),mod$ParticleNumber)
      Corrections$weights = c(Corrections$weights,t)
    }

    # check me: break if I got NA weight
    if (any(is.na(w[t,]))| sum(w[t,])==0 ){
      Corrections$weights2 = c(Corrections$weights2,t)
      if(mod$verbose){
        message(sprintf('WARNING: At t=%.0f some of the weights are either too small or sum to 0.\n',t))
      }
      return(10^8)
    }

    # Resample the particles using common random numbers
    old_state1 = get_rand_state()
    Znew = ResampleParticles(mod, w, t, Znew)
    set_rand_state(old_state1)

    # save the current particle
    Z[t,]   = Znew

    #print(t)
    #print(Z[t,])
    # print(nloglik)
    # update likelihood
    nloglik = nloglik - log(mean(w[t,]))
  }

  # report messages
  if(mod$verbose){
    ReportDiagnostics(Corrections, Parms, mod)
  }
  # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
  # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))

  # if (nloglik==Inf | is.na(nloglik)){
  #   nloglik = 10^8
  # }

  return(nloglik)
}

# handling case when C==1 not in log space
ComputeLimits_MisSpecOld = function(mod, Parms, t, Zhat, Rt){
  # a and b are the arguments in the two normal cdfs in the 4th line in equation (19) in JASA paper
  Lim = list()
  # fix me: this will be ok for AR or MA models but for ARMA? is it p+q instead of max(p,q)
  # fix me: why do I have t(Parms$MargParms)?
  index = min(t, max(mod$ARMAModel))

  # add the following for White Noise models
  if(max(mod$ARMAModel)==0) index=1

  # check me: this may need to be surfaced in the wrapper. It is small positive constant I will
  # use to get away from zero or one. I need to think this more rigorously.
  epsilon = 10^(-16)

  # Initialize flags that track whether the limits are set to 7 - 10^(-11)
  Lim$Correction_C_xt_1 = FALSE
  Lim$Correction_C_xt = FALSE

  # Compute C_xt and C_xt_1 in relation (19), Jia et al. (2021).
  if(mod$nreg==0){
    C_1 = mod$mycdf(mod$DependentVar[t]-1,t(Parms$MargParms))
    C   = mod$mycdf(mod$DependentVar[t],t(Parms$MargParms))
  }else{
    C_1 = mod$mycdf(mod$DependentVar[t]-1,Parms$ConstMargParm, Parms$DynamMargParm[t,])
    C   = mod$mycdf(mod$DependentVar[t],Parms$ConstMargParm, Parms$DynamMargParm[t,])
  }

  if (C_1==1){
    C_1_upper = mod$mycdf(mod$DependentVar[t]-1,t(Parms$MargParms),lower.tail = FALSE)
  }
  if (C==1){
    C_upper = mod$mycdf(mod$DependentVar[t],t(Parms$MargParms),lower.tail = FALSE)
  }
  # if only C_1 is zero and C>0 we add an epsilon to C_1, but make sure that C_1 < C.
  # if(C_1==0 & C>0) {
  #   C_1 = 0 + min(epsilon, 0.9*C)
  #   Lim$Correction_C_xt_1 = TRUE
  # }

  # if C==0 then C_1 should also be equal to zero. Then both need to be corrected but I still need C_1 < C.
  # if(C==0 ) {
  #   C_1 = 0 + epsilon/2
  #   C   = 0 + epsilon
  #   Lim$Correction_C_xt_1 = TRUE
  #   Lim$Correction_C_xt = TRUE
  # }

  # If C=1 and C_1 is smaller than 1 then subtract an epsilon from C, but we still need C_1 < C.
  # if(C==1 & C_1<1) {
  #   diff = C-C_1
  #   C = 1 - min(epsilon,diff/2)
  #   Lim$Correction_C_xt = TRUE
  # }

  # If C_1 = 1 then C should also be equal to 1. Then both need to be corrected but I still need C_1 < C.
  # if(C_1==1) {
  #   C_1 = 1 - epsilon
  #   C = 1 - epsilon/2
  #   Lim$Correction_C_xt_1 = TRUE
  #   Lim$Correction_C_xt = TRUE
  # }

  # compute the limits
  if(C_1<1){
    Lim$a = as.numeric((qnorm(C_1,0,1)) - Zhat)/Rt[index]
  }else{
    Lim$a = as.numeric((-qnorm(C_1_upper ,0,1)) - Zhat)/Rt[index]
  }

  if(C<1){
    Lim$b = as.numeric((qnorm(C  ,0,1)) - Zhat)/Rt[index]
  }else{
    Lim$b = as.numeric((-qnorm(C_upper  ,0,1)) - Zhat)/Rt[index]
  }

  #Lim$a = as.numeric((safe_qnorm(C_1)) - Zhat)/Rt[index]
  #Lim$b = as.numeric((safe_qnorm(C)) - Zhat)/Rt[index]
  return(Lim)
}

#' Computes the quantiles (inverse CDF) of the two-component mixed Poisson distribution
#' using a brute-force loop-based method.
# qmixpoisOld = function(y, lam1, lam2, p){
#   yl = length(y)  # Length of the input vector y
#   x  = rep(0, yl) # Initialize a vector of zeros to store results
#
#   for (n in 1:yl) {
#     lambda1 = ifelse(length(lam1) > 1, lam1[n], lam1)  # Use indexed lambda1 or constant
#     lambda2 = ifelse(length(lam2) > 1, lam2[n], lam2)  # Use indexed lambda2 or constant
#
#     # Increment x[n] until the CDF is greater than y[n]
#     while (pmixpois1(x[n], lambda1, lambda2, p) <= y[n]) {
#       x[n] = x[n] + 1
#     }
#   }
#
#   return(x)
# }



#############################################################################################
#-------------------------------------------------------------------------------------------#
# functions we wrote, and may use in the future
# poisson cdf using incomplete gamma and its derivative wrt to lambda

# myppois = function(x, lambda){
#   # compute poisson cdf as the ratio of an incomplete gamma function over the standard gamma function
#   # I will also compute the derivative of the poisson cdf wrt lambda
#   X  = c(lambda,x+1)
#   v1 = gammainc(X)
#   v2 = gamma(x+1)
#
#   # straight forward formula from the definition of incomplete gamma integral
#   v1_d = -lambda^x*exp(-lambda)
#   v2_d = 0
#
#   z  = v1/v2
#   z_d = (v1_d*v2 - v2_d*v1)/v2^2
#   return(c(z,z_d))
# }


# PF likelihood with resampling for AR(p)
# ParticleFilter_Res_AR_old = function(theta, mod){
#   #--------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling
#   #           to approximate the likelihood of the
#   #           a specified count time series model with an underlying AR(p)
#   #           dependence structure.
#   #
#   #
#   # INPUTS:
#   #    theta: parameter vector
#   #      mod: a list containing all t he information for the model, such as
#   #           count distribution. ARMA model, etc
#   # OUTPUT:
#   #    loglik: approximate log-likelihood
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias
#   # DATE:    July  2020
#   #--------------------------------------------------------------------------#
#
#   # keep track of the random seed to use common random numbers
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # Retrieve parameters ans save them in a li
#   Parms = RetrieveParameters(theta,mod)
#
#   # check for causality
#   if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))
#
#   # Initialize the negative log likelihood computation
#   nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
#                    - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1])))
#
#   # Compute the theoretical covariance for the AR model for current estimate
#   gt    = ARMAacf(ar = Parms$AR, ma = Parms$MA,lag.max = mod$n)
#
#   # Compute the best linear predictor coefficients and errors using Durbin Levinson
#   Phi = list()
#   for (t in 2:mod$n){
#     CurrentDL     = DLAcfToAR(gt[2:t])
#     Phi[[t-1]] = CurrentDL[,1]
#   }
#   Rt           = c(1,sqrt(as.numeric(CurrentDL[,3])))
#
#   # allocate memory for particle weights and the latent Gaussian Series particles, check me: do I weights for 1:T or only 2?
#   w     = matrix(0, mod$n, mod$ParticleNumber)
#   Z     = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)
#
#   #======================   Start the SIS algorithm   ======================#
#   # Initialize the weights
#   w[1,] = rep(1,mod$ParticleNumber)
#
#   # Compute the first integral limits Limit$ a and Limit$b
#   Limit = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#   # Initialize particles from truncated normal distribution
#   #Z[1,] = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#   Z[1,] = SampleTruncNormParticles(mod, Limit, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#
#   # =================== Loop over t ===================== #
#   for (t in 2:mod$n){
#     # Compute the latent Gaussian predictions Zhat_t using Innovations Algorithm
#
#     if(t==2 || mod$ParticleNumber==1){
#       Zhat  =         Z[1:(t-1),] %*% Phi[[t-1]]
#     }else{
#       Zhat  = colSums(Z[1:(t-1),] %*% Phi[[t-1]])
#     }
#
#     # Compute integral limits
#     Limit = ComputeLimits(mod, Parms, t, Zhat, Rt)
#
#     # Sample truncated normal particles
#     #Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b,t, Zhat, Rt)
#     Znew  = SampleTruncNormParticles(mod, Limit, t, Zhat, Rt)
#
#     # update weights
#     #w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])
#     w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])
#
#     # check me: break if I got NA weight
#     if (any(is.na(w[t,]))| sum(w[t,])==0 ){
#       message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
#       return(10^8)
#     }
#
#     # Resample the particles using common random numbers
#     old_state1 = get_rand_state()
#     Znew = ResampleParticles(mod, w, t, Znew)
#     set_rand_state(old_state1)
#
#     # Combine current particles, with particles from previous iterations
#     Z = rbind(Znew, as.matrix(Z[1:(t-1),]))
#     # if (mod$ARMAModel[1]>1){
#     #   Z = rbind(Znew, Z[1:( min(t,max(mod$ARMAModel)) -1),])
#     # }else {
#     #   Z[1,]=Znew
#     # }
#
#     # update log-likelihood
#     nloglik = nloglik - log(mean(w[t,]))
#   }
#
#   return(nloglik)
# }
#
# # PF likelihood with resampling for MA(q)
# ParticleFilter_Res_MA_old = function(theta, mod){
#   #------------------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
#   #           of the a specified count time series model with an underlying MA(1)
#   #           dependence structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paper can be found at:
#   #           https://arxiv.org/abs/1811.00203
#   #           2. This function is very similar to LikSISGenDist_ARp but here
#   #           I have a resampling step.
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             data
#   #    ParticleNumber:   number of particles to be used.
#   #    Regressor:        independent variable
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    July 2020
#   #------------------------------------------------------------------------------------#
#
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # Retrieve parameters and save them in a list called Parms
#   Parms = RetrieveParameters(theta,mod)
#
#   # check for causality
#   if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))
#
#   # Initialize the negative log likelihood computation
#   nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
#                    - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1])))
#
#   # Compute covariance up to lag n-1
#   gt    = as.vector(ARMAacf(ar = Parms$AR, ma = Parms$MA, lag.max = mod$n))
#
#   # Compute coefficients of Innovations Algorithm see 5.2.16 and 5.3.9 in in Brockwell Davis book
#   IA    = innovations.algorithm(gt)
#   Theta = IA$thetas
#   Rt    = sqrt(IA$v)
#
#   # allocate matrices for weights, particles and innovations which are equal to Z-Zhat
#   w     = matrix(0, mod$n, mod$ParticleNumber)
#   Z     = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)
#   Inn   = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)
#
#   # particle filter weights
#   w[1,]   = rep(1,mod$ParticleNumber)
#
#   # Compute the first integral limits Limit$ a and Limit$b
#   Limit = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#   # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#   #Z[1,]   = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#   Z[1,]   = SampleTruncNormParticles(mod, Limit, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#
#   # Compute the first innovation (Zhat_1=0)
#   Inn[1,] = Z[1,]
#
#
#   for (t in 2:mod$n){
#
#     # Compute the latent Gaussian predictions Zhat_t using Innovations Algorithm - see 5.3.9 in Brockwell Davis book
#     if(t==2 || mod$ParticleNumber==1){
#       Zhat  =         Inn[1:(min(t-1,mod$nMA)),] %*% Theta[[t-1]][1:(min(t-1,mod$nMA))]
#     }else{
#       Zhat  = colSums(Inn[1:(min(t-1,mod$nMA)),] %*% Theta[[t-1]][1:(min(t-1,mod$nMA))])
#     }
#
#     # Compute integral limits
#     Limit = ComputeLimits(mod, Parms, t, Zhat, Rt)
#
#
#     # Sample truncated normal particles
#     #Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b, t, Zhat, Rt)
#     Znew  = SampleTruncNormParticles(mod, Limit, t, Zhat, Rt)
#
#     # update weights
#     #w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])
#     w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])
#
#
#     # check me: break if I got NA weight
#     if (any(is.na(w[t,]))| sum(w[t,])==0 ){
#       message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
#       return(10^8)
#     }
#
#     # Resample the particles using common random numbers
#     old_state1 = get_rand_state()
#     Znew = ResampleParticles(mod, w, t, Znew)
#     set_rand_state(old_state1)
#
#     # Compute the new Innovation
#     InnNew = Znew - Zhat
#
#     # Combine current particles, with particles from previous iterations
#     Inn = rbind(InnNew, as.matrix(Inn[1:min(t-1,mod$nMA),]))
#
#     # update likelihood
#     nloglik = nloglik - log(mean(w[t,]))
#
#   }
#
#   # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#   # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#
#   # if (nloglik==Inf | is.na(nloglik)){
#   #   nloglik = 10^8
#   # }
#
#
#   return(nloglik)
# }

# PF likelihood with resampling
# ParticleFilter_Res = function(theta, mod){
#   #--------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling
#   #           to approximate the likelihood of the
#   #           a specified count time series model with an underlying AR(p)
#   #           dependence structure or MA(q) structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paer can be found at:
#   #           https://arxiv.org/abs/1811.00203
#
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             dependent variable
#   #    Regressor:        independent variables
#   #    ParticleNumber:   number of particles to be used in likelihood approximation
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    July 2020
#   #--------------------------------------------------------------------------#
#
#
#   # Pure AR model
#   if(mod$ARMAModel[1]>0 && mod$ARMAModel[2]==0) loglik = ParticleFilter_Res_AR(theta, mod)
#   # Pure MA model or White noise
#   if(mod$ARMAModel[1]==0&& mod$ARMAModel[2]>=0) loglik = ParticleFilter_Res_MA(theta, mod)
#   return(loglik)
# }
#
# # innovations algorithm code - use in previous likelihood implementations
# innovations.algorithm <- function(gamma){
#   n.max=length(gamma)-1
#   # Found this online need to check it
#   # http://faculty.washington.edu/dbp/s519/R-code/innovations-algorithm.R
#   thetas <- vector(mode="list",length=n.max)
#   v <- rep(gamma[1],n.max+1)
#   for(n in 1:n.max){
#     thetas[[n]] <- rep(0,n)
#     thetas[[n]][n] <- gamma[n+1]/v[1]
#     if(n>1){
#       for(k in 1:(n-1)){
#         js <- 0:(k-1)
#         thetas[[n]][n-k] <- (gamma[n-k+1] - sum(thetas[[k]][k-js]*thetas[[n]][n-js]*v[js+1]))/v[k+1]
#       }
#     }
#     js <- 0:(n-1)
#     v[n+1] <- v[n+1] - sum(thetas[[n]][n-js]^2*v[js+1])
#   }
#   v = v/v[1]
#   return(structure(list(v=v,thetas=thetas)))
# }
#
# # PF likelihood with resampling for AR(p) - written more concicely
# ParticleFilter_Res_AR = function(theta, mod){
#   #--------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling
#   #           to approximate the likelihood of the
#   #           a specified count time series model with an underlying AR(p)
#   #           dependence structure.
#   #
#   #
#   # INPUTS:
#   #    theta: parameter vector
#   #      mod: a list containing all t he information for the model, such as
#   #           count distribution. ARMA model, etc
#   # OUTPUT:
#   #    loglik: approximate log-likelihood
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias
#   # DATE:    July  2020
#   #--------------------------------------------------------------------------#
#
#   # keep track of the random seed to use common random numbers
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # Retrieve parameters ans save them in a li
#   Parms = RetrieveParameters(theta,mod)
#
#   # check for causality
#   if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))
#
#   # Initialize the negative log likelihood computation
#   nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
#                    - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1])))
#
#   # Compute the theoretical covariance for the AR model for current estimate
#   gt    = ARMAacf(ar = Parms$AR, ma = Parms$MA,lag.max = mod$n)
#
#   # Compute the best linear predictor coefficients and errors using Durbin Levinson
#   Phi = list()
#   for (t in 2:mod$n){
#     CurrentDL     = DLAcfToAR(gt[2:t])
#     Phi[[t-1]] = CurrentDL[,1]
#   }
#   Rt           = c(1,sqrt(as.numeric(CurrentDL[,3])))
#
#   # allocate memory for particle weights and the latent Gaussian Series particles, check me: do I weights for 1:T or only 2?
#   w     = matrix(0, mod$n, mod$ParticleNumber)
#   Z     = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)
#
#   #======================   Start the SIS algorithm   ======================#
#   # Initialize the weights
#   w[1,] = rep(1,mod$ParticleNumber)
#
#   # Compute the first integral limits Limit$ a and Limit$b
#   Limit = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#   # Initialize particles from truncated normal distribution
#   #Z[1,] = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#   Z[1,] = SampleTruncNormParticles(mod, Limit, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#   # =================== Loop over t ===================== #
#   for (t in 2:mod$n){
#     # Compute the latent Gaussian predictions Zhat_t using Innovations Algorithm
#     Zhat  =         Phi[[t-1]]%*%Z[1:(t-1),]
#
#     # Compute integral limits
#     Limit = ComputeLimits(mod, Parms, t, Zhat, Rt)
#
#     # Sample truncated normal particles
#     #Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b,t, Zhat, Rt)
#     Znew  = SampleTruncNormParticles(mod, Limit, t, Zhat, Rt)
#
#     # update weights
#     #w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])
#     w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])
#
#     # check me: break if I got NA weight
#     if (any(is.na(w[t,]))| sum(w[t,])==0 ){
#       message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
#       return(10^8)
#     }
#
#     # Resample the particles using common random numbers
#     old_state1 = get_rand_state()
#     Znew = ResampleParticles(mod, w, t, Znew)
#     set_rand_state(old_state1)
#
#     # Combine current particles, with particles from previous iterations
#     Z = rbind(Znew, matrix(Z[1:(t-1),],ncol = mod$ParticleNumber))
#
#     # update log-likelihood
#     nloglik = nloglik - log(mean(w[t,]))
#   }
#
#   return(nloglik)
# }
#
# # PF likelihood with resampling for MA(q)
# ParticleFilter_Res_MA = function(theta, mod){
#   #------------------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
#   #           of the a specified count time series model with an underlying MA(1)
#   #           dependence structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paper can be found at:
#   #           https://arxiv.org/abs/1811.00203
#   #           2. This function is very similar to LikSISGenDist_ARp but here
#   #           I have a resampling step.
#   #
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             data
#   #    ParticleNumber:   number of particles to be used.
#   #    Regressor:        independent variable
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    July 2020
#   #------------------------------------------------------------------------------------#
#
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # Retrieve parameters and save them in a list called Parms
#   Parms = RetrieveParameters(theta,mod)
#
#   # check for causality
#   if( CheckStability(Parms$AR,Parms$MA) ) return(10^(8))
#
#   # Initialize the negative log likelihood computation
#   nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
#                    - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1])))
#
#   # Compute covariance up to lag n-1
#   gt    = as.vector(ARMAacf(ar = Parms$AR, ma = Parms$MA, lag.max = mod$n))
#
#   # Compute coefficients of Innovations Algorithm see 5.2.16 and 5.3.9 in in Brockwell Davis book
#   IA    = innovations.algorithm(gt)
#   Theta = IA$thetas
#   Rt    = sqrt(IA$v)
#
#   # allocate matrices for weights, particles and innovations which are equal to Z-Zhat
#   w     = matrix(0, mod$n, mod$ParticleNumber)
#   Z     = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)
#   Inn   = matrix(0, max(mod$ARMAModel), mod$ParticleNumber)
#
#   # particle filter weights
#   w[1,]   = rep(1,mod$ParticleNumber)
#
#   # Compute the first integral limits Limit$ a and Limit$b
#   Limit = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#   # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#   #Z[1,]   = SampleTruncNormParticles(mod, Limit$a, Limit$b, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#   Z[1,]   = SampleTruncNormParticles(mod, Limit, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#   # Compute the first innovation (Zhat_1=0)
#   Inn[1,] = Z[1,]
#
#
#   for (t in 2:mod$n){
#
#     # Compute the latent Gaussian predictions Zhat_t using Innovations Algorithm - see 5.3.9 in Brockwell Davis book
#     if(mod$ParticleNumber==1){
#       if(t==2){
#         Zhat  =         Inn[1:(min(t-1,mod$nMA)),] %*% Theta[[t-1]][1:(min(t-1,mod$nMA))]
#       }else{
#         Zhat  = colSums(Inn[1:(min(t-1,mod$nMA)),] %*% Theta[[t-1]][1:(min(t-1,mod$nMA))])
#       }
#     }else{
#       if(t==2){
#         Zhat  =         Inn[1:(min(t-1,mod$nMA)),] * Theta[[t-1]][1:(min(t-1,mod$nMA))]
#       }else{
#         Zhat  = colSums(Inn[1:(min(t-1,mod$nMA)),] * Theta[[t-1]][1:(min(t-1,mod$nMA))])
#       }
#     }
#
#     # Compute integral limits
#     Limit = ComputeLimits(mod, Parms, t, Zhat, Rt)
#
#     # Sample truncated normal particles
#     #Znew  = SampleTruncNormParticles(mod, Limit$a, Limit$b, t, Zhat, Rt)
#     Znew  = SampleTruncNormParticles(mod, Limit, t, Zhat, Rt)
#
#     # update weights
#     #w[t,] = ComputeWeights(mod, Limit$a, Limit$b, t, w[(t-1),])
#     w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])
#
#
#     # check me: break if I got NA weight
#     if (any(is.na(w[t,]))| sum(w[t,])==0 ){
#       message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
#       return(10^8)
#     }
#
#     # Resample the particles using common random numbers
#     old_state1 = get_rand_state()
#     Znew = ResampleParticles(mod, w, t, Znew)
#     set_rand_state(old_state1)
#
#     # Compute the new Innovation
#     InnNew = Znew - Zhat
#
#     # Combine current particles, with particles from previous iterations
#     Inn[1:min(t,mod$nMA),] = rbind(matrix(InnNew,ncol=mod$ParticleNumber), matrix(Inn[1:min(t-1,mod$nMA-1),],ncol = mod$ParticleNumber))
#
#     # update likelihood
#     nloglik = nloglik - log(mean(w[t,]))
#
#   }
#
#   # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#   # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#
#   # if (nloglik==Inf | is.na(nloglik)){
#   #   nloglik = 10^8
#   # }
#
#
#   return(nloglik)
# }
#
# # older PF likelihood with resampling for ARMA(p,q)
# ParticleFilter_Res_ARMA_old = function(theta, mod){
#   #------------------------------------------------------------------------------------#
#   # PURPOSE:  Use particle filtering with resampling to approximate the likelihood
#   #           of the a specified count time series model with an underlying MA(1)
#   #           dependence structure.
#   #
#   # NOTES:    1. See "Latent Gaussian Count Time Series Modeling" for  more
#   #           details. A first version of the paper can be found at:
#   #           https://arxiv.org/abs/1811.00203
#   #           2. This function is very similar to LikSISGenDist_ARp but here
#   #           I have a resampling step.
#   #           3. The innovations algorithm here is the standard one, likely it will
#   #           be slower.
#   # INPUTS:
#   #    theta:            parameter vector
#   #    data:             data
#   #    ParticleNumber:   number of particles to be used.
#   #    Regressor:        independent variable
#   #    CountDist:        count marginal distribution
#   #    epsilon           resampling when ESS<epsilon*N
#   #
#   # OUTPUT:
#   #    loglik:           approximate log-likelihood
#   #
#   #
#   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#   # DATE:    July 2020
#   #------------------------------------------------------------------------------------#
#
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # Retrieve parameters and save them in a list called Parms
#   Parms = RetrieveParameters(theta,mod)
#   #print(Parms$AR)
#   # check for causality and invertibility
#   if( CheckStability(Parms$AR,Parms$MA) ){
#     mod$ErrorMsg = sprintf('WARNING: The ARMA polynomial must be causal and invertible.')
#     warning(mod$ErrorMsg)
#     return(mod$loglik_BadValue1)
#   }
#
#   # Initialize the negative log likelihood computation
#   nloglik = ifelse(mod$nreg==0,  - log(mod$mypdf(mod$DependentVar[1],Parms$MargParms)),
#                    - log(mod$mypdf(mod$DependentVar[1], Parms$ConstMargParm, Parms$DynamMargParm[1,])))
#
#   # retrieve AR, MA orders and their max
#   m = max(mod$ARMAModel)
#   p = mod$ARMAModel[1]
#   q = mod$ARMAModel[2]
#
#   # Compute ARMA covariance up to lag n-1
#   a        = list()
#   if(!is.null(Parms$AR)){
#     a$phi = Parms$AR
#   }else{
#     a$phi = 0
#   }
#   if(!is.null(Parms$MA)){
#     a$theta = Parms$MA
#   }else{
#     a$theta = 0
#   }
#   a$sigma2 = 1
#   gamma    = itsmr::aacvf(a,mod$n)
#
#   # Compute coefficients of Innovations Algorithm see 5.2.16 and 5.3.9 in in Brockwell Davis book
#   #IA       = InnovAlg(Parms, gamma, mod)
#   IA       = innovations.algorithm(gamma)
#   Rt       = sqrt(IA$v)
#
#   # allocate matrices for weights, particles and predictions of the latent series
#   w        = matrix(0, mod$n, mod$ParticleNumber)
#   Z        = matrix(0, mod$n, mod$ParticleNumber)
#   Zhat     = matrix(0, mod$n, mod$ParticleNumber)
#
#   # initialize particle filter weights
#   w[1,]    = rep(1,mod$ParticleNumber)
#
#   # Compute the first integral limits Limit$a and Limit$b
#   Limit    = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#   # Initialize the particles using N(0,1) variables truncated to the limits computed above
#   Z[1,]    = SampleTruncNormParticles(mod, Limit, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
#
#
#   for (t in 2:mod$n){
#
#     # compute Zhat_t
#     Zhat[t,] = ComputeZhat_t(mod, IA, Z, Zhat,t, Parms)
#
#     # Compute integral limits
#     Limit = ComputeLimits(mod, Parms, t, Zhat[t,], Rt)
#
#     # Sample truncated normal particles
#     Znew  = SampleTruncNormParticles(mod, Limit, t, Zhat[t,], Rt)
#
#     # update weights
#     w[t,] = ComputeWeights(mod, Limit, t, w[(t-1),])
#
#     # check me: break if I got NA weight
#     if (any(is.na(w[t,]))| sum(w[t,])==0 ){
#       #print(t)
#       #print(w[t,])
#       message(sprintf('WARNING: Some of the weights are either too small or sum to 0'))
#       return(mod$loglik_BadValue2)
#     }
#
#     # Resample the particles using common random numbers
#     old_state1 = get_rand_state()
#     Znew = ResampleParticles(mod, w, t, Znew)
#     set_rand_state(old_state1)
#
#     # save the current particle
#     Z[t,]   = Znew
#
#     # update likelihood
#     nloglik = nloglik - log(mean(w[t,]))
#
#   }
#
#   # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#   # nloglik = nloglik- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#
#   # if (nloglik==Inf | is.na(nloglik)){
#   #   nloglik = 10^8
#   # }
#
#
#   return(nloglik)
# }
