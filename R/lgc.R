#---------------------------------------------------------------------------------------------------------#
# PURPOSE: Main Wrapper for the lgc package.
#
#
#
#
# Authors: Stefanos Kechagias, Jiajie Kong, James Livsey, Robert Lund, Vladas Pipiras
#---------------------------------------------------------------------------------------------------------#

#' Main Wrapper for Fitting Latent Gaussian Count Time Series Models
#'
#' Fits latent Gaussian count models using particle filtering methods for various count distributions
#' (e.g., Poisson, Negative Binomial, Zero-Inflated Poisson) and optional ARMA dependencies and regressors.
#' The function supports tasks such as model evaluation, optimization, and simulation.
#'
#' @param formula An object of class \code{formula} (e.g., \code{y ~ x1 + x2}). Specifies the response and regressors.
#' @param data A data frame containing the variables used in the formula.
#' @param EstMethod Character. Estimation method to use. Default is \code{"PFR"} (particle filter with resampling).
#' @param CountDist Character. Distribution of the count variable. Supported: \code{"Poisson"}, \code{"Negative Binomial"}, \code{"ZIP"}, etc.
#' @param ARMAModel Numeric vector of length 2 specifying AR and MA orders (e.g., \code{c(1,1)} for ARMA(1,1)).
#' @param ParticleNumber Number of particles to use in the particle filter.
#' @param epsilon Numeric. Smoothing parameter for the resampling step (if applicable).
#' @param initialParam Optional. Numeric vector of initial parameter values. If \code{NULL}, estimates are generated internally.
#' @param TrueParam Optional. True parameter vector used when \code{Task == "Simulation"}.
#' @param Task Character. Specifies the task to perform: \code{"Evaluation"}, \code{"Optimization"}, or \code{"Simulation"}.
#' @param SampleSize Integer. Used only when \code{Task == "Simulation"}.
#' @param nsim Integer. Number of simulations to run (if \code{Task == "Simulation"}).
#' @param no_cores Number of CPU cores to use in parallel computing (only for simulation).
#' @param OptMethod Optimization method passed to \code{optimx} (e.g., \code{"L-BFGS-B"}).
#' @param OutputType Output format: \code{"list"} (default) or \code{"wide"}.
#' @param ParamScheme Optional. Parameterization scheme name (if applicable).
#' @param maxdiff Numeric. Convergence threshold for the Innovations algorithm.
#' @param ntrials Integer. Required if \code{CountDist == "Binomial"} to specify number of trials.
#' @param verbose Logical. If \code{TRUE} (default), informative messages are printed during execution.
#' Set to \code{FALSE} to suppress messages.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' If \code{Task} is \code{"Evaluation"} or \code{"Optimization"}, returns a list of class \code{"lgc"} containing:
#' \itemize{
#'   \item \code{ParamEstimates}, \code{StdErrors}, \code{FitStatistics}, \code{OptimOutput}
#'   \item \code{residuals}: Model residuals
#'   \item \code{Model}, \code{Task}, \code{EstMethod}, \code{SampleSize}, etc.
#' }
#'
#' If \code{Task == "Simulation"}, returns a list of \code{lgc} objects (length = \code{nsim}).
#'
#' @details
#' This function serves as the main interface for the \code{lgc} package. It parses inputs, prepares model specifications,
#' fits the model via particle filtering with optional optimization, and computes residuals using results from
#' \href{https://doi.org/10.1080/01621459.2021.1944874}{Jia et al. (2021)}.
#'
#' @references
#' Jia, Y., Kechagias, S., Livsey, J., Lund, R., & Pipiras, V. (2021). Latent Gaussian Count Time Series.
#' \emph{Journal of the American Statistical Association}, 118(541), 596–606.
#' \doi{10.1080/01621459.2021.1944874}
#'
#' @examples
#' CountDist      = "Generalized Poisson"
#' alpha          = 1
#' b0             = 0.5
#' b1             = 2
#' MargParm       = c(b0,b1,alpha)
#' ARParm         = 0.75
#' MAParm         = NULL
#' ARMAModel      = c(length(ARParm),length(MAParm))
#' SampleSize     = 50
#' Regressor      = rbinom(SampleSize,1,0.1)
#' Intercept      = TRUE
#'
#' # simulate data with the old lgc function
#' set.seed(1)
#' DependentVar  = sim_lgc(SampleSize, CountDist, MargParm, ARParm, MAParm, Regressor,Intercept)
#'
#' # save the data in a data frame
#' df = data.frame(DependentVar, Regressor)
#'
#' # specify the regression model
#' formula = DependentVar~Regressor
#'
#' # specify task
#' Task = "Evaluation"
#'
#' # specify parameters to evaluate the log-likelihood
#' initialParam = c(MargParm, ARParm, MAParm)
#'
#' # call the wrapper function with less arguments
#' mylgc = lgc(formula   = formula,
#'             data      = df,
#'             CountDist = CountDist,
#'             ARMAModel = ARMAModel,
#'            Task = Task,
#'             initialParam = initialParam)
#'
#' # check the results
#' summary(mylgc)
#' @export
# Final wrapper function
lgc = function(formula        = NULL,
               data           = NULL,
               EstMethod      = "PFR",
               CountDist      = NULL,
               ARMAModel      = NULL,
               ParticleNumber = 5,
               epsilon        = 0.5,
               initialParam   = NULL,
               TrueParam      = NULL,
               Task           = 'Evaluation',
               SampleSize     = NULL,
               nsim           = NULL,
               no_cores       = 1,
               OptMethod      = "L-BFGS-B",
               OutputType     = "list",
               ParamScheme    = NULL,
               maxdiff        = 10^(-8),
               ntrials        = NULL,
               verbose        = TRUE,...){

  # parse the regression formula
  parsed_formula <- parse_formula(formula)

  # if task is Simulation or Synthesis the data frame will not have a dependent variable
  if (Task %in% c("Simulation", "Synthesis")){
    DependentVar = NULL
  }else{
    # retrieve the Dependent variable
    DependentVar = data[parsed_formula$DependentVar]
  }

  # retrieve the Regressors variable
  if(is.null(parsed_formula$Regressor)){
    Regressor = NULL
  } else{
    Regressor =   data[parsed_formula$Regressor]
  }

  # retrieve intercept
  Intercept = parsed_formula$intercept

  # add a column of ones in the Regressors if Intercept is present
  if (!is.null(Regressor) && Intercept){
    Regressor = cbind(rep(1,dim(data)[1]),Regressor)
    names(Regressor)[1] = "Intercept"
  }

  # parse all the parameters and the data into a list called mod
  mod = ModelScheme(DependentVar   = DependentVar,
                    Regressor      = Regressor,
                    Intercept      = Intercept,
                    EstMethod      = EstMethod,
                    ARMAModel      = ARMAModel,
                    CountDist      = CountDist,
                    ParticleNumber = ParticleNumber,
                    epsilon        = epsilon,
                    initialParam   = initialParam,
                    TrueParam      = TrueParam,
                    Task           = Task,
                    SampleSize     = SampleSize,
                    OptMethod      = OptMethod,
                    OutputType     = OutputType,
                    ParamScheme    = ParamScheme,
                    maxdiff        = maxdiff,
                    ntrials        = ntrials,
                    verbose        = verbose,
                    nsim           = nsim)


  # if simulation task has been chosen simulate the data and compute initial estimates
  # check me how fast is this?
  if(Task=='Simulation'){

    # retrieve the parameters
    Parms = RetrieveParameters(TrueParam,mod)

    AllSimulatedSeries <- vector(mode='list', length=nsim)
    AllInitialParam    <- vector(mode='list', length=nsim)
    for (i in 1:nsim) {
      set.seed(i)
      AllSimulatedSeries[[i]] = mod$DependentVar =sim_lgc(SampleSize, mod$CountDist, Parms$MargParms, Parms$AR, Parms$MA,
                                                          mod$Regressor,mod$Intercept)
      AllInitialParam[[i]]    = InitialEstimates(mod)
    }

    # renew the task (after we simulated we need to fit)
    mod$Task = 'Optimization'

    # we need some cores in order to fit the data - we ll use all but one
    if(is.null(no_cores)) no_cores = detectCores() - 1

    # initiate and register the cluster
    cl <- makeCluster(no_cores)

    #clusterSetRNGStream(cl, 1001) #make the bootstrapping exactly the same as above to equate computation time
    registerDoParallel(cl)

    # run foreach
    # fix me: need to be very careful here with packages - and run tests for all distributions with and without regressors

    # .packages = c("ltsa", "optimx", 'tictoc', 'countsFun', 'itsmr',
    #               'doParallel','numDeriv','VGAM','iZID','extraDistr','devtools',
    #               'parallel','MASS','mixtools', 'optextras')

    SimResults = foreach(ForEachIndex = 1:nsim,
                .packages = c("optimx", 'countsFun'),.export= c("FitMultiplePF_Res"))  %dopar%  {
                  mod$DependentVar =  AllSimulatedSeries[[ForEachIndex]]
                  theta  = mod$initialParam = AllInitialParam[[ForEachIndex]]
                  #FitMultiplePF_Res(theta,mod)

                  # Fit the model
                  fit_result <- FitMultiplePF_Res(theta, mod)

                  # Attach initial parameters to the result
                  fit_result$initialParam <- theta

                  # Return the full object
                  fit_result
                }

    stopCluster(cl)

    # revert the task back to Simulation
    mod$Task = 'Simulation'

  }

  if(Task %in% c('Evaluation', 'Optimization')){
    # compute initial parameters if they haven't been provided
    if (is.null(mod$initialParam)){
      mod$initialParam = InitialEstimates(mod)
    }
    theta  = mod$initialParam
    FitResults = FitMultiplePF_Res(theta, mod)

    # gather the input information and the Fit Results in one output structure
    out = PrepareOutput(mod, FitResults)

    # compute residuals
    resid = ComputeResiduals(out)
    out$residuals = resid

    # create an lgc object
    class(out) <- "lgc"
  }

  if(Task=='Simulation'){
    # gather the input information and the Fit Results in one output structure
    out = PrepareOutput(mod, SimResults)

  }


return(out)

}

