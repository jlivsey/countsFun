NegBin_PF = function(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon,  no_cores, UseDEOptim, nHC, useTrueInit){

  #---------------------------------------------------------------------------------------------#
  # PURPOSE:        Fit many realizations of synthetic Negative Binomial data with an
  #                 underlying ARMA series using particle filtering.
  #
  # INPUTS
  #   trueParam     true parameters
  #   p,q           orders of underlying ARMA model
  #   LB            lower bound for parameters
  #   UB            upper bound for parameters
  #   n             sample size
  #   nsim          number of realizations to be fitted
  #   Particles     number of particles
  #   epsilon       resampling when ESS<epsilon*N
  #   no_cores      number of cores to be used
  #   nHC           number of HC used for initial estimation
  #
  # OUTPUT
  #   df            parameter estimates, true values, standard errors and likelihood value
  #
  # NOTES           no standard errors in cases where maximum is achieved at the boundary
  #
  # AUTHORS:        James Livsey, Stefanos Kechagias, Vladas Pipiras
  #
  # DATE:           April 2020
  #
  # R version 3.6.3
  #---------------------------------------------------------------------------------------------#

  # ---- Load libraries ----
  library(parallel)
  library(doParallel)
  library(countsFun)

  # distribution and ARMA order
  CountDist = "Negative Binomial"
  ARMAorder = c(p,q)

    # list with true ARMA parameters
  ARMAmodel = list(NULL,NULL)
  if(p>0){ARMAmodel[[1]] = trueParam[3:(2+p)]}
  if(q>0){ARMAmodel[[2]] = trueParam[(3+p):length(trueParam)]}

  # Generate all the data and save it in a list. Also compute all initial points
  l <- list()
  initParam <- list()
  for(r in 1:nsim){
    set.seed(r)
    l[[r]] = sim_negbin(n, ARMAmodel, MargParm[1], MargParm[2] )
    if(!useTrueInit){
      initParam[[r]] = ComputeInitNegBinMA(l[[r]],n,nHC)
    }else{
      initParam[[r]] = trueParam
    }
  }


  # initiate and register the cluster
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  # fit the Particle Filter likelihood using foreach
  all = foreach(index = 1:nsim,
                .combine = rbind,
                .packages = c("FitAR","countsFun")) %dopar%
    FitMultiplePFMA1New(initParam[[index]], l[[index]], CountDist, Particles, LB, UB, ARMAorder, epsilon, UseDEOptim)


  stopCluster(cl)

  # Prepare results for the plot.
  df = data.frame(matrix(ncol = 11, nrow = nsim))


  names(df) = c('r.est',
                'p.est',
                'theta.est',
                'estim.method',
                'n',
                'theta.true',
                'theta.se',
                'r.true',
                'r.se',
                'p.true',
                'p.se' )

  df[,1:3] = all[,1:3]
  df[,4]   = 'particle'
  df[,5]   = n
  df[,6]   = trueParam[3]
  df[,7]   = all[,6]
  df[,8]   = MargParm[1]
  df[,9]   = all[,4]
  df[,10]  = MargParm[2]
  df[,11]  = all[,5]



  return(df)
}
