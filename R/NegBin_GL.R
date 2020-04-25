NegBin_GL = function(initParam, trueParam, p, q, LB, UB, MaxCdf, nHC, n, nsim, no_cores){

  #---------------------------------------------------------------------------------------------#
  # PURPOSE:        Fit Gaussian likelihood to many realizations of synthetic Poisson data with
  #                 an underlying ARMA series.
  #
  # INPUTS
  #   initParam     initial parameters
  #   trueParam     true parameters
  #   p,q           orders of underlying ARMA model
  #   LB            lower bound for parameters
  #   UB            upper bound for parameters
  #   MaxCdf        cdf will be computed up to this number (for light tails cdf=1 fast)
  #   n             sample size
  #   nsim          number of realizations to be fitted
  #   no_cores      number of cores to be used
  #   nHC           number of HC to be computed
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


  # list with true ARMA parameters
  ARMAmodel = list()
  if(p>0){ARMAmodel[[1]] = trueParam[3:(2+p)]}
  if(q>0){ARMAmodel[[2]] = trueParam[(3+p):length(trueParam)]}

  # Generate all the data and save it in a list
  l <- list()
  for(r in 1:nsim){
    set.seed(r)
    l[[r]] = sim_negbin(n, ARMAmodel, MargParm[1], MargParm[2] )
  }


  # initiate and register the cluster
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  # fit the gaussian log likelihood using foreach
  all = foreach(index = 1:nsim,
                .combine = rbind,
                .packages = c("countsFun")) %dopar%
    FitGaussianLikNB(initParam, l[[index]], LB, UB, c(p,q), MaxCdf, nHC)


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
  df[,4]   = 'gaussianLik'
  df[,5]   = n
  df[,6]   = MAParm
  df[,7]   = all[,6]
  df[,8]   = MargParm[1]
  df[,9]   = all[,4]
  df[,10]  = MargParm[2]
  df[,11]  = all[,5]



  return(df)
}
