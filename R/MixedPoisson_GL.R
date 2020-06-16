MixedPoisson_GL = function(trueParam, p, q, LB, UB, MaxCdf, nHC, n, nsim, no_cores){

  #---------------------------------------------------------------------------------------------#
  # PURPOSE:        Fit Gaussian likelihood to many realizations of synthetic MixPois data with
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
  library(optimx)

  # list with true ARMA parameters
  ARMAmodel = list(NULL,NULL)
  if(p>0){ARMAmodel[[1]] = trueParam[4:(3+p)]}
  if(q>0){ARMAmodel[[2]] = trueParam[(4+p):length(trueParam)]}

  # Generate all the data and save it in a list
  l <- list()
  initParam <- list()
  for(r in 1:nsim){
    set.seed(r)
    l[[r]] = sim_mixedPoisson(n, ARMAmodel, MargParm[1], MargParm[2], MargParm[3] )
    # initParam[[r]] = ComputeInitNegBinMA(l[[r]],n,nHC)
    initParam[[r]] = trueParam
  }

  # initiate and register the cluster
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  # fit the gaussian log likelihood using foreach
  all = foreach(index = 1:nsim,
                .combine = rbind,
                .packages = c("countsFun","optimx")) %dopar%
    FitGaussianLikMP(initParam[[index]], l[[index]], LB, UB, c(p,q), MaxCdf, nHC)


  stopCluster(cl)

  # Prepare results for the plot.
  df = data.frame(matrix(ncol = 14, nrow = nsim))


  names(df) = c('p.true',
                'p.est',
                'p.se',
                'lam1.true',
                'lam1.est',
                'lam1.se',
                'lam2.true',
                'lam2.est',
                'lam2.se',
                'phi.true',
                'phi.est',
                'phi.se',
                'estim.method',
                'n')

  # fix me generalize this
  # true values
  df[,1] =  trueParam[1]
  df[,4] =  trueParam[2]
  df[,7] =  trueParam[3]
  df[,10] =  trueParam[4]

  # estimated values
  df[,2] = all[,1]
  df[,5] = all[,2]
  df[,8] = all[,3]
  df[,11] = all[,4]

  # standard errors
  df[,3] = all[,5]
  df[,6] = all[,6]
  df[,9] = all[,7]
  df[,12] = all[,8]


  df[,13]   = 'gaussianLik'
  df[,14]   = n

  return(df)
}



#
#   for(index in 1 :3){
#     FitGaussianLikMP(initParam[[index]], l[[index]], LB, UB, c(p,q), MaxCdf, nHC)
#   }
#
#   FitGaussianLikMP(initParam[[34]], l[[34]], LB, UB, c(p,q), MaxCdf, nHC)
#
#
#
i = 181
  optimx.output <- optim(par       = initParam[[i]],
                        fn        = GaussLogLikMP,
                        data      = l[[i]],
                        ARMAorder = c(p,q),
                        MaxCdf    = MaxCdf,
                        nHC       = nHC,
                        method    = "L-BFGS-B",
                        hessian   = TRUE,
                        lower     = LB,
                        upper     = UB
  )
#
#   # save estimates, loglik and standard errors
#   ParmEst  = c(optimx.output$p1,optimx.output$p2,optimx.output$p3,optimx.output$p4)
#
#   ARMAorder = c(1,0)
#   # compute hessian
#   H = gHgen(par   = ParmEst,
#             fn    = GaussLogLikMP,
#             data  = l[[i]],
#             ARMAorder = ARMAorder,
#             MaxCdf    = MaxCdf,
#             nHC       = nHC
#   )$Hn
#
#   se      = sqrt(abs(diag(solve(H))))
