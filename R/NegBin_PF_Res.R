NegBin_PF_Res = function(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon,  no_cores, nHC, useTrueInit){

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
  #library(optimx)
  # add nloptr package
  #library(nloptr)
  #library(numDeriv)

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
      initParam[[r]] = ComputeInitNegBinMA(l[[r]],n,nHC, LB, UB)
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
                .packages = c("FitAR","countsFun", "optimx")) %dopar%
    FitMultiplePFMA1Res(initParam[[index]], l[[index]], CountDist, Particles, LB, UB, ARMAorder, epsilon)


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


# # select series that has an issue
# i=1
# x0  = initParam[[i]]
# X = l[[i]]
# ParticleNumber = Particles
#
#
# # this is what I currently run but with optim instead of optimx
# optim.output <- optim(par            = x0,
#                       fn             = ParticleFilterMA1_Res,
#                       data           = X,
#                       ARMAorder      = ARMAorder,
#                       ParticleNumber = ParticleNumber,
#                       CountDist      = CountDist,
#                       epsilon        = epsilon,
#                       lower          = LB,
#                       upper          = UB,
#                       hessian        = TRUE,
#                       method         = "L-BFGS-B")
#
# # get parameter estimates
# ParmEst  = c(optim.output$p1,optim.output$p2,optim.output$p3)
#
# # get Hessian (different than optim! although the estimates are the same)
# H = gHgen(par            = ParmEst,
#           fn             = ParticleFilterMA1_Res,
#           data           = X,
#           ARMAorder      = ARMAorder,
#           CountDist      = CountDist,
#           ParticleNumber = ParticleNumber,
#           epsilon        = epsilon
# )
#
#
#
#
# # wrapper for my function and the derivative
# myfun = function(theta)ParticleFilterMA1_Res(theta, X, ARMAorder, ParticleNumber, CountDist, 1)
# mygrad = function(theta)grad(myfun, theta)
#
#
# #=============================================================================================#
# # nloptr version of what I am doing
# x0 = trueParam
# nlopt = nloptr(x0     = x0 ,
#        eval_f = myfun,
#        eval_grad_f = mygrad,
#            lb = LB,
#            ub = UB,
#          opts = list("algorithm" = "NLOPT_LD_LBFGS_NOCEDAL",
#                      "xtol_rel"  = 1.0e-16))
#
#
# #=============================================================================================#
#
# # nloptr NLOPT_LN_BOBYQA solver
# x0 = trueParam
# UB[1] =10
# nlopt = nloptr(x0     = x0 ,
#                eval_f = myfun,
#                lb = LB,
#                ub = UB,
#                opts = list("algorithm" = "NLOPT_LN_BOBYQA",
#                            "ftol_rel"  = 1.0e-16))
#
#
#
# #=============================================================================================#
# x0 = trueParam
# UB[1] =10
# nlopt = nloptr(x0     = x0 ,
#                eval_f = myfun,
#                lb = LB,
#                ub = UB,
#                opts = list("algorithm" = "NLOPT_LN_NEWUOA_BOUND",
#                            "xtol_rel"  = 1.0e-16))
#
# # Call:
# #   nloptr(x0 = x0, eval_f = myfun, eval_grad_f = mygrad, lb = LB,     ub = UB, opts = list(algorithm = "NLOPT_LN_NEWUOA_BOUND",         xtol_rel = 1e-16))
# #
# #
# # Minimization using NLopt version 2.4.2
# #
# # NLopt solver status: 5 ( NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached. )
# #
# # Number of Iterations....: 100
# # Termination conditions:  xtol_rel: 1e-16
# # Number of inequality constraints:  0
# # Number of equality constraints:    0
# # Current value of objective function:  1229.86345482372
# # Current value of controls: 3.006447 0.7998288 -0.8219214
#
#
#
# #=============================================================================================#
# # starting from the truth this gives me a better value!
# x0 = trueParam
# nlopt = nloptr(x0     = x0 ,
#                eval_f = ParticleFilterMA1_Res,
#                data = X,
#                ARMAorder      = ARMAorder,
#                ParticleNumber = ParticleNumber,
#                CountDist      = CountDist,
#                epsilon        = epsilon,
#                lb = LB,
#                ub = UB,
#                opts = list("algorithm" = "NLOPT_LN_COBYLA",
#                            "xtol_rel"  = 1.0e-16))

# Call:
#   nloptr(x0 = x0, eval_f = ParticleFilterMA1_Res, lb = LB, ub = UB,     opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-16),
#          data = X, ARMAorder = ARMAorder, ParticleNumber = ParticleNumber,     CountDist = CountDist, epsilon = epsilon)
#
#
# Minimization using NLopt version 2.4.2
#
# NLopt solver status: 5 ( NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached. )
#
# Number of Iterations....: 100
# Termination conditions:  xtol_rel: 1e-16
# Number of inequality constraints:  0
# Number of equality constraints:    0
# Current value of objective function:  1231.41322889749
# Current value of controls: 3.013989 0.7991799 -0.7638762
