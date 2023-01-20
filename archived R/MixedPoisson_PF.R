# MixedPoisson_PF = function(trueParam, p, q, LB, UB, n, nsim, Particles, epsilon,  no_cores, UseDEOptim, nHC, useTrueInit){
#
#   #---------------------------------------------------------------------------------------------#
#   # PURPOSE:        Fit many realizations of synthetic Mixed Poisson data with an
#   #                 underlying ARMA series using particle filtering.
#   #
#   # INPUTS
#   #   trueParam     true parameters
#   #   p,q           orders of underlying ARMA model
#   #   LB            lower bound for parameters
#   #   UB            upper bound for parameters
#   #   n             sample size
#   #   nsim          number of realizations to be fitted
#   #   Particles     number of particles
#   #   epsilon       resampling when ESS<epsilon*N
#   #   no_cores      number of cores to be used
#   #   nHC           number of HC used for initial estimation
#   #
#   # OUTPUT
#   #   df            parameter estimates, true values, standard errors and likelihood value
#   #
#   # NOTES           no standard errors in cases where maximum is achieved at the boundary
#   #
#   # AUTHORS:        James Livsey, Stefanos Kechagias, Vladas Pipiras
#   #
#   # DATE:           June 2020
#   #
#   # R version 3.6.3
#   #---------------------------------------------------------------------------------------------#
#
#   # ---- Load libraries ----
#   library(parallel)
#   library(doParallel)
#   library(countsFun)
#   library(MixtureInf)
#   library(FitAR)
#
#   # distribution and ARMA order
#   CountDist = "Mixed Poisson"
#   ARMAorder = c(p,q)
#
#   # list with true ARMA parameters
#   ARMAmodel = list(NULL,NULL)
#   if(p>0){ARMAmodel[[1]] = trueParam[4:(3+p)]}
#   if(q>0){ARMAmodel[[2]] = trueParam[(4+p):length(trueParam)]}
#
#   # Generate all the data and save it in a list. Also compute all initial points
#   l <- list()
#   initParam <- list()
#   for(r in 1:nsim){
#     set.seed(r)
#     l[[r]] = sim_mixedPoisson(n, ARMAmodel, MargParm[1], MargParm[2],MargParm[3] )
#     if(!useTrueInit){
#       initParam[[r]] = ComputeInitMixedPoissonAR(l[[r]],n,nHC,LB,UB)
#     }else{
#       initParam[[r]] = trueParam
#     }
#   }
#
#
#   # initiate and register the cluster
#   cl <- makeCluster(no_cores)
#   registerDoParallel(cl)
#
#   # fit the Particle Filter likelihood using foreach
#   all = foreach(index = 1:nsim,
#                 .combine = rbind,
#                 .packages = c("FitAR","countsFun","optimx")) %dopar%
#     FitMultiplePFRes(initParam[[index]], l[[index]], CountDist, Particles, LB, UB, ARMAorder, epsilon)
#
#
#   # Prepare results for the plot.
#   df = data.frame(matrix(ncol = 18, nrow = nsim))
#
#
#   names(df) = c('p.true',
#                 'p.est',
#                 'p.se',
#                 'lam1.true',
#                 'lam1.est',
#                 'lam1.se',
#                 'lam2.true',
#                 'lam2.est',
#                 'lam2.se',
#                 'phi.true',
#                 'phi.est',
#                 'phi.se',
#                 'estim.method',
#                 'n',
#                 'loglik',
#                 'convcode',
#                 'kkt1',
#                 'kkt2')
#
#   # fix me generalize this
#   # true values
#   df[,1] =  trueParam[1]
#   df[,4] =  trueParam[2]
#   df[,7] =  trueParam[3]
#   df[,10] =  trueParam[4]
#
#   # estimated values
#   df[,2] = all[,1]
#   df[,5] = all[,2]
#   df[,8] = all[,3]
#   df[,11] = all[,4]
#
#   # standard errors
#   df[,3] = all[,5]
#   df[,6] = all[,6]
#   df[,9] = all[,7]
#   df[,12] = all[,8]
#
#
#   df[,13]   = 'particle'
#   df[,14]   = n
#   df[,15]   = all[,9]
#   df[,16]   = all[,10]
#   df[,17]   = all[,11]
#   df[,18]   = all[,12]
#
#
#
#
#   return(df)
# }
#
#
#
#
# # i = 60
# # x0 = initParam[[i]]
# # X  = l[[i]]
# # ARMAorder = c(1,0)
# # ParticleNumber = Particles
# #
# # optim.output <- optimx(par           = x0,
# #                        fn             = ParticleFilterRes,
# #                        data           = X,
# #                        ARMAorder      = ARMAorder,
# #                        ParticleNumber = ParticleNumber,
# #                        CountDist      = CountDist,
# #                        epsilon        = epsilon,
# #                        lower          = LB,
# #                        upper          = UB,
# #                        hessian        = TRUE,
# #                        method         = "L-BFGS-B",
# #                        control        = list(dowarn = FALSE))
# #
# #
# # ParmEst = c(optim.output$p1,optim.output$p2,optim.output$p3,optim.output$p4)
# # loglik   = optim.output$value
# #
# #
# # H = gHgen(par       = ParmEst,
# #           fn        = ParticleFilterRes,
# #           data      = X,
# #           ARMAorder = ARMAorder,
# #           CountDist      = CountDist,
# #           ParticleNumber = ParticleNumber,
# #           epsilon        = epsilon
# # )
# # ifelse(H$hessOK,H$Hn,NA)
#
#
#
#
#
#
#
#
#
#
# # i = 60
# # x0 = initParam[[i]]
# # X  = l[[i]]
# # ARMAorder = c(1,0)
# # ParticleNumber = Particles
# # ParticleFilterRes(x0, X, ARMAorder, ParticleNumber, CountDist, epsilon)
# # # #
# # # #
# # # FitMultiplePFRes(initParam[[1]], l[[1]], CountDist, Particles, LB, UB, ARMAorder, epsilon)
# # optim.output <- optimx(par            = x0,
# #                        fn             = ParticleFilterRes,
# #                        data           = X,
# #                        ARMAorder      = ARMAorder,
# #                        ParticleNumber = ParticleNumber,
# #                        CountDist      = CountDist,
# #                        epsilon        = epsilon,
# #                        lower          = LB,
# #                        upper          = UB,
# #                        method         = "L-BFGS-B",
# #                        control        = (dowarn = FALSE))
#
# #
# # ParmEst = c(optim.output$p1,optim.output$p2,optim.output$p3,optim.output$p4)
# # loglik   = optim.output$value
# #
# # # compute hessian
# # H = gHgen(par       = ParmEst,
# #           fn        = ParticleFilterRes,
# #           data      = X,
# #           ARMAorder = ARMAorder,
# #           ParticleNumber = ParticleNumber,
# #           CountDist      = CountDist,
# #           epsilon        = epsilon
# # )$Hn
# #
