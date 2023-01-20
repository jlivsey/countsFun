# MixedPoisson_IYW = function(trueParam, p, q, LB, UB, n, nsim, nHC){
#
#   #---------------------------------------------------------------------------------------------#
#   # PURPOSE:        Fit many realizations of synthetic Mixed Poisson data with an
#   #                 underlying ARMA series using IYW.
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
#
#   # Prepare results for the plot.
#   df = data.frame(matrix(ncol = 14, nrow = nsim))
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
#                 'n')
#
#   ParamEst = matrix(ncol=4, nrow = nsim)
#   # Generate data and fit the model
#   for(r in 1:nsim){
#     set.seed(r)
#     l = sim_mixedPoisson(n, ARMAmodel, MargParm[1], MargParm[2],MargParm[3] )
#     ParamEst[r,] = ComputeInitMixedPoissonAR(l,n,nHC,LB,UB)
#   }
#
#
#
#   # fix me generalize this
#   # true values
#   df[,1]  =  trueParam[1]
#   df[,4]  =  trueParam[2]
#   df[,7]  =  trueParam[3]
#   df[,10] =  trueParam[4]
#
#   # estimated values
#   df[,2]  = ParamEst[,1]
#   df[,5]  = ParamEst[,2]
#   df[,8]  = ParamEst[,3]
#   df[,11] = ParamEst[,4]
#
#
#   df[,13]   = 'IYW'
#   df[,14]   = n
#
#
#   return(df)
# }
#
# #
# # i = 1
# # x0 = initParam[[i]]
# # X  = l[[i]]
# # ARMAorder = c(1,0)
# # ParticleNumber = Particles
# # # ParticleFilterRes(theta, data, ARMAorder, ParticleNumber, CountDist, epsilon)
# # #
# # #
# # FitMultiplePFRes(initParam[[1]], l[[1]], CountDist, Particles, LB, UB, ARMAorder, epsilon)
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
# #                        method         = "L-BFGS-B")
# #
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
