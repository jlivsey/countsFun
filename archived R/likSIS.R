# #-----likSIS function used in initial implementation of LCG for initial paper version#
# likSIS = function(theta, data,  ParticleNumber, CountDist){
#   old_state <- get_rand_state()
#   on.exit(set_rand_state(old_state))
#
#   # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments
#
#   # retrieve indices of marginal distribution parameters
#   MargParmIndices = switch(CountDist,
#                            "Poisson"             = 1,
#                            "Negative Binomial"   = 1:2,
#                            "Mixed Poisson"       = 1:3,
#                            "Generalized Poisson" = 1:2,
#                            "Binomial"            = 1:2)
#
#   # retrieve marginal cdf
#   cdf = switch(CountDist,
#                  "Poisson"             = ppois,
#                  "Negative Binomial"   = function(x, theta){ pnbinom (q = x, size = theta[1], prob = 1-theta[2])},
#                  "Mixed Poisson"       = function(x, theta){ pmixpois(x, p = theta[1], lam1 = theta[2], lam2 = theta[3])},
#                  "Generalized Poisson" = pGenPoisson,
#                  "Binomial"            = pbinom
#   )
#
#   # retrieve marginal pdf
#   pdf = switch(CountDist,
#                  "Poisson"             = dpois,
#                  "Negative Binomial"   = function(x, theta){ dnbinom (x, size = theta[1], prob = 1-theta[2]) },
#                  "Mixed Poisson"       = function(x, theta){ dmixpois(x, p = theta[1], lam1 = theta[2], lam2 = theta[3])},
#                  "Generalized Poisson" = dGenPoisson,
#                  "Binomial"            = dbinom
#   )
#
#
#   # retrieve marginal distribution parameters
#   # retrieve marginal distribution parameters
#   MargParms  = theta[MargParmIndices]
#   nMargParms = length(MargParms)
#   nparms     = length(theta)
#
#
#   # retrieve ARMA parameters
#   if(ARMAorder[1]>0){
#     AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
#   }else{
#     AR = NULL
#   }
#
#   if(ARMAorder[2]>0){
#     MA = theta[ (length(theta) - ARMAorder[2]) : length(theta)]
#   }else{
#     MA = NULL
#   }
#
#   # retrieve AR order
#   ARorder = ARMAorder[1]
#   phi = AR
#   theta1 = MargParms
#
#
#   xt = data
#   T1 = length(xt)
#   prt = matrix(0,N,T1) # to collect all particles
#   wgh = matrix(0,N,T1) # to collect all particle weights
#
#   a = qnorm(cdf(xt[1]-1,theta1),0,1)
#   b = qnorm(cdf(xt[1],theta1),0,1)
#   a = rep(a,N)
#   b = rep(b,N)
#   zprev = z.rest(a,b)
#   zhat = phi*zprev
#   prt[,1] = zhat
#
#   wprev = rep(1,N)
#   wgh[,1] = wprev
#
#   for (t in 2:T1)
#   {
#     rt = sqrt(1-phi^2)
#     a = (qnorm(cdf(xt[t]-1,theta1),0,1) - phi*zprev)/rt
#     b = (qnorm(cdf(xt[t],theta1),0,1) - phi*zprev)/rt
#     err = z.rest(a,b)
#     znew = phi*zprev + rt*err
#     zhat = phi*znew
#     prt[,t] = zhat
#     zprev = znew
#
#     wgh[,t] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
#     wprev = wgh[,t]
#   }
#
#   lik = pdf(xt[1],theta1)*mean(wgh[,T1])
#   nloglik = (-2)*log(lik)
#
#   out = if (is.na(nloglik)) Inf else nloglik
#   return(out)
# }
