# likSISRMA1Old = function(theta, data,ARMAorder,ParticleNumber, CountDist, epsilon){
#
#
#   MargParmIndices = switch(CountDist,
#                            "Poisson"             = 1,
#                            "Negative Binomial"   = 1:2,
#                            "Mixed Poisson"       = 1:3,
#                            "Generalized Poisson" = 1:2,
#                            "Binomial"            = 1:2)
#
#   # retrieve marginal cdf
#   cdf = switch(CountDist,
#                "Poisson"                       = ppois,
#                "Negative Binomial"             = function(x, theta){ pnbinom(q = x, size = theta[1], prob = 1-theta[2]) },
#                "Mixed Poisson"                 = pMixedPoisson,
#                "Generalized Poisson"           = pGenPoisson,
#                "Binomial"                      = pbinom
#   )
#
#   # retrieve marginal pdf
#   pdf = switch(CountDist,
#                "Poisson"                       = dpois,
#                "Negative Binomial"             = function(x, theta){ dnbinom(x, size = theta[1], prob = 1-theta[2]) },
#                "Mixed Poisson"                 = dMixedPoisson,
#                "Generalized Poisson"           = dGenPoisson,
#                "Binomial"                      = dbinom
#   )
#
#   #set.seed(1)
#   theta1.idx = MargParmIndices
#   theta2.idx = ARMAorder[1]
#
#   #set.seed(1)
#   theta1 = theta[theta1.idx]
#   n.theta1.idx = theta1.idx[length(theta1.idx)] # num params in theta1
#   theta2.idx = (n.theta1.idx + 1):(n.theta1.idx + 1)
#   tht = theta[theta2.idx]
#   xt = data
#   T1 = length(xt)
#   N = ParticleNumber # number of particles
#   prt = matrix(0,N,T1) # to collect all particles
#   wgh = matrix(0,N,T1) # to collect all particle weights
#
#   a = qnorm(cdf(xt[1]-1,theta1),0,1)
#   b = qnorm(cdf(xt[1],theta1),0,1)
#   a = rep(a,N)
#   b = rep(b,N)
#   zprev = z.rest(a,b)
#   rt0 = 1+tht^2
#   zhat = tht*zprev/rt0
#   prt[,1] = zhat
#
#
#   nloglik <- 0
#   for (t in 2:T1)
#   {
#     rt0 = 1+tht^2-tht^2/rt0 # This is based on p. 173 in BD book
#     rt = sqrt(rt0/(1+tht^2))
#     a = (qnorm(cdf(xt[t]-1,theta1),0,1) - zhat)/rt
#     b = (qnorm(cdf(xt[t],theta1),0,1) - zhat)/rt
#     err = z.rest(a,b)
#     znew = zhat + rt*err
#     wgh <- pnorm(b,0,1) - pnorm(a,0,1)
#     if (any(is.na(wgh))) # see apf code below for explanation
#     {
#       #nloglik <- NaN
#       nloglik <- Inf
#       break
#     }
#     if (sum(wgh)==0)
#     {
#       #nloglik <- NaN
#       nloglik <- Inf
#       break
#     }
#     wghn <- wgh/sum(wgh)
#     ind <- rmultinom(1, N, wghn)
#     znew <- rep(znew,ind)
#
#     zhat = tht*(znew-zhat)/rt0
#     prt[,t] = zhat
#
#     nloglik <- nloglik -2*log(mean(wgh))
#   }
#
#   nloglik <- nloglik - 2*log(pdf(xt[1],theta1))
#
#
#   out = nloglik
#
#
#   if (out==Inf | is.na(out)){
#     out = 10^8
#   }
#
#   return(out)
# }
