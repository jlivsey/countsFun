#' #' Title
#' #'
#' #' @param theta
#' #' @param data
#' #'
#' #' @return
#' #' @export
#' #'
#' #'
#' LikSISGenDist_ARp = function(theta, data, ARMAorder, ParticleNumber, CountDist){
#'   ##########################################################################
#'   # PURPOSE:  Use particle filtering to approximate the likelihood of the
#'   #           a specified count time series model with an underlying AR(p)
#'   #           dependence structure.
#'   #
#'   # NOTES:    See "Latent Gaussian Count Time Series Modeling" for  more
#'   #           details. A first version of the paer can be found at:
#'   #           https://arxiv.org/abs/1811.00203
#'   #
#'   # INPUTS:
#'   #    theta:            parameter vector
#'   #    data:             data
#'   #    ParticleNumber:   number of particles to be used.
#'   #    CountDist:        count marginal distribution
#'   #
#'   # OUTPUT:
#'   #    loglik:           approximate log-likelihood
#'   #
#'   #
#'   # AUTHORS: James Livsey, Vladas Pipiras, Stefanos Kechagias,
#'   # DATE:    November 2019
#'   ##########################################################################
#'
#'   old_state <- get_rand_state()
#'   on.exit(set_rand_state(old_state))
#'
#'   # FIX ME: When I add the function in LGC the following can happen in LGC and pass here as arguments
#'
#'   # retrieve indices of marginal distribution parameters
#'   MargParmIndices = switch(CountDist,
#'                            "Poisson"             = 1,
#'                            "Negative Binomial"   = 1:2,
#'                            "Mixed Poisson"       = 1:3,
#'                            "Generalized Poisson" = 1:2,
#'                            "Binomial"            = 1:2)
#'
#'   # retrieve marginal cdf
#'   mycdf = switch(CountDist,
#'                  "Poisson"             = ppois,
#'                  "Negative Binomial"   = function(x, theta){ pnbinom (q = x, size = theta[1], prob = 1-theta[2])},
#'                  "Mixed Poisson"       = function(x, theta){ pmixpois(x, p = theta[1], lam1 = theta[2], lam2 = theta[3])},
#'                  "Generalized Poisson" = pGenPoisson,
#'                  "Binomial"            = pbinom
#'   )
#'
#'   # retrieve marginal pdf
#'   mypdf = switch(CountDist,
#'                  "Poisson"             = dpois,
#'                  "Negative Binomial"   = function(x, theta){ dnbinom (x, size = theta[1], prob = 1-theta[2]) },
#'                  "Mixed Poisson"       = function(x, theta){ dmixpois(x, p = theta[1], lam1 = theta[2], lam2 = theta[3])},
#'                  "Generalized Poisson" = dGenPoisson,
#'                  "Binomial"            = dbinom
#'   )
#'
#'
#'   # retrieve marginal distribution parameters
#'   MargParms  = theta[MargParmIndices]
#'   nMargParms = length(MargParms)
#'   nparms     = length(theta)
#'
#'
#'   # retrieve ARMA parameters
#'   if(ARMAorder[1]>0){
#'     AR = theta[(nparms-ARMAorder[1]+1):(nMargParms + ARMAorder[1])  ]
#'   }else{
#'     AR = NULL
#'   }
#'
#'   if(ARMAorder[2]>0){
#'     MA = theta[ (length(theta) - ARMAorder[2]) : length(theta)]
#'   }else{
#'     MA = NULL
#'   }
#'
#'   # retrieve AR order
#'   ARorder = ARMAorder[1]
#'
#'
#'   if (prod(abs(polyroot(c(1,-AR))) > 1)){ # check if the ar model is causal
#'
#'     xt = data
#'     T1 = length(xt)
#'     N = ParticleNumber          # number of particles
#'     prt = matrix(0,N,T1)        # to collect all particles
#'     wgh = matrix(0,T1,N)        # to collect all particle weights
#'
#'     # allocate memory for zprev
#'     ZprevAll = matrix(0,p,N)
#'
#'     # Compute integral limits
#'     a = rep( qnorm(mycdf(xt[1]-1,t(theta1)),0,1), N)
#'     b = rep( qnorm(mycdf(xt[1],t(theta1)),0,1), N)
#'
#'     # Generate N(0,1) variables restricted to (ai,bi),i=1,...n
#'     zprev = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#'
#'     # save the currrent normal variables
#'     ZprevAll[1,] = zprev
#'
#'     # initial estimate of first AR coefficient as Gamma(1)/Gamma(0) and corresponding error
#'     phit = TacvfAR(AR)[2]/TacvfAR(AR)[1]
#'     rt = as.numeric(sqrt(1-phit^2))
#'
#'     # particle filter weights
#'     wprev = rep(1,N)
#'     wgh[1,] = wprev
#'
#'     #t0 = proc.time()
#'     # First p steps:
#'     if (p>=2){
#'       for (t in 2:p){
#'
#'         # best linear predictor is just Phi1*lag(Z,1)+...+phiP*lag(Z,p)
#'         if (t==2) {
#'           ZpreviousTimesPhi = ZprevAll[1:(t-1),]*phit
#'         } else{
#'           ZpreviousTimesPhi = colSums(ZprevAll[1:(t-1),]*phit)
#'         }
#'
#'         # Recompute integral limits
#'         a = (qnorm(mycdf(xt[t]-1,theta1),0,1) - ZpreviousTimesPhi)/rt
#'         b = (qnorm(mycdf(xt[t],theta1),0,1) - ZpreviousTimesPhi)/rt
#'
#'         # compute random errors from truncated normal
#'         err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#'
#'         # compute the new Z and add it to the previous ones
#'         znew = rbind(ZpreviousTimesPhi + rt*err, ZprevAll[1:(t-1),])
#'         ZprevAll[1:t,] = znew
#'
#'         # recompute weights
#'         wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
#'         wprev = wgh[t,]
#'
#'         # use YW equation to compute estimates of phi and of the erros
#'         Gt = toeplitz(TacvfAR(AR)[1:t])
#'         gt = TacvfAR(AR)[2:(t+1)]
#'         phit = as.numeric(solve(Gt) %*% gt)
#'         rt =  as.numeric(sqrt(1 - gt %*% solve(Gt) %*% gt/TacvfAR(AR)[1]))
#'
#'       }
#'     }
#'
#'     # From p to T1 I dont need to estimate phi anymore
#'     for (t in (p+1):T1){
#'        # compute phi_1*Z_{t-1} + phi_2*Z_{t-2} for all particles
#'        if(p>1){# colsums doesnt work for 1-dimensional matrix
#'          ZpreviousTimesPhi = colSums(ZprevAll*AR)
#'         }else{
#'           ZpreviousTimesPhi=ZprevAll*AR
#'         }
#'
#'         # compute limits of truncated normal distribution
#'         a = as.numeric(qnorm(mycdf(xt[t]-1,theta1),0,1) - ZpreviousTimesPhi)/rt
#'         b = as.numeric(qnorm(mycdf(xt[t],theta1),0,1) -   ZpreviousTimesPhi)/rt
#'
#'         # draw errors from truncated normal
#'         err = qnorm(runif(length(a),0,1)*(pnorm(b,0,1)-pnorm(a,0,1))+pnorm(a,0,1),0,1)
#'
#'         # Update the underlying Gaussian series (see step 3 in SIS section in the paper)
#'         znew = ZpreviousTimesPhi + rt*err
#'         if (p>1){
#'           ZprevAll = rbind(znew, ZprevAll[1:(p-1),])
#'         }else {
#'           ZprevAll[1,]=znew
#'         }
#'
#'         # update weights
#'         wgh[t,] = wprev*(pnorm(b,0,1) - pnorm(a,0,1))
#'         wprev = wgh[t,]
#'     }
#'
#'     # likelihood approximation
#'     lik = mypdf(xt[1],theta1)*mean(na.omit(wgh[T1,]))
#'
#'     # for log-likelihood we use a bias correction--see par2.3 in Durbin Koopman, 1997
#'     nloglik = -log(lik)
#'     #- (1/(2*N))*(var(na.omit(wgh[T1,]))/mean(na.omit(wgh[T1,])))/mean(na.omit(wgh[T1,]))
#'
#'     out = (if (is.na(nloglik) | lik==0) 10^8 else nloglik)
#'
#'   }else{
#'     out = 10^8 # for noncasusal AR
#'   }
#'
#'   return(out)
#' }
