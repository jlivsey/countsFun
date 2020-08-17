#======================================================================================================#
#Purpose:   Function to evalutate Sandwich Std Errors for Count Models for NegBinom-AR(1) model
#
# Author:   Stefanos Kechagias and James Livsey
# Team:     Vladas Pipiras, James Livsey, Stefanos Kechagias, Robert Lund, Yisu Jia
# Date:     July 2020
#=====================================================================================================#

# library(countsFun)
# library(numDeriv)
# library(sandwich)
#
#
# # standard erros with sandwhich method following notation of Freedman (2006)
# sand <- function(theta, data, Regressor, mod){
#
#   # # vvvv Temporary inputs to test function vvvv
#   # ARMAorder <- c(1, 0)
#   # theta <- c(5, 1/4, 1/2)
#   #
#   # r <- theta[1]
#   # p <- theta[2]
#   # phi <- theta[3]
#   #
#   # #data <- sim_negbin(n = 200, ARMAmodel = list(phi, 0), r = r, p = p)
#   #
#   # n <- length(data)
#   #
#   # theta.est <- c(jitter(theta[1]),
#   #                jitter(theta[2]),
#   #                jitter(theta[3]))
#   #
#   # # ^^^^ END of temporary inputs ^^^^
#   # GaussianLogLik = function(theta, data, Regressor, mod)
#
#   # Calulate numerical Hessian
#   h <- gHgen(fn         = GaussianLogLik,
#               par       = theta,
#               data      = data,           # additional arg for GaussLogLikNB
#               Regressor = Regressor,      # additional arg for GaussLogLikNB
#               mod       = mod)            # additional arg for GaussLogLikNB
#
#
#   if(!(h$hessOK && det(h$Hn)>10^(-8))){
#     SE.sand = rep(NA, mod$nparms)
#   }else{
#     gi <- matrix(NA, length(data), length(theta))
#
#     for(k in 1:(length(data)-1)){
#       print(k)
#       gi[k, ] <- grad(func = logf_i, x = theta, mod = mod, i = k)
#     }
#     gi <- gi[-length(data), ]
#
#     # Calculate Cov matrix
#     A <- h
#     B <- t(gi) %*% gi
#     V.sand  <- solve(-A) %*% B %*% solve(-A)
#     SE.sand <- sqrt(diag(V.sand))
#
#     # standard errors usual way using hessian
#     #SE.hess <- sqrt(diag(solve(h)))
#   }
#
#   return(SE.sand)
# }
#
# #--------------log of density at observation i, Freedman (2006)
# logf_i <- function(theta, mod, i){
#
#   # retrieve ARMA parameters
#   AR = NULL
#   if(mod$ARMAorder[1]>0) AR = theta[(mod$nMargParms+1):(mod$nMargParms + mod$ARMAorder[1])]
#
#   MA = NULL
#   if(mod$ARMAorder[2]>0) MA = theta[(mod$nMargParms+mod$ARMAorder[1]+1) : (mod$nMargParms + mod$ARMAorder[1] + mod$ARMAorder[2]) ]
#
#   # retrieve marginal parameters
#   MargParms        = theta[mod$MargParmIndices]
#
#   # retrieve regressor parameters
#   if(mod$nreg>0){
#     beta  = MargParms[1:(mod$nreg+1)]
#     m     = exp(Regressor%*%beta)
#   }
#
#   # retrieve GLM type  parameters
#   if(mod$CountDist == "Negative Binomial" && mod$nreg>0){
#     ConstMargParm  = 1/MargParms[mod$nreg+2]
#     DynamMargParm  = MargParms[mod$nreg+2]*m/(1+MargParms[mod$nreg+2]*m)
#   }
#
#   if(mod$CountDist == "Generalized Poisson" && mod$nreg>0){
#     ConstMargParm  = MargParms[mod$nreg+2]
#     DynamMargParm  = m
#   }
#
#   if(mod$CountDist == "Poisson" && mod$nreg>0){
#     ConstMargParm  = NULL
#     DynamMargParm  = m
#   }
#
#
#   # retrieve mean
#   if(mod$nreg>0){
#     MeanValue = m
#   }else{
#     MeanValue = switch(mod$CountDist,
#                        "Poisson"             = MargParms[1],
#                        "Negative Binomial"   = MargParms[1]*MargParms[2]/(1-MargParms[2]),
#                        "Generalized Poisson" = MargParms[2])
#   }
#
#   # Compute truncation of relation (21)
#   if(mod$nreg>0){
#     N = sapply(unique(DynamMargParm),function(x)which(mod$mycdf(1:mod$MaxCdf, ConstMargParm, x)>=1-1e-7)[1])-1
#     N[is.na(N)] = mod$MaxCdf
#   }else{
#     N <- which(round(mod$mycdf(1:mod$MaxCdf, MargParms), 7) == 1)[1]
#     if(length(N)==0 |is.na(N) ){
#       N =mod$MaxCdf
#     }
#   }
#
#   # compute hermitte coef
#   HC = HermCoef(MargParms, ConstMargParm, DynamMargParm, N, mod$nHC, mod$mycdf, mod$nreg)
#
#
#
#   # ARMA autocorrelation function
#   if(!length(AR) && length(MA)){arma.acf <- ARMAacf(ma = MA, lag.max = n)}
#   if(!length(MA) && length(AR)){arma.acf <- ARMAacf(ar = AR, lag.max = n)}
#   if(length(AR)  && length(MA)){arma.acf <- ARMAacf(ar = AR, ma = MA, lag.max = n)}
#
#
#
#   # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
#   if(is.null(AR) && is.null(MA)){
#     gamma = c(as.numeric(HC^2%*%factorial(1:nHC)), rep(0,n-1))
#   }
#   else{
#     gamma = CountACVF(h = 0:(n-1), myacf = arma.acf, g = HC)
#   }
#
#   # DL algorithm
#   if(mod$nreg==0){
#     DLout <- DLalg(data, gamma, MeanValue)
#     ei <- DLout$e[i]
#     vi <- DLout$v[i]
#   }
#
#   # else{
#   #   INAlgout = innovations.algorithm(gamma)
#   #   Theta    = INAlgout$thetas
#   #   vi <- INAlgout$vs[i]
#   #
#   #   ei <- sum(data - Theta*data)
#   #
#   #
#   #
#   #   Theta = ia$thetas
#   #   # first stage of Innovations
#   #   v0    = ia$vs[1]
#   #   zhat = -Theta[[1]][1]*zprev
#   #
#   # }
#
#
#   return(-(log(vi) + ei^2/vi)/2)
# }
#
#



