# #----------log of density at observation i, Freedman (2006)
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
#
#   # Autocovariance of count series--relation (9) in https://arxiv.org/pdf/1811.00203.pdf
#   GAMMA =  CountCovariance(mod$n, MargParms, ConstMargParm, DynamMargParm, AR, MA, N, mod$nHC, mod$mycdf, mod$nreg, 0)
#
#
#   # DL algorithm
#   if(mod$nreg==0){
#     DLout <- DLalg(data, GAMMA, MeanValue)
#     ei <- DLout$e[i]
#     vi <- DLout$v[i]
#   }else{
#     IAout <- Innalg(data, GAMMA)
#     ei <- IAout$e[i]
#     vi <- IAout$v[i]
#   }
#
#
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
