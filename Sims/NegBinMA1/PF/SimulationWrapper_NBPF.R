# PURPOSE: Call the NegBinAR1_GL function that performs the NegBin-MA(1) GL simulations for our paper.
#          Here we call it multiple times for different simulation schemes

# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# set directory to  save results
setwd("C:/Users/Stef/Desktop/countsFun/Sims/NegBinMA1/PF/RData")

# ---- Load libraries ----
library(parallel)
library(doParallel)
library(countsFun)

# fixed parameters across all simulation schemes
CountDist       = "Negative Binomial"
r               = 3
prob            = 0.8
MargParm        = c(r,prob)
nsim            = 200
no_cores <- detectCores()/2-1
LB = c(0.001, 0.001, -0.999)
UB = c(Inf, 0.999, 0.999)
p = 0
q = 1
ARMAorder = c(p,q)
Particles = 1000
epsilon = 1
#-----------------------------------------------Positive MA parameter--------------------------------------------------#
MAParm = 0.75
ThetaSign = ifelse(MAParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param
trueParam = c(4, 0.5, 0.5)
initParam = trueParam

n = 100
df1 = NegBin_PF(initParam, trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0)
save(df1, file = sprintf("NegBin%s_%sMA%s_PF_N%s_NS%s_Part%s_Theta%s_p.RData", MargParm[1], MargParm[2], 1, n, nsim,Particles, ThetaSign, p*10))

n = 200
df2 = NegBin_PF(initParam, trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0)
save(df2, file = sprintf("NegBin%s_%sMA%s_PF_N%s_NS%s_Part%s_Theta%s_p.RData", MargParm[1], MargParm[2], 1, n, nsim,Particles, ThetaSign, p*10))

n = 400
df3 = NegBin_PF(initParam, trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0)
save(df3, file = sprintf("NegBin%s_%sMA%s_PF_N%s_NS%s_Part%s_Theta%s_p.RData", MargParm[1], MargParm[2], 1, n, nsim,Particles, ThetaSign, p*10))



#-----------------------------------------------Negative MA parameter--------------------------------------------------#
MAParm = -0.75
ThetaSign = ifelse(MAParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param
trueParam = c(4, 0.5, -0.5)
initParam = trueParam

n = 100
df4 = NegBin_PF(initParam, trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0)
save(df4, file = sprintf("NegBin%s_%sMA%s_PF_N%s_NS%s_Part%s_Theta%s_p.RData", MargParm[1], MargParm[2], 1, n, nsim,Particles, ThetaSign, p*10))

n=200
df5 = NegBin_PF(initParam, trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0)
save(df5, file = sprintf("NegBin%s_%sMA%s_PF_N%s_NS%s_Part%s_Theta%s_p.RData", MargParm[1], MargParm[2], 1, n, nsim,Particles, ThetaSign, p*10))

n = 400
df6 = NegBin_PF(initParam, trueParam, p, q, LB, UB, n, nsim, Particles, epsilon, no_cores, 0)
save(df6, file = sprintf("NegBin%s_%sMA%s_PF_N%s_NS%s_Part%s_Theta%s_p.RData", MargParm[1], MargParm[2], 1, n, nsim,Particles, ThetaSign, p*10))







#
# ARMAmodel = list(NULL,NULL)
# if(p>0){ARMAmodel[[1]] = trueParam[3:(2+p)]}
# if(q>0){ARMAmodel[[2]] = trueParam[(3+p):length(trueParam)]}
#
# # Generate all the data and save it in a list
# l <- list()
# for(r in 1:nsim){
#   set.seed(r)
#   l[[r]] = sim_negbin(n, ARMAmodel, MargParm[1], MargParm[2] )
# }
#
# initParam = c(3.1088904 ,0.3869446, -0.1541726)
#
#
#
#
# optim.output<- DEoptim::DEoptim(par            = initParam,
#                                 fn             = ParticleFilterMA1_Res,
#                                 lower          = LB,
#                                 upper          = UB,
#                                 data           = X,
#                                 ARMAorder      = ARMAorder,
#                                 ParticleNumber = ParticleNumber,
#                                 CountDist      = CountDist,
#                                 epsilon        = epsilon,
#                                 control        = DEoptim::DEoptim.control(trace = 10, itermax = 200, steptol = 50, reltol = 1e-5))
#
#
#
#
#
#
#
# k1 = FitMultiplePFNew(initParam, l[[1]], CountDist, Particles, LB, UB, ARMAorder, epsilon, 1)
# k2 = FitMultiplePFNew(initParam, l[[1]], CountDist, Particles, LB, UB, ARMAorder, epsilon, 0)
#
#
# # ParticleFilterMA1_Res(initParam, l[[6]], C(0,1), Particles, CountDist, 0.5)
# # likSISRMA1(initParam, l[[6]], ARMAorder, CountDist)/2
# #




