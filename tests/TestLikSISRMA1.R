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
library(latentGaussCounts)

# fixed parameters across all simulation schemes
CountDist       = "Negative Binomial"
r               = 3
prob            = 0.8
MargParm        = c(r,prob)
nsim            = 5
no_cores <- detectCores()-4
LB = c(0.001, 0.001, -0.999)
UB = c(200, 0.999, 0.999)
p = 0
q = 1
ARMAorder = c(p,q)
Particles = 1000

#-----------------------------------------------Positive MA parameter--------------------------------------------------#
MAParm = -0.75
ThetaSign = ifelse(MAParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param
trueParam = c(MargParm, MAParm)
initParam = trueParam


ARMAmodel = list(NULL,NULL)
if(p>0){ARMAmodel[[1]] = trueParam[3:(2+p)]}
if(q>0){ARMAmodel[[2]] = trueParam[(3+p):length(trueParam)]}
l <- list()
n=400
for(r in 1:nsim){
  set.seed(r)
  l[[r]] = sim_negbin(n, ARMAmodel, MargParm[1], MargParm[2] )
}

data = l[[1]]
theta = initParam
ParticleNumber = 1000
epsilon = 1


theta1 = c(3, 2, 0.390)





# my function
ParticleFilterMA1(theta1, l[[2]], ARMAorder, Particles, CountDist)

# Vladas function
likSISMA1(theta1, l[[2]],ARMAorder, CountDist)/2




# my function
ParticleFilterMA1_Res(theta1, l[[2]], ARMAorder, Particles, CountDist)

# Vladas function
likSISRMA1(theta1, l[[2]],ARMAorder, CountDist)/2





















X = l[[2]]
optim.outputR1 <- optim(par      = initParam,
                      fn        = ParticleFilterMA1_Res,
                      data      = X,
                      ARMAorder = ARMAorder,
                      epsilon   = 1,
                      CountDist = CountDist,
                      ParticleNumber = Particles,
                      method    = "L-BFGS-B",
                      hessian   = TRUE,
                      lower     = LB,
                      upper     = UB
)



optim.outputR <- optim(par       = initParam,
                      fn        = likSISRMA1,
                      data      = l[[2]],
                      ARMAorder = ARMAorder,
                      CountDist = CountDist,
                      method    = "L-BFGS-B",
                      hessian   = TRUE,
                      lower     = LB,
                      upper     = UB
)

optim.output<- DEoptim::DEoptim(fn             = likSISRMA1,
                                lower          = LB,
                                upper          = UB,
                                data           = l[[2]],
                                ARMAorder      = ARMAorder,
                                CountDist      = CountDist,
                                control        = DEoptim::DEoptim.control(trace = 10, itermax = 200, steptol = 50, reltol = 1e-5))





