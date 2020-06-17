# Compare CountsFun PF versus LCG Particle Filter
library(countsFun)
#library(latentGaussCounts)
library(FitAR)

# setup parameters for Poisson(lam)-AR(p) series
lam1 = 2
lam2 = 5
prob = 0.25
phi = 0.75
n = 100             # sample size
ARMAorder = c(1,0)
CountDist = "Mixed Poisson"
ParticleNumber = 100
p = 1
q = 0
nHC         = 30
LB          = c(0.001, 0.01, 0.01, -0.995)
UB          = c(0.499, Inf, Inf, 0.995)
# list with true ARMA parameters
ARMAmodel = list(NULL,NULL)
trueParam = c(prob, lam1,lam2,phi)
if(p>0){ARMAmodel[[1]] = trueParam[4:(3+p)]}
if(q>0){ARMAmodel[[2]] = trueParam[(4+p):length(trueParam)]}

# generate data
set.seed(3)
x=sim_mixedPoisson(n, ARMAmodel,prob,lam1,lam2 )

# initial parameters
initial.param = trueParam


# to match with old code I need to multiply with 2 and remove the bias correction
# LikSISGenDist_ARp_Res(initial.param,x,ParticleNumber,CountDist)
#ParticleFilterRes(initial.param, x, ARMAorder, ParticleNumber, CountDist, 1)

#likSISRAR1(initial.param,x,ARMAorder,CountDist)/2
data = x


theta = ComputeInitMixedPoissonAR(data,n,nHC,LB,UB)


optim.output <- optimx(par            = theta,
                      fn             = likSISRAR1,
                      data           = data,
                      ARMAorder      = ARMAorder,
                      CountDist      = CountDist,
                      lower          = LB,
                      upper          = UB,
                      hessian        = TRUE,
                      method         = "L-BFGS-B")
optim.output


#x0 = c(0.001000,1.864859,5.255515,0.995000)
#ParticleFilterRes(x0, x, ARMAorder, ParticleNumber, CountDist, 1)









