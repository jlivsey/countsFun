# Compare CountsFun PF versus LCG Particle Filter
library(countsFun)
library(latentGaussCounts)
library(FitAR)

# setup parameters for Poisson(lam)-AR(p) series
lam = 3
phi = 0.5
n = 100             # sample size

# generate data
x=sim_pois_ar(n, phi, lam )

# initial parameters
initial.param = c(lam, phi)
theta = initial.param
data = x
CountDist = "Poisson"
ParticleNumber = 1000
ARMAorder = c(0,1)
CountDist = "Poisson"



# to match with old code I need to multiply with 2 and remove the bias correction
# LikSISGenDist_ARp_Res(initial.param,x,ParticleNumber,CountDist)
ParticleFilterMA1_Res(initial.param, x, ARMAorder, ParticleNumber, CountDist, 0)

likSISRMA1(initial.param,x,ARMAorder,CountDist)/2
