ARMAorder = c(1,0)
MaxCdf = 1000
nHC = 30
Regressor = cbind(rep(1,length(Buy)),Buy)
data = MOVE
b0 = 2
b1 = 1
p = 0.6
phi = -0.4
theta = c(b0,b1,p,phi)
CountDist = "Negative Binomial"
Particles = 1000
epsilon =0.5
LB = c(-100, -100, 0.001, -Inf)
UB = c( 100,  100, 0.999, Inf)
Model = 1
OptMethod = "bobyqa"
#GaussLogLikNB_Reg2(theta, data, Regressor, ARMAorder, MaxCdf, nHC, CountDist)


mod = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB,
                               ARMAorder, epsilon, MaxCdf, nHC, Model, OptMethod)


LB = c(-100, -100, 0.001, -Inf)
UB = c( 100,  100, Inf  , Inf)
mod0 = FitGaussianLikNB_Reg(theta, MOVE, Regressor, CountDist, Particles, LB, UB,
                           ARMAorder, epsilon, MaxCdf, nHC, 0, OptMethod)
