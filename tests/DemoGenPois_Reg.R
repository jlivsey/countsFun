# Compare Runtimes of the NegBinomial likelihood evaluation with and without covariates
# I need to vectorize the covariate version. It is too slow.

# libraries
library(countsFun)
library(tictoc)
library(lavaSearch2)
library(optimx)

symmetrize <- lavaSearch2:::symmetrize
# sample size
n = 200

# generate a regressor dummy variable
Regressor  = cbind(rep(1,n),rbinom(n, 1 ,0.4))

# specify parameters
b0 = 1
b1 = 2
b = rbind(b0,b1)
lambda = exp(Regressor%*%b)
omega = 0.5
phi = -0.5
MaxCdf = 100
nHC = 30

# list with true ARMA parameters
ARMAmodel = list(NULL,NULL)
ARMAmodel[[1]] = phi
ARMAorder = c(1,0)

# parameters to be passed in likelihood
theta = c(b,omega,phi)

# generate Gen Pois series
y = sim_genpois(n, ARMAmodel, lambda, omega)


# test likelihood
# tic()
# GaussLogLikGP_Reg(theta = theta,
#                 data = y,
#                 Regressor = Regressor,
#                 ARMAorder = ARMAorder,
#                 MaxCdf = MaxCdf,
#                 nHC = nHC)
# toc()


LB = c(-Inf, -Inf, 0.001, -0.999)
UB = c(Inf, Inf, 0.999, 0.999)

optim.output <- optimx(par      = theta,
                      fn        = GaussLogLikGP_Reg,
                      data      = y,
                      Regressor = Regressor,
                      ARMAorder = ARMAorder,
                      MaxCdf    = MaxCdf,
                      nHC       = nHC,
                      hessian   = TRUE,
                      lower     = LB,
                      upper     = UB,
                      control = list(all.methods = TRUE,trace = 2))

xt <- cbind(y,Regressor[,2])


out = LGC(xt,count.family = "GenPoisson",
          gauss.series = "AR",
          estim.method = "gaussianLik", p=1)


