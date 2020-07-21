

# ---- Load libraries ----
library(countsFun)
library(MixtureInf)


# fixed parameters across all simulation schemes
CountDist       = "NegBin"
r               = 3
prob            = 0.8
MargParm        = c(r,prob)
nsim            = 2
LB = c(0.001, 0.001, -0.999)
UB = c(Inf, 0.999, 0.999)
nsim            = 2
p = 0
q = 1
ARParm = 0.75
trueParam = c(MargParm, ARParm)
MaxCdf = 5000
nHC = 30

n = 100

# list with true ARMA parameters
ARMAmodel = list(NULL,NULL)
if(p>0){ARMAmodel[[1]] = trueParam[3:(2+p)]}
if(q>0){ARMAmodel[[2]] = trueParam[(3+p):length(trueParam)]}

# Generate all the data and save it in a list
x = sim_negbin(n, ARMAmodel, MargParm[1], MargParm[2] )
initParam = ComputeInitNegBinMA(x,n,nHC,LB,UB)
initParam

#====================================================================================================================#






