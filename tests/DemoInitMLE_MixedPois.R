

# ---- Load libraries ----
library(countsFun)
library(MixtureInf)


# fixed parameters across all simulation schemes
nsim            = 2
p = 1
q = 0
ARParm = 0.75
lam1            = 2
lam2            = 5
prob            = 0.25
MargParm        = c(prob, lam1, lam2)
trueParam = c(MargParm, ARParm)
LB = c(0.001, 0.01, 0.01, -0.995)
UB = c(0.499, Inf, Inf, 0.995)
MaxCdf = 5000
nHC = 30

n =400

# list with true ARMA parameters
ARMAmodel = list(NULL,NULL)
if(p>0){ARMAmodel[[1]] = trueParam[4:(3+p)]}
if(q>0){ARMAmodel[[2]] = trueParam[(4+p):length(trueParam)]}

# Generate all the data and save it in a list
l <- list()
initParam <- list()
for(r in 1:nsim){
  set.seed(r)
  l[[r]] = sim_mixedPoisson(n, ARMAmodel, MargParm[1], MargParm[2], MargParm[3] )
}
#====================================================================================================================#

X = l[[2]]
x0 = ComputeInitMixedPoissonAR(X,n,nHC,LB,UB)
x0




ARMAorder = c(1,0)


optim.output <- optim(par       = x0,
                      fn        = GaussLogLikMP,
                      data      = X,
                      ARMAorder = ARMAorder,
                      MaxCdf    = MaxCdf,
                      nHC       = nHC,
                      method    = "L-BFGS-B",
                      hessian   = TRUE,
                      lower     = LB,
                      upper     = UB
)

# save estimates, loglik and standard errors
ParmEst   = c(optim.output$p1,optim.output$p2,optim.output$p3,optim.output$p4)
loglik               = optim.output$value

# compute hessian
H = gHgen(par   = ParmEst,
          fn    = GaussLogLikMP,
          data  = X,
          ARMAorder = ARMAorder,
          MaxCdf    = MaxCdf,
          nHC       = nHC
)$Hn

sqrt(abs(diag(solve(H))))

GaussLogLikMP(trueParam, X, ARMAorder, MaxCdf, nHC)




