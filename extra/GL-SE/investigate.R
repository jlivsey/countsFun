# Purpose: look at 1-dimensional slices of the GL for mixed-Poisson models
# ------------------------------------------------------------------------------

# parameters
p         <- 1/4
lam1      <- 2
lam2      <- 10
phi       <- .5
n         <- 200
ARMAorder <- c(1, 0)
MaxCdf    <- 5000
nHC       <- 20
theta <- c(p, lam1, lam2, phi)
ARMAmodel <- list(phi, NULL)

# simulate data
data <- sim_mixedPoisson(n, ARMAmodel, p, lam1, lam2)


# ---- Plot as a function of p ------------------------------------------------
f.p <- rep(NA, 1000)
p.seq <- seq(.01, .49, length.out = 1000)
for(i in 1:1000){
  theta <- c(p.seq[i], lam1, lam2, phi)
  f.p[i] <- GaussLogLikMP(theta, data, ARMAorder = c(1, 0), MaxCdf = 5000, nHC = 20)
}
plot(f.p ~ p.seq)

# ---- Plot as a function of lam1 ------------------------------------------------
f.lam1 <- rep(NA, 1000)
lam1.seq <- seq(.01, 10, length.out = 1000)
for(i in 1:1000){
  if(i < 10 | i%%100==0) print(i)
  theta <- c(p, lam1.seq[i], lam2, phi)
  f.lam1[i] <- GaussLogLikMP(theta, data, ARMAorder = c(1, 0), MaxCdf = 5000, nHC = 20)
}
plot(f.lam1 ~ lam1.seq)
lam1.seq[which.min(f.lam1)]


# ---- Plot as a function of lam2 ------------------------------------------------
f.lam2 <- rep(NA, 1000)
lam2.seq <- seq(2, 20, length.out = 1000)
for(i in 1:1000){
  if(i < 10 | i%%100==0) print(i)
  theta <- c(p, lam1, lam2.seq[i], phi)
  f.lam2[i] <- GaussLogLikMP(theta, data, ARMAorder = c(1, 0), MaxCdf = 5000, nHC = 20)
}
plot(f.lam2 ~ lam2.seq)
lam2.seq[which.min(f.lam2)]

# ---- Plot as a function of phi ------------------------------------------------
f.phi <- rep(NA, 1000)
phi.seq <- seq(.05, .95, length.out = 1000)
for(i in 1:1000){
  if(i < 10 | i%%100==0) print(i)
  theta <- c(p, lam1, lam2, phi.seq[i])
  f.phi[i] <- GaussLogLikMP(theta, data, ARMAorder = c(1, 0), MaxCdf = 5000, nHC = 20)
}
plot(f.phi ~ phi.seq)
phi.seq[which.min(f.phi)]


