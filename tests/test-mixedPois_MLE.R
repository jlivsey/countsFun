# Purpose - small simulation to test mixedPois_MLE proceedue to get estimates
#           of the marginal parameters during IYW (mixedPoisson) estimation


# True parameters
n <- 300
phi <- 0
prob <- 1/4
lam1 <- 2
lam2 <- 10

nsim <- 100

# Set initial conditions
lam1.init <- jitter(lam1, 0)
lam2.init <- jitter(lam2, 0)
prob.init <- jitter(prob, 0)
iv <- c(log(lam1.init), log(lam2.init), log(prob.init))


# storage
out <- matrix(ncol = 3, nrow = nsim)


# main simulation loop
for(i in 1:nsim){
  x <- sim_mixpois_ar1(Tt = n, phi = phi, p = prob, lam1 = lam1, lam2 = lam2)
  est <- mixedPois_MLE(x = x, inital.value = c(log(lam1), log(lam2), prob))
  out[i, ] <- est
}

# look at results
boxplot(out[, 1]); abline(h = lam1); summary(out[, 1])
boxplot(out[, 2]); abline(h = lam2); summary(out[, 2])
boxplot(out[, 3]); abline(h = prob); summary(out[, 3])

