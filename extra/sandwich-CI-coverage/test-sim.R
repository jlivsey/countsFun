library(ggplot2)

# Source simulation function
source('~/github/countsFun/extra/sandwich-CI-coverage/sim-function.R')

# Simulation parameters
lam_tru <- 1
phi_tru <- .5
nsim    <- 200

# Run simulation
A <- sim_sandy(lam = lam_tru, phi = phi_tru, n = 100, nsim = 200)

# Extract values from simulation
lamhat <- A$P[, 1]
phihat <- A$P[, 2]

lamSeHess <- A$H[, 1]
lamSeSand <- A$S[, 1]

phiSeHess <- A$H[, 2]
phiSeSand <- A$S[, 2]

# ---- Lambda with Hessian ----
# data
hessCI <- data.frame(lower = lamhat - 1.96 * lamSeHess,
                     upper = lamhat + 1.96 * lamSeHess,
                     est   = lamhat,
                     simnum = 1:nsim)
# plot
g1 <- ggplot(hessCI, aes(x = simnum, ymin = lower, ymax = upper, y = est)) +
  geom_errorbar(color = 'blue') +
  geom_hline(yintercept = lam_tru, color = 'red') +
  ggtitle("Lambda with Hessian")

# number missed
(numMiss_lamHess <- sum(hessCI$lower > lam_tru | hessCI$upper < lam_tru))

# ---- Lambda with Sandwich ----
# data
sandCI <- data.frame(lower = lamhat - 1.96 * lamSeSand,
                     upper = lamhat + 1.96 * lamSeSand,
                     est   = lamhat,
                     simnum = 1:nsim)
# plot
g2 <- ggplot(sandCI, aes(x = simnum, ymin = lower, ymax = upper, y = est)) +
        geom_errorbar(color = 'blue') +
        geom_hline(yintercept = lam_tru, color = 'red') +
        ggtitle('Lambda with Sandwich')

# number missed
(numMiss_lamSand <- sum(sandCI$lower > lam_tru | sandCI$upper < lam_tru))


# ---- Phi with Hessian ----
# data
hessCI <- data.frame(lower = phihat - 1.96 * phiSeHess,
                     upper = phihat + 1.96 * phiSeHess,
                     est   = phihat,
                     simnum = 1:nsim)
# plot
g3 <- ggplot(hessCI, aes(x = simnum, ymin = lower, ymax = upper, y = est)) +
  geom_errorbar(color = 'blue') +
  geom_hline(yintercept = phi_tru, color = 'red') +
  ggtitle('Phi with Hessian')


# number missed
(numMiss_phiHess <- sum(hessCI$lower > phi_tru | hessCI$upper < phi_tru))

# ---- Phi with Sandwich ----
# data
sandCI <- data.frame(lower = phihat - 1.96 * phiSeSand,
                     upper = phihat + 1.96 * phiSeSand,
                     est   = phihat,
                     simnum = 1:nsim)
# plot
g4 <- ggplot(sandCI, aes(x = simnum, ymin = lower, ymax = upper, y = est)) +
  geom_errorbar(color = 'blue') +
  geom_hline(yintercept = phi_tru, color = 'red') +
  ggtitle('Phi with Sandwich')

# number missed
(numMiss_phiSand <- sum(sandCI$lower > phi_tru | sandCI$upper < phi_tru))

# Plot Results
gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2, nrow = 2)

# numerical results
numMiss_lamHess
numMiss_lamSand
numMiss_phiHess
numMiss_phiSand
