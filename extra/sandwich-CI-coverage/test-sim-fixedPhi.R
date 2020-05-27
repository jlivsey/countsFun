library(ggplot2)

# Source simulation function
source('~/github/countsFun/extra/sandwich-CI-coverage/function-sim-fixedPhi.R')

# Simulation parameters
lam_tru <- 1
phi_tru <- 0 # FIXED - This parameter is fixed in this script

# Run simulation
B <- sim_sandy_fixedPhi(lam = lam_tru, phi = phi_tru, n = 100, nsim = 200)

# Extract values from simulation
lamhat <- B$P[, 1]

lamSeHess <- B$H[, 1]
lamSeSand <- B$S[, 1]

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

# Plot Results
gridExtra::grid.arrange(g1, g2, ncol = 2)

# numerical results
numMiss_lamHess
numMiss_lamSand
