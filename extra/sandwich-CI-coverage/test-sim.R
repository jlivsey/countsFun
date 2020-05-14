library(ggplot2)

# Run simulation
A <- sim_sandy(5, .5)

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
                     simnum = 1:100)
# plot
g1 <- ggplot(hessCI, aes(x = simnum, ymin = lower, ymax = upper, y = est)) +
  geom_errorbar(color = 'blue') +
  geom_hline(yintercept = 5, color = 'red') +
  ggtitle("Lambda with Hessian")

# number missed
sum(hessCI$lower > 5 | hessCI$upper < 5)

# ---- Lambda with Sandwich ----
# data
sandCI <- data.frame(lower = lamhat - 1.96 * lamSeSand,
                     upper = lamhat + 1.96 * lamSeSand,
                     est   = lamhat,
                     simnum = 1:100)
# plot
g2 <- ggplot(sandCI, aes(x = simnum, ymin = lower, ymax = upper, y = est)) +
        geom_errorbar(color = 'blue') +
        geom_hline(yintercept = 5, color = 'red') +
        ggtitle('Lambda with Sandwich')

# number missed
sum(sandCI$lower > 5 | sandCI$upper < 5)


# ---- Phi with Hessian ----
# data
hessCI <- data.frame(lower = phihat - 1.96 * phiSeHess,
                     upper = phihat + 1.96 * phiSeHess,
                     est   = phihat,
                     simnum = 1:100)
# plot
g3 <- ggplot(hessCI, aes(x = simnum, ymin = lower, ymax = upper, y = est)) +
  geom_errorbar(color = 'blue') +
  geom_hline(yintercept = .5, color = 'red') +
  ggtitle('Phi with Hessian')
# number missed
sum(hessCI$lower > .5 | hessCI$upper < .5)

# ---- Phi with Sandwich ----
# data
sandCI <- data.frame(lower = phihat - 1.96 * phiSeSand,
                     upper = phihat + 1.96 * phiSeSand,
                     est   = phihat,
                     simnum = 1:100)
# plot
g4 <- ggplot(sandCI, aes(x = simnum, ymin = lower, ymax = upper, y = est)) +
  geom_errorbar(color = 'blue') +
  geom_hline(yintercept = .5, color = 'red') +
  ggtitle('Phi with Sandwich')
# number missed
sum(sandCI$lower > .5 | sandCI$upper < .5)

# Plot Results
gridExtra::grid.arrange(g1, g2, g3, g4, nrow = 2, ncol = 2)
