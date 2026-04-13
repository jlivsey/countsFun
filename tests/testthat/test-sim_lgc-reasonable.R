#================================================================================#
# PURPOSE: Check that sim_lgc returns reasonable moments
# Author: Stefanos Kechagias
#================================================================================#

test_that("sim_lgc matches several Negative Binomial marginal summaries", {

  n         <- 200
  nsim      <- 50
  ARParm    <- NULL
  MAParm    <- NULL
  CountDist <- "Negative Binomial"
  MargParm  <- c(5, 0.3)

  size <- MargParm[1]
  prob <- MargParm[2]

  RegModel <- DependentVar ~ 1

  SUMM <- matrix(NA_real_, nrow = nsim, ncol = 5)
  colnames(SUMM) <- c("mean", "var", "p0", "q50", "q90")

  for (i in seq_len(nsim)) {
    set.seed(i)
    DependentVar <- as.numeric(
      sim_lgc(n, CountDist, MargParm, ARParm, MAParm, RegModel)
    )

    SUMM[i, "mean"] <- mean(DependentVar)
    SUMM[i, "var"]  <- var(DependentVar)
    SUMM[i, "p0"]   <- mean(DependentVar == 0)
    SUMM[i, "q50"]  <- unname(quantile(DependentVar, 0.50, type = 1))
    SUMM[i, "q90"]  <- unname(quantile(DependentVar, 0.90, type = 1))
  }

  TRUE_VALS <- c(
    mean = size * (1 - prob) / prob,
    var  = size * (1 - prob) / prob^2,
    p0   = dnbinom(0, size = size, prob = prob),
    q50  = qnbinom(0.50, size = size, prob = prob),
    q90  = qnbinom(0.90, size = size, prob = prob)
  )

  AVG <- colMeans(SUMM)

  tol_mean <- 3 * sd(SUMM[, "mean"]) / sqrt(nsim)
  tol_var  <- 3 * sd(SUMM[, "var"])  / sqrt(nsim)
  tol_p0   <- 3 * sqrt(TRUE_VALS["p0"] * (1 - TRUE_VALS["p0"]) / (n * nsim))

  expect_equal(AVG["mean"], TRUE_VALS["mean"], tolerance = tol_mean)
  expect_equal(AVG["var"],  TRUE_VALS["var"],  tolerance = tol_var)
  expect_true(abs(AVG["p0"] - TRUE_VALS["p0"]) <= tol_p0)

  expect_true(abs(AVG["q50"] - TRUE_VALS["q50"]) <= 1)
  expect_true(abs(AVG["q90"] - TRUE_VALS["q90"]) <= 2)
})
