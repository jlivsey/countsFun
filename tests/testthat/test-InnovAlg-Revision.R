#----------------------------------------------------------------------------------------------------#
# Purpose: Test the reworked version of the InnovAlg function.
#
# Notes:   In September 2025 I reworked the InnovAlg function. the stopping condition wasnt implemented
#          very well for ARMA(p,q) with p>0 and q>0. Sin ce the stopping condition essentially controls
#          the approximation that takes place, as expected the new implementation  resulted in a slighlthy
#          different loglike value in the test-lgc-NegBin_ARMA_2_4-Evaluation.R. I am showing this diff in
#          the last test below.
#
# Author:  Stefanos Kechagias
# Date:    Sep 2025
#----------------------------------------------------------------------------------------------------#
test_that("InnovAlg and InnovAlgOld return the same theta sequence - ARMA(1,1)", {
  # Set seed and parameters
  set.seed(123)
  n     <- 1000
  phi   <- 0.6
  theta <- 0.4

  # Simulate ARMA(1,1) data
  arma_process <- arima.sim(n = n, list(ar = phi, ma = theta))

  # Compute autocovariances
  gamma <- acf(arma_process, plot = FALSE, lag.max = 20)$acf[, , 1]

  # Set parameter and model list
  Parms <- list(AR = phi, MA = theta)
  mod <- list(nAR = 1, nMA = 1, maxdiff = 1e-6)

  # Run both algorithms
  out_new <- InnovAlg(Parms, gamma, mod)
  out_old <- InnovAlgOld(Parms, gamma, mod)

  # Compare common parts of theta vectors
  n_common <- min(length(out_old$thetas), length(out_new$thetas))

  for (i in 1:n_common) {
    expect_equal(out_old$thetas[[i]], out_new$thetas[[i]], tolerance = 1e-6)
  }
})

test_that("InnovAlg and InnovAlgOld return the same theta sequence for ARMA(2,2)", {
  # Set seed and ARMA(2,2) parameters
  set.seed(123)
  n     <- 1000
  phi   <- c(0.5, -0.3)
  theta <- c(0.4, 0.2)

  # Simulate ARMA(2,2) data
  arma_process <- arima.sim(n = n, list(ar = phi, ma = theta))

  # Compute sample autocovariances
  gamma <- acf(arma_process, plot = FALSE, lag.max = 20)$acf[, , 1]

  # Define parameters and model structure
  Parms <- list(AR = phi, MA = theta)
  mod <- list(nAR = 2, nMA = 2, maxdiff = 1e-6)

  # Run both versions of the algorithm
  out_new <- InnovAlg(Parms, gamma, mod)
  out_old <- InnovAlgOld(Parms, gamma, mod)

  # Compare only the overlapping theta estimates
  n_common <- min(length(out_old$thetas), length(out_new$thetas))

  for (i in 1:n_common) {
    expect_equal(out_old$thetas[[i]], out_new$thetas[[i]], tolerance = 1e-6)
  }
})

test_that("InnovAlg and InnovAlgOld return the same theta sequence for ARMA(2,4)", {
  # Set seed and define ARMA(2,4) parameters
  set.seed(123)
  n     <- 1000
  phi   <- c(0.4, 0.4)
  theta <- c(0.3, 0.5, -0.2, 0.1)

  # Simulate ARMA(2,4) time series
  arma_process <- arima.sim(n = n, list(ar = phi, ma = theta))

  # Compute sample autocovariances
  gamma <- acf(arma_process, plot = FALSE, lag.max = 20)$acf[, , 1]

  # Define parameters and model structure
  Parms <- list(AR = phi, MA = theta)
  mod <- list(nAR = 2, nMA = 4, maxdiff = 1e-6)

  # Run both algorithm versions
  out_new <- InnovAlg(Parms, gamma, mod)
  out_old <- InnovAlgOld(Parms, gamma, mod)

  # Determine overlapping length
  n_common <- min(length(out_old$thetas), length(out_new$thetas))

  # Compare theta estimates with tolerance
  for (i in 1:n_common) {
    expect_equal(out_old$thetas[[i]], out_new$thetas[[i]], tolerance = 1e-12)
  }
})

test_that("New and old particle filter implementations give slightly un-equal results", {
  # Parameters
  n         <- 200
  CountDist <- "Negative Binomial"
  MargParm  <- c(3, 0.8)
  ARParm    <- c(0.4, 0.4)
  MAParm    <- c(0.3, 0.5, -0.2, 0.1)
  ARMAModel <- c(length(ARParm), length(MAParm))

  # specify the regression formula (no regressors here)
  RegModel       = DependentVar ~ 1

  # Simulate data
  set.seed(2)
  DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm, RegModel)
  df = data.frame(DependentVar)

  # populate a list with the model characteristics
  mod = ModelScheme(RegModel = RegModel,
                    df = df,
                    CountDist = CountDist,
                    ARMAModel = ARMAModel)

  # Get initial parameter estimates
  theta1 <- InitialEstimates(mod)

  # Run particle filters
  set.seed(1)
  a1 <- ParticleFilter(theta1, mod)
  a2 <- ParticleFilter_Res_ARMAOld(theta1, mod)

  # Compare results -
  expect_gt(abs(a1-a2),  0.00001)

  # Compare results
  expect_equal(a1, a2, tolerance = 1e-5)
})








