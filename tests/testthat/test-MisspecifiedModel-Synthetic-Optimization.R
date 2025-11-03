
# test_that("The likelihood function ParticleFilter_Res_ARMA works for a 'misspecified'
#           model on synthetic data.", {
#
#
# # Generate Mixed Poisson model
# set.seed(2)
# n              = 50
# CountDist      = "Mixed Poisson"
# MargParm       = c(2, 5, 0.7)
# ARParm         = 0.75
# MAParm         = NULL
# DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm)
#
# # generate a zero inflated (p =0.9) poisson 50 regressor
# Regressor      = ifelse(rbinom(n, size = 1, prob = 0.9) > 0, 0, rpois(n, lambda = 50))
# DependentVar    = DependentVar + Regressor
# #plot(1:n,DependentVar, type="l")
#
# # Set parameters for model to fit
# CountDist      = "Poisson"
# MargParm       = 55
# ARParm         = 0.75
# MAParm         = -0.3
# ARMAModel      = c(length(ARParm),length(MAParm))
# Task           = 'Optimization'
# OptMethod      = 'L-BFGS-B'
# initialParam   = c(MargParm,ARParm,MAParm)
# theta          = initialParam
# verbose        = FALSE
# # call the wrapper function with less arguments
# mod = ModelScheme(DependentVar   = DependentVar,
#                   Regressor      = NULL,
#                   CountDist      = CountDist,
#                   ARMAModel      = ARMAModel,
#                   OptMethod      = OptMethod,
#                   initialParam   = initialParam,
#                   verbose        = verbose)
#
# # Run the Optim function for now, later I ll do the test with the lgc wrapper
# optim.output <- optimx(
#   par     = theta,
#   fn      = ParticleFilter,
#   lower   = mod$LB,
#   upper   = mod$UB,
#   hessian = TRUE,
#   method  = mod$OptMethod,
#   mod     = mod)
#
# expect_equal(optim.output[[1]], 54.99979, tolerance = 10^(-4))
# expect_equal(optim.output[[2]], 0.7803059, tolerance = 10^(-4))
# expect_equal(optim.output[[3]], -0.2923531, tolerance = 10^(-4))
# expect_equal(optim.output[[4]], 654.9023, tolerance = 10^(-3))
# })



test_that("The likelihood function ParticleFilter works for a 'misspecified'
          model on synthetic data - Using the new misspec version of the loglik function.", {


            # Generate Mixed Poisson model
            set.seed(2)
            n              = 50
            CountDist      = "Mixed Poisson"
            MargParm       = c(2, 5, 0.7)
            ARParm         = 0.75
            MAParm         = NULL
            DependentVar   = sim_lgc(n, CountDist, MargParm, ARParm, MAParm)

            # generate a zero inflated (p =0.9) poisson 50 regressor
            Regressor      = ifelse(rbinom(n, size = 1, prob = 0.9) > 0, 0, rpois(n, lambda = 50))
            DependentVar    = DependentVar + Regressor
            #plot(1:n,DependentVar, type="l")

            # Set parameters for model to fit
            CountDist      = "Poisson"
            MargParm       = 55
            ARParm         = 0.75
            MAParm         = -0.3
            ARMAModel      = c(length(ARParm),length(MAParm))
            Task           = 'Optimization'
            OptMethod      = 'L-BFGS-B'
            initialParam   = c(MargParm,ARParm,MAParm)
            theta          = initialParam
            verbose        = FALSE
            # call the wrapper function with less arguments
            mod = ModelScheme(DependentVar   = DependentVar,
                              Regressor      = NULL,
                              CountDist      = CountDist,
                              ARMAModel      = ARMAModel,
                              OptMethod      = OptMethod,
                              initialParam   = initialParam,
                              verbose        = verbose)

            # Run the Optim function for now, later I ll do the test with the lgc wrapper
            optim.output <- optimx(
              par     = theta,
              fn      = ParticleFilter,
              lower   = mod$LB,
              upper   = mod$UB,
              hessian = TRUE,
              method  = mod$OptMethod,
              mod     = mod)

            expect_equal(optim.output[[1]], 54.99974, tolerance = 10^(-4))
            expect_equal(optim.output[[2]], 0.8010119, tolerance = 10^(-4))
            expect_equal(optim.output[[3]], -0.290344, tolerance = 10^(-4))
            expect_equal(optim.output[[4]], 631.9192, tolerance = 10^(-4))
          })

