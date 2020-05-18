# Simulation function
simulation_mixedPoissonAR1_GL <- function(n, lam1, lam2, phi, prob, nsim){

  # PURPOSE: Produce Mixed-Poisson AR(1) GL estimates
  #
  #
  # AUTHORS: Stefanos Kechagias, James Livsey
  #
  # DATE:    May 2020
  #
  # R version 3.6.3

  # ---- Load libraries ----
  library(countsFun)

  #Create columns lam.est, phi.est, estim.method, n, phi, phi.se, lam, lam.se
  dfCols    = c('estim.method', 'n',
                'phi.true',  'phi.est',  'phi.se',
                'lam1.true', 'lam1.est', 'lam1.se',
                'lam2.true', 'lam2.est', 'lam2.se',
                'prob.true', 'prob.est', 'prob.se')

  # Prepare results for the plot.
  df = data.frame(matrix(ncol = length(dfCols), nrow = nsim))
  names(df) <- dfCols

  for(r in 1:nsim){

    set.seed(r)
    x <- sim_MixPois_ar(n    = n,
                        phi  = phi,
                        lam1 = lam1,
                        lam2 = lam2,
                        prob = prob)

    fit <- fit_GL_mixedPois_AR1(x)

    # estimates from this replication to data.frame
    df[r,1] = 'GL'
    df[r,2] = n
    df[r,3]  = phi
    df[r,4] = fit[4]
    df[r,5] = NA
    df[r,6] = lam1
    df[r,7] = fit[1]
    df[r,8] = NA
    df[r,9] = lam2
    df[r,10] = fit[2]
    df[r,11] = NA
    df[r,12] = prob
    df[r,13] = fit[3]
    df[r,14] = NA

  }

  return(df)
}
