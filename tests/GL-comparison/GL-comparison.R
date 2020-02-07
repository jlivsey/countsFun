#' ---
#' title: "GL-comparison"
#' author: "James_Livsey"
#' output: pdf_document
#' ---


#' This is a comparison of three methods.
#'
#'    1. Gaussian Liklihood as done in paper
#'    2. Gaussian likelihood done as 2-stage proceedure.
#'       First estimate $\lambda$ with $\bar{X}$ and then only use GL to estimate $\phi$
#'    3. Approximate GL as in Vladas note W-12-23-19-pseudo likelihood-2.pdf


library(countsFun)

#' ---- setup parameters for Poisson(lam)-AR(p) series ----
lam  <- 2
phi  <- .75
n    <- 100     # sample size
nsim <- 50  # number of realizations

#' ---- Variables needed that are functions of manual inputs ----
p <- length(phi)  # AR order
nparms <- p+1     # total number of parameters

#' # 1. Gaussian Lik'd as in the paper

#' ---- allocate memory to save output: ----
ParmEst1 = matrix(0, nrow=nsim, ncol=nparms)

#' for-loop for each realization
for (r in 1:nsim){

  # generate Poisson-AR(p) data
  x <- sim_pois_ar(n, phi, lam )

  # select initial parameters as true ones
  initial.param <- c(lam, phi)

  # run optimization for our model
  optim.output <- FitGaussianLik(initialParam = initial.param, x = x)

  # Store parameter estimates
  ParmEst1[r, ] <- optim.output[1:2]

  # reset the seed
  set.seed(r)
}





