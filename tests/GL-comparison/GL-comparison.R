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

#-------------------------------------------------------------------------------
#' # Gaussian Lik'd as in paper

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

  # Print 1st 10 iterations numbers and then every 20th
  if(r < 10) print(r)
  if(r %% 20 == 0) print(r)

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

#-------------------------------------------------------------------------------
#' # 2-stage estimation Fixed Lambda

# ---- Variables needed that are functions of manual inputs ----
p = length(phi)  # AR order
nparms = p+1     # total number of parameters
HC <- HermCoef(lam)

# ---- allocate memory to save the following: ----
ParmEst2 = matrix(0, nrow=nsim, ncol=nparms)

# for-loop for each realization
for (r in 1:nsim){

  # Print 1st 10 iterations numbers and then every 20th
  if(r < 10) print(r)
  if(r %% 20 == 0) print(r)

  # generate Poisson-AR(p) data
  x <- sim_pois_ar(n, phi, lam )

  # Estimate lambda
  lam.est <- mean(x)

  # select initial parameters as true ones
  initial.param <- phi

  # run optimization for our model
  optim.output <-
    FitGaussianLik_fixedLambda(initialParam = initial.param,
                               x = x,
                               lam = lam.est,
                               HC = HC)

  # Store parameter estimates
  ParmEst2[r, ] <- c(lam.est, optim.output[1])

  # reset the seed
  set.seed(r)
}

#-------------------------------------------------------------------------------
#' # Approximate GL as in Vladas note


#' ---- Variables needed that are functions of manual inputs ----
p <- length(phi)  # AR order
nparms <- p+1     # total number of parameters

#' ---- allocate memory to save output: ----
ParmEst3 = matrix(0, nrow=nsim, ncol=nparms)

#' for-loop for each realization
for (r in 1:nsim){

  # Print 1st 10 iterations numbers and then every 20th
  if(r < 10) print(r)
  if(r %% 20 == 0) print(r)

  # generate Poisson-AR(p) data
  x <- sim_pois_ar(n, phi, lam )

  # select initial parameters as true ones
  initial.param <- c(lam, phi)

  # run optimization for our model
  optim.output <- FitGaussianLikApprox(initialParam = initial.param, x = x)

  # Store parameter estimates
  ParmEst3[r, ] <- optim.output[1:2]

  # reset the seed
  set.seed(r)
}


#' -----------------------------------------------------------------------------
#' # Plot Results

# Plot phi from all 3 methods
allPhi <- cbind(ParmEst1[,2], ParmEst2[,2], ParmEst3[,2])
boxplot(allPhi, main = "phi estimates")
abline(h = phi)

# Plot lambda from methods 1 and 2

allLam <- cbind(ParmEst1[, 1], ParmEst2[, 1], ParmEst3[, 1])
boxplot(allLam, main = "lambda estimates")
abline(h = lam)
