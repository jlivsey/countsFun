#' @title Shared Parameters for LGC Model Functions
#' @description This function exists only to hold shared parameter documentation.
#' @param mod A list with model specifications, including:
#'   \itemize{
#'     \item \code{DependentVar}: observed count series
#'     \item \code{CountDist}: count distribution (e.g., "Poisson", "Negative Binomial")
#'     \item \code{ARMAModel}: vector of length 2 specifying AR and MA orders
#'     \item \code{ParticleNumber}: number of particles to use
#'     \item \code{nreg}: number of regressors
#'     \item \code{mypdf}: function to compute marginal PDF
#'     \item \code{mycdf}: function to compute marginal CDF
#'     \item \code{n}: sample size
#'     \item \code{maxdiff}: convergence threshold for Innovations Algorithm
#'   }
#' @keywords internal
mod_shared_params <- function(mod) NULL
