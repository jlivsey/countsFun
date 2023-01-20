# # Calculates link coeficients mapping gam_X --> gam_Z
# # Inputs: hermite coefficients and marginal variance
# link_coefs <- function(g_coefs, gamx0){
#   K <- length(g_coefs)
#   out <- factorial(1:K) * g_coefs^2 / gamx0
#   return(out)
# }
#
# # test comment
# reversion <- function(A){
#   # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   #   %
#   # % INPUT
#   # %
#   # % A : K:1 vector !!!!!! (not 1:K vector)
#   # %
#   # %
#   # % OUTPUT
#   # %
#   # % a : 1:K sequence obtained by reversion from A
#   # %
#   # % Note: this use calculations of reversion through matrix on p. 47 in
#   # % Henrici (1974, volume 1)
#   # %
#   # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   # function [a] = reversion(A) ;
#
#   # K = length(A) ;
#   K = length(A)
#   # A_matrix = A*ones(1,K) ;
#   A_matrix = A %*% matrix(1, 1, K)
#   # A_matrix = [A_matrix; zeros(1,K)] ;
#   A_matrix = rbind(A_matrix, rep(0, K))
#   # A_matrix = reshape(A_matrix,K,K+1) ;
#   dim(A_matrix) = c(K, K+1)
#   # A_matrix = tril(A_matrix) ;
#   A_matrix[upper.tri(A_matrix)] = 0
#   # A_matrix = A_matrix' ;
#   A_matrix = t(A_matrix)
#   # A_matrix = A_matrix(1:K,:) ;
#   A_matrix = A_matrix[1:K, ]
#
#   # a_matrix = A' ;
#   a_matrix = t(A)
#   # A_matrix_power = A_matrix ;
#   A_matrix_power = A_matrix
#   # for k=2:K
#   for(k in 2:K){
#     # A_matrix_power = A_matrix_power*A_matrix ;
#     A_matrix_power = A_matrix_power %*% A_matrix
#     # a_matrix = [a_matrix;
#     #            [zeros(1,k-1) A_matrix_power(1,1:K-k+1)]] ;
#     a_matrix = rbind(a_matrix, c(rep(0, k-1), A_matrix_power[1, 1:(K-k+1)]))
#     # end ;
#   }
#
#   # a_matrix = a_matrix^(-1) ;
#   a_matrix = solve(a_matrix)
#   # a = a_matrix(1,:) ;
#   a = a_matrix[1, ]
#
#   return(a)
# }
#
# # Power series f(x) = \sum_k coef_k * x^k
# temp_power_series <- function(x, coef){
#   K = length(coef)
#   return(sum(coef * x^{1:K}))
# }
# power_series = Vectorize(temp_power_series, vectorize.args = "x")
#
#
# # Implied Yule-Walker for Poisson-AR(1) model
# fit_IWY <- function(x, p){
#   # Estimate lambda values from data
#   lam.est <- mean(x)
#
#   # Calculate Hermite coefficients
#   g.coefs <- HermCoef(lam.est)
#
#   # Link Coefficients of el: gam.x --> gam.z (autocorrelations--relation 29 May 1 version)
#   link.coefs <- link_coefs(g.coefs, lam.est)
#
#   # Inverse Link coefficients of f^-1: gam.z --> gam.x
#   #   (Perform reversion)
#   inv.link.coefs <- reversion(link.coefs)
#
#   # sample acf of count data
#   gam.x <- acf(x = x, lag.max = p, plot = FALSE, type = "correlation")$acf
#
#   # compute gamma Z thru reversion
#   gam.z <- power_series(gam.x[,,1], inv.link.coefs)
#   gam.z[1] <- 1
#
#   # Construct block toeplitz matrix from given ACF function and AR(P)
#   R = toeplitz(gam.z[1:p])
#
#   # stack the gammas
#   r = gam.z[2:(p+1)]
#
#   # multivariate Yule-Walker
#   my.yw = r %*% solve(R)
#
#   return(c(lam.est, my.yw))
#
# }
#
#
# # Negative Binomial parameter converters
# rp2musig <- function(r, p){
#   mu <- p*r / (1-p)
#   sig2 <- p * r / (1-p)
#   return(c(mu, sqrt(sig2)))
# }
#
# # Negative Binomial parameter converters
# musig2rp <- function(mu, sig){
#   r <- mu^2 / (sig^2 - mu)
#   p <- (sig^2 - mu) / sig^2
#   return(c(r, p))
# }
#
# # Function Implied-Yule-Walker for mixed-Poisson-ARmodel
# fit_IWY_mixedPois_AR1 <- function(x, ar.order = 1){
#
#   n <- length(x)
#
#   prob.est <- 1/4
#   n1 <- floor(n * prob.est)
#   x.sort <- sort(x)
#   lam1.est <- mean(x.sort[1:n1])
#   lam2.est <- mean(x.sort[(n1+1):n])
#
#   # mleFit <- mixedPois_MLE(x, c(log(lam1.est), log(lam2.est), prob.est))
#   # lam1.est <- mleFit[1]
#   # lam2.est <- mleFit[2]
#   # prob.est <- mleFit[3]
#
#   # Calculate Hermite coefficients
#   g.coefs <- HermCoefMixedPois(lam1 = lam1.est,
#                                lam2 = lam2.est,
#                                prob = prob.est)
#
#   # Link Coefficients of el: gam.x --> gam.z (autocorrelations--relation 29 May 1 version)
#   link.coefs <- link_coefs(g.coefs,
#                            # gamx0 = prob.est * lam1.est + (1-prob.est)*lam2.est)
#                            gamx0 = var(x))
#
#   # Inverse Link coefficients of f^-1: gam.z --> gam.x
#   #   (Perform reversion)
#   inv.link.coefs <- reversion(link.coefs)
#
#   # sample acf of count data
#   gam.x <- acf(x = x, lag.max = ar.order, plot = FALSE, type = "correlation")$acf
#
#   # compute gamma Z thru reversion
#   gam.z <- power_series(gam.x[,,1], inv.link.coefs)
#   gam.z[1] <- 1
#
#   # Construct block toeplitz matrix from given ACF function and AR(P)
#   R = toeplitz(gam.z[1:ar.order])
#
#   # stack the gammas
#   r = gam.z[2:(ar.order+1)]
#
#   # multivariate Yule-Walker
#   phi.est = r %*% solve(R)
#
#   return(c(lam1.est, lam2.est, prob.est, phi.est))
#
# }
#
