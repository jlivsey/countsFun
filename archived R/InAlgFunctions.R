#
# # innovations algorithm code
# innovations.algorithm <- function(acvf,n.max=length(acvf)-1){
#   # Found this onlinbe need to check it
#   # http://faculty.washington.edu/dbp/s519/R-code/innovations-algorithm.R
#   thetas <- vector(mode="list",length=n.max)
#   vs <- rep(acvf[1],n.max+1)
#   for(n in 1:n.max){
#     thetas[[n]] <- rep(0,n)
#     thetas[[n]][n] <- acvf[n+1]/vs[1]
#     if(n>1){
#       for(k in 1:(n-1)){
#         js <- 0:(k-1)
#         thetas[[n]][n-k] <- (acvf[n-k+1] - sum(thetas[[k]][k-js]*thetas[[n]][n-js]*vs[js+1]))/vs[k+1]
#       }
#     }
#     js <- 0:(n-1)
#     vs[n+1] <- vs[n+1] - sum(thetas[[n]][n-js]^2*vs[js+1])
#   }
#   return(structure(list(vs=vs,thetas=thetas)))
# }
#
#
# # Nonstationary innovations algorithm
# Innalg <- function(data, GAMMA){
#   N <- length(data)
#   x.hat <- numeric(N)
#   v <- numeric(N)
#   e <- numeric(N)
#   theta <- matrix(0, N, N)
#
#   x.hat[1] <- 0
#   v[1] <- GAMMA[1, 1]
#   e[1] <- data[1]
#
#   for (n in 1:(N-1)){
#     for (k in 0:(n-1)){
#       a <- 0
#       if (k > 0) {
#         a <- sum(theta[k, 1:k] * theta[n, 1:k] * v[1:k])
#       }
#
#       theta[n, k+1] <- (1/v[k+1]) * (GAMMA[n+1, k+1] - a)
#     }
#     if(n<N){
#       x.hat[n+1] <- sum(theta[n, 1:n] * (data[1:n] - x.hat[1:n]))
#       v[n+1] <- GAMMA[n+1, n+1] - sum(theta[n, 1:n]^2 * v[1:n])
#       e[n+1] <- data[n+1] - x.hat[n+1]
#     }
#   }
#
#   return(list(x.hat = x.hat,
#               theta = theta,
#               v     = v    ,
#               e     = e     ))
# }
