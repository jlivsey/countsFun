n <- 100
lam.tru <- 3
phi.tru <- 0.5
lag.max <- 20

library(orthopolynom)
polys <- hermite.he.polynomials(100)

myV <- function(lam.tru, phi.tru, lam.free, phi.free){

  lam.free <- 3
  phi.tru <- 1/2

  hc <- HermCoef(lam.tru)
  gam.z <- ARMAacf(ar = phi.tru, lag.max = n)
  gam.x <- CountACVF(h = 1:n, myacf = gam.z, g = hc)
  gam.x[1] <- 1

  gam.x.tru <- gam.x

  # STOPPED HERE ^^^^

  gvec = rep(NA, lag.max+1) # storage
  for(i in 0:lag.max) gvec[i+1] <- g(lam = lam.free, k = i, polys = polys)
  gam.z <- ARMAacf(ar = phi.free, lag.max = n)
  gam.x <- gam_x(h = 0:n, lam = lam.free, myacf = gam.z, gvec = gvec)
  gam.x[1] <- 1

  gam.x.free <- gam.x

  gam0 <- gam.x.tru[-1]
  gam1 <- gam.x.free[-1]
  G0 <- toeplitz(gam.x.tru[-(n+1)])
  G1 <- toeplitz(gam.x.free[-(n-1)])

  term1 <- lam.tru
  term2 <- -2 * lam.tru * t(gam0) %*% solve(G1) %*% gam1
  term3 <- lam.tru * t(gam1) %*% solve(G1) %*% G0 %*% solve(G1) %*% gam1

  V <- term1 + term2 + term3

  return(V)
}


phi.seq <- seq(.4, .6, length.out = 50)
V <- c()
for(i in 1:length(phi.seq)){
  phi <- phi.seq[i]
  V[i] <- myV(lam.tru = 1, phi.tru = .5, lam.free = 1, phi.free = phi)
}
plot(V~phi.seq)

phi.seq[which.min(V)]


lam.seq <- seq(.8, 1.2, length.out = 50)
V <- c()
for(i in 1:length(lam.seq)){
  lam <- lam.seq[i]
  V[i] <- myV(lam.tru = 1, phi.tru = .5, lam.free = lam, phi.free = .5)
}
plot(V~lam.seq)

lam.seq[which.min(V)]


M <- matrix(ncol = 3)
colnames(M) <- c("lam", "phi", "V")
V.tru <- myV(lam.tru = 1, phi.tru = .5, lam.free = 1, phi.free = .5)
M[1, ] <- c(1, .5, V.tru)
for(i in 1:length(lam.seq)){
  print(i)
  for(j in 1:length(phi.seq)){
    lam <- lam.seq[i]
    phi <- phi.seq[j]
    V <- myV(lam.tru = 1, phi.tru = .5, lam.free = lam, phi.free = phi)
    M <- rbind(M, c(lam, phi, V))
  }
}

M[which.min(M[, 3]), ]




lam.seq <- seq(.8, 1.2, length.out = 25)
phi.seq <- seq(-.6, -.4, length.out = 25)
M <- matrix(ncol = 3)
colnames(M) <- c("lam", "phi", "V")
V.tru <- myV(lam.tru = 1, phi.tru = -.5, lam.free = 1, phi.free = -.5)
M[1, ] <- c(1, -.5, V.tru)
for(i in 1:length(lam.seq)){
  print(i)
  for(j in 1:length(phi.seq)){
    lam <- lam.seq[i]
    phi <- phi.seq[j]
    V <- myV(lam.tru = 1, phi.tru = -.5, lam.free = lam, phi.free = phi)
    M <- rbind(M, c(lam, phi, V))
  }
}
M[which.min(M[, 3]), ]

library(plotly)
p <- plot_ly(x = M[,1], y = M[,2], z = M[,3]) %>% add_surface()
p








