n <- 100
lam.tru <- 3
phi.tru <- 0.5
lag.max <- 2


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
