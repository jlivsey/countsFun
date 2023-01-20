# Purpose: impliment the ideas in note by Vladas
#          W-12-23-19-psuedo likelihood.pdf

n <- 100
k.tru <- 3
p.tru <- 1/4
phi.tru <- 0.5
lag.max <- 20


# only let phi be a free parameter
phi.seq <- seq(.4, .6, length.out = 50)
V <- c()
for(i in 1:length(phi.seq)){
  phi <- phi.seq[i]
  V[i] <- myV_negBinom(size.tru = 3,
                       prob.tru = 1/4,
                       phi.tru = 0.5,
                       size.free = 3,
                       prob.free = 1/4,
                       phi.free = phi)
}
plot(V~phi.seq)

phi.seq[which.min(V)]

# Only let size be a free parameter
size.seq <- seq(2.8, 3.2, length.out = 50)
V <- c()
for(i in 1:length(size.seq)){
  size <- size.seq[i]
  V[i] <- myV_negBinom(size.tru = 3,
                       prob.tru = 1/4,
                       phi.tru = 0.5,
                       size.free = size,
                       prob.free = 1/4,
                       phi.free = 0.5)
}
plot(V~size.seq)

size.seq[which.min(V)]


# Do for all Free parameters concurrently
phi.seq <- seq(.4, .6, length.out = 10)
size.seq <- seq(2.8, 3.2, length.out = 10)
prob.seq <- seq(.2, .3, length.out = 10)

M <- matrix(ncol = 4)
colnames(M) <- c("size", "prob", "phi", "V")
V.tru <- myV_negBinom(
  size.tru = 3,
  prob.tru = 1/4,
  phi.tru = 0.5,
  size.free = 3,
  prob.free = 1/4,
  phi.free = 0.5
)
M[1, ] <- c(3, 1/4, 0.5, V.tru)
for(k in 1:length(size.seq)){
  print(k)
  for(i in 1:length(prob.seq)){
    for(j in 1:length(phi.seq)){
      size <- size.seq[k]
      prob <- prob.seq[i]
      phi <- phi.seq[j]
      V <- myV_negBinom(
        size.tru = 3,
        prob.tru = 1/4,
        phi.tru = 0.5,
        size.free = size,
        prob.free = prob,
        phi.free = phi
      )
      M <- rbind(M, c(size, prob, phi, V))
    }
  }}

M[which.min(M[, 4]), ]









