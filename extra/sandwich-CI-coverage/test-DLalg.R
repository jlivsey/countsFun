N <- 100
x <- arima.sim(list(ar = c(1.7, -1.2, .4)), n = N+1)
x <- x / sd(x)
g <- acf(x, lag.max = N, plot = FALSE)$acf[2:(N+1)]
out <- DLalg(x = x[1:N], g = g)
out2 <- DLalg2xhat(out)
plot(x)
lines(out2$xhat, col = 2)
out$v
