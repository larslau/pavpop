dat <- read.csv('~/Dropbox/vækst/boys.csv')
dat <- dat[order(dat$Age),]


pers <- unique(dat$Person)

n <- length(pers)

y_height <- t_height <- list()
y <- t <- list()

k <- 1
for (p in pers) {
  idx <- dat$Person == p
  y_height[[k]] <- dat$Heightmean[idx]
  ind <- !is.na(y_height[[k]])
  y_height[[k]] <- y_height[[k]][ind]
  t_height[[k]] <- dat$Age[idx][ind]

  m <- length(t_height[[k]])
  t[[k]] <- (t_height[[k]][1:(m - 1)] + t_height[[k]][2:m]) / 2
  y[[k]] <- (y_height[[k]][2:m] - y_height[[k]][1:(m-1)]) / (t_height[[k]][2:m] - t_height[[k]][1:(m-1)])
  k <- k + 1
}

# Remove single observations
ind <- c()
for (i in 1:n) {
  if (length(y[[i]]) < 10) ind <- c(ind, i)
}

y[ind] <- t[ind] <- NULL


n <- length(y)
m <- sapply(y, length)


# Normalize t (should be done more careful)
t_range <- range(sapply(t, range)) + c(-1, 1)
t <- lapply(t, function(x) (x - t_range[1]) / diff(t_range))


mw <- 2
nknots <- 5
tw <- seq(0, 1, length = mw + 2)[2:(mw + 1)]
kts <- seq(0, 1, length = nknots + 2)[2:(nknots + 1)]

# Plot data
par(mfrow = c(1, 2))
plot(0, 0, xlim = c(0, 1), ylim = range(sapply(y, range)), type = "n")
for (i in 1:n) lines(t[[i]], y[[i]], lwd = 0.2)

scale <- 1
range <- 0.1
smoothness <- 0.5
tau <- 1


amp_cov_par <- c(scale, range)
amp_cov_fct <- function(t, param) Matern_cov(t, c(param, smoothness), noise = TRUE)

warp_cov_par <- tau
warp_cov_fct <- function(t, param) Brownian_cov(t, param, type = 'bridge')

res <- estimate_generic(y, t, amp_cov_par, amp_cov_fct, warp_cov_par, warp_cov_fct, kts, tw, iter = c(5, 2), like_optim_control = list(lower = c(1e-3, 1e-5, 1e-4), upper = rep(1e5, 3, 10)))


#
# Plot
#
t_p <- seq(0, 1, length = 100)

w <- array(0, dim = c(mw, n))
S <- Ainv <- list()
for (i in 1:n) {
  S[[i]] <- amp_cov_fct(t[[i]], amp_cov_par)
  Ainv[[i]] <- chol2inv(chol(S[[i]]))
}
basis <- bs(x = t_p, knots = kts)

par(mfrow = c(1, 2))
plot(0, 0, xlim = c(0, 1), ylim = range(sapply(y, range)), type = "n")
for (i in 1:n) lines(t[[i]], y[[i]], lwd = 0.2)

c <- spline_weights(y, t, w, tw, bdiag(Ainv), kts)
lines(t_p, basis %*% c, col = 'red', lwd = 2)



plot(0, 0, xlim = c(0, 1), ylim = range(sapply(y, range)), type = "l", lwd = 2, lty = 2)
for (i in 1:n) lines(v(res$w[,i], t[[i]], tw), y[[i]], lwd = 0.2)
lines(t_p, basis %*% c, col = 'red', lwd = 2)
lines(t_p, basis %*% res$c, col = 'green', lwd = 2)


par(mfrow = c(1, 1))
plot(t_p, t_p, type = 'l', lty = 2)
for (i in 1:n) lines(t[[i]], v(res$w[, i], t[[i]], tw), lwd = 0.5)
