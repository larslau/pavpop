rm(list = ls())

#
# Make package
#
#
# library(devtools)
#
# setwd('~/Dropbox/registration/')
# document()
#
# setwd("..")
# install("registration")
#
#
# setwd('~/Dropbox/registration/')
# require(registration)
#
#
# Set up data
#

n <- 50
m <- 100
nw <- 10
nknots <- 10


t <- seq(0, 1, length = m)
tw <- seq(0, 1, length = nw + 2)[2:(nw + 1)]
knots <- seq(0, 1, length = nknots + 2)[2:(nknots + 1)]

basis <- bs(t, knots = knots, Boundary.knots = c(0, 1))


# Simulate data

sigma <- 0.2
tau <- 0.5
scale <- 1
range <- 0.5

c_true <- rnorm(ncol(basis))

theta_t <- basis %*% c_true

# Covariance

C <- Brownian_cov(tw, type = 'bridge')
Cinv <- solve(C)
Csqrt <- t(chol(sigma^2 * tau^2 * C))

S <- Matern_cov(t, scale, range, 2)
Ssqrt <- t(chol(sigma^2 * S))

w <- replicate(n, Csqrt %*% rnorm(nw))[, , 1:n]


y <- array(NA, dim = c(m, n))

for (i in 1:n) {
  y[, i] <- bs(v(w[, i], t, tw), knots = knots, Boundary.knots = c(0, 1)) %*% c_true + Ssqrt %*% rnorm(m, sd = 1) +  rnorm(m, sd = sigma)
}


# Plot data
par(mfrow = c(1,2))
plot(t, theta_t, ylim = c(min(y), max(y)), type = "l", lwd = 2, lty = 2)
for (i in 1:n) lines(t, y[, i], lwd = 0.2)


tm <- yl <- list()
for (i in 1:n) {
  tm[[i]] <- t
  yl[[i]] <- y[,i]
}

# Estimate

tw <- seq(0, 1, length = 12)[2:11]
kts <- knots
tau <- scale <- range <- 1
smoothness <- 2


res <- estimate(yl, tm, kts, tw, tau, scale, range, 1)

cat('sigma:\t', sigma, '\t\t', res$sigma, '\n', 'scale:\t', scale, '\t\t', res$scale, '\n', 'range:\t', range, '\t\t', res$range, '\n', 'tau:\t', tau, '\t\t', res$tau, sep = '')


#traceback()


plot(t, theta_t, ylim = c(min(y), max(y)), type = "l", lwd = 2, lty = 2)
for (i in 1:n) lines(v(res$w[,i], t, tw), y[, i], lwd = 0.2)
lines(t, basis %*% res$c, col = 'green', lwd = 2)


# plow warps

par(mfrow = c(1,2))

plot(t, t, type = 'n')
for (i in 1:n) lines(v(w[, i], t, tw), t)

plot(t, t, type = 'n')
for (i in 1:n) lines(v(w[, i], t, tw), v(res$w[,i], t, tw))
