# rm(list = ls())
require(Matrix)

# set.seed(seed = 12123)

n <- 50
m <- 50
mw <- 10
nknots <- 10


t <- seq(0, 1, length = m)
tw <- seq(0, 1, length = mw + 2)[2:(mw + 1)]
kts <- seq(0, 1, length = nknots + 2)[2:(nknots + 1)]

basis <- bs(t, knots = kts, Boundary.knots = c(0, 1))


# Simulate data

sigma <- 0.1
tau <- 0.1 / sigma
scale <- 0.1 / sigma^2
range <- 0.3
smoothness <- 2


c_true <- rnorm(ncol(basis))

theta_t <- basis %*% c_true

# Covariance

C <- Brownian_cov(tw, tau, type = 'bridge')
Cinv <- solve(C)
Csqrt <- t(chol(sigma^2 * C))

S <- Matern_cov(t, c(scale, range, smoothness), noise = TRUE)
Ssqrt <- t(chol(sigma^2 * S))

w <- replicate(n, Csqrt %*% rnorm(mw))[, , 1:n]
for (i in 1:n) w[, i] <- make_homeo(w[,i], tw)

y <- array(NA, dim = c(m, n))

for (i in 1:n) {
  y[, i] <- bs(v(w[, i], t, tw), knots = kts, Boundary.knots = c(0, 1)) %*% c_true + Ssqrt %*% rnorm(m, sd = 1)
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

t <- tm
y <- yl


#Estimate mean curve
basis <- bs(x = t[[1]], knots = kts)
c <- spline_weights(y, t, array(0, dim = c(mw, n)), tw, kronecker(diag(n), solve(S)), kts)

lines(t[[1]], basis %*% c, lwd = 2, col = 'red')


#res <- estimate(yl, tm, kts, tw, tau, scale, range, 1)
#TODO: PROBLEM WITH SMOOTHNESS >= 2!!!
#res <- estimate_generic(y, t, c(scale, range, smoothness), function(t, param) Matern_cov(t, param, noise = TRUE), tau, function(t, param) Brownian_cov(t, param, type = 'bridge'), kts, tw, iter = c(10, 2), like_optim_control = list(lower = c(1e-3, 1e-3, 0.5, 1e-3), upper = c(1e3, 3, 4, 20), ndeps = c(1e-3, 1e-3, 0.1, 1e-3)))
system.time({
res <- estimate_generic(y, t, c(10, 1, 1), function(t, param) Matern_cov(t, param, noise = TRUE), 1, function(t, param) Brownian_cov(t, param, type = 'bridge'), kts, tw, iter = c(3, 3), use_warp_gradient = FALSE, homeomorphism = 'soft', like_optim_control = list(lower = c(1e-3, 1e-3, 0.5, 1e-3), upper = c(1e3, 3, 4, 20), ndeps = c(1e-3, 1e-3, 0.1, 1e-3)))
})
#cat('sigma:\t', sigma, '\t\t', res$sigma, '\n', 'scale:\t', scale, '\t\t', res$scale, '\n', 'range:\t', range, '\t\t', res$range, '\n', 'tau:\t', tau, '\t\t', res$tau, sep = '')
cat('sigma:\t', sigma, '\t\t', res$sigma, '\n', 'scale:\t', scale, '\t\t', res$amp_cov_par[1], '\n', 'range:\t', range, '\t\t', res$amp_cov_par[2], '\n', 'smoothness:\t', smoothness, '\t\t', res$amp_cov_par[3], '\n', 'tau:\t', tau, '\t\t', res$warp_cov_par[1], sep = '')

#traceback()

t_p <- t[[1]]
plot(t_p, theta_t, ylim = range(sapply(y, range)), type = "l", lwd = 2, lty = 2)
for (i in 1:n) lines(v(res$w[,i], t[[i]], tw), y[[i]], lwd = 0.2)
lines(t_p, basis %*% res$c, col = 'green', lwd = 2)


# plot warps

par(mfrow = c(1,2))

plot(t_p, t_p, type = 'n')
for (i in 1:n) lines(v(w[, i], t[[i]], tw), t[[i]])

plot(t_p, t_p, type = 'n')
for (i in 1:n) lines(v(w[, i], t[[i]], tw), v(res$w[,i], t[[i]], tw))

cat(mean((res$w-w)^2), '\n', sep = '\t')






# plot(t_p, t_p, type = 'n')
# for (i in 1:n) lines(t[[i]], v(res$w[,i], t[[i]], tw), lwd = 0.4)


#
# Predict warps
#
warp_cov_par <- tau
warp_cov_fct <- function(t, param) Brownian_cov(t, param, type = 'bridge')
amp_cov_par <- c(scale, range, smoothness)
amp_cov_fct <- function(t, param) Matern_cov(t, param, noise = FALSE)
# Build amplitude covariances and inverse covariances
S <- Ainv <- list()
for (i in 1:n) {
  S[[i]] <- amp_cov_fct(tm[[i]], amp_cov_par)
  Ainv[[i]] <- chol2inv(chol(S[[i]]))
}


for (iter in c(1, 10, 20)) {
  for (beta0 in c(5, 10, 15)) {
    w_pred <- predict_warp_pyramid(array(0, dim = c(mw, n)), yl, Ainv, tm, tw, kts, warp_cov_par, warp_cov_fct, iter = c(iter, 10), plevels = 10,  beta0 = beta0)
    cat(beta0, iter, mean((w_pred-w)^2), '\n', sep = '\t')
    plot(w, w_pred, xlim = range(w, w_pred), ylim = range(w, w_pred))
    abline(a = 0, b = 1)
  }
}


w_pred <- predict_warp_pyramid(array(0, dim = c(mw, n)), yl, Ainv, tm, tw, kts, warp_cov_par, warp_cov_fct, iter = c(10, 10), plevels = 5,  beta0 = 15)
cat(15, 10, mean((w_pred-w)^2), '\n', sep = '\t')
plot(w, w_pred, xlim = range(w, w_pred), ylim = range(w, w_pred))
abline(a = 0, b = 1)

plot(w_pred[,1], ylim = range(w_pred), type = 'n')
for (i in 1:n) lines(w_pred[,i], col = rainbow(n)[i])


par(mfrow = c(1,3))
plot(t, theta_t, ylim = c(min(y), max(y)), type = "l", lwd = 2, lty = 2)
for (i in 1:n) lines(t, y[, i], lwd = 0.2)

plot(t, theta_t, ylim = c(min(y), max(y)), type = "l", lwd = 2, lty = 2)
for (i in 1:n) lines(v(res$w[,i], t, tw), y[, i], lwd = 0.2)
lines(t, basis %*% res$c, col = 'green', lwd = 2)

plot(t, theta_t, ylim = c(min(y), max(y)), type = "l", lwd = 2, lty = 2)
for (i in 1:n) lines(v(w_pred[,i], t, tw), y[, i], lwd = 0.2)

par(mfrow = c(1,4))

plot(t, t, type = 'n')
for (i in 1:n) lines(v(w[, i], t, tw), t)

plot(t, t, type = 'n')
for (i in 1:n) lines(v(w[, i], t, tw), v(make_homeo(res$w[,i], tw), t, tw))

plot(t, t, type = 'n')
for (i in 1:n) lines(v(w[, i], t, tw), v(res$w[,i], t, tw))


plot(t, t, type = 'n')
for (i in 1:n) lines(v(w[, i], t, tw), v(w_pred[,i], t, tw))






