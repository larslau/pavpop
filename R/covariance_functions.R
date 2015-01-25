# TODO:
# EXPLICITLY GENERATE INVERSE OF BROWNIAN COVARIANCES (+OTHERS?!)



#' Generate Brownian covariances
#'
#' This function generates a Brownian motion/bridge covariance matrix corresponding to specified evaluation points.
#' @param t evaluation points.
#' @param tau scale parameter.
#' @param type type of covariance, either 'motion' or 'bridge'.
#' @keywords covariance
#' @export
#' @examples
#' t <- seq(0, 1, length = 10)
#' Brownian_cov(t, 1)

Brownian_cov <- function(t, tau, type = 'motion') {
	m <- length(t)
	C <- matrix(NA, m, m)
	is_bridge <- type == 'bridge'
	for (i in 1:m) {
		for (j in i:m) {
			C[i, j] <- C[j, i] <- tau^2 * (min(t[i], t[j]) - is_bridge * t[i] * t[j])
		}
	}
	return(C)
}

#' Generate Brownian motion covariances
#'
#' This function generates a Brownian motion covariance matrix corresponding to specified evaluation points.
#' @param t evaluation points.
#' @param tau scale parameter.
#' @keywords covariance
#' @export
#' @examples
#' t <- seq(0, 1, length = 10)
#' Brownian_motion_cov_fast(t)

#TODO: UPDATE

Brownian_motion_cov_fast <- function(t, tau = 1) {
	m <- length(t)
	C <- matrix(NA, m, m)

	for (i in m:1) {
		C[i, 1:m] <- C[1:m, i] <- t[i]
	}
	return(C)
}

# require(Rmpfr)
# x <- t[[2]]
# jn(2, x)
# yn(2, x)
# BesselK(x, as.integer(2))
#
# Matern <- function (d, scale = 1, range = 1, alpha = 1/range, smoothness = 0.5,
#           nu = smoothness, phi = scale)
# {
#   if (any(d < 0))
#     stop("distance argument must be nonnegative")
#   d <- d * alpha
#   d[d == 0] <- 1e-10
#   con <- (2^(nu - 1)) * gamma(nu)
#   con <- 1/con
#   return(phi * con * (d^nu) * BesselK(d, nu, expon.scaled = FALSE))
# }

#' Generate Matern plus measurement noise covariances
#'
#' This function generates a Matern motion covariance matrix corresponding to specified evaluation points.
#' @param t evaluation points.
#' @param param parameter vector consisting of scale, range and smoothness.
#' @param noise logical, should a diagonal matrix be added to the Matern covariance?
#' @keywords covariance
#' @export
#' @examples
#' t <- seq(0, 1, length = 10)
#' Matern_cov(t, param = c(1, 1, 1/2))

Matern_cov <- function(t, param = c(scale = 1, range = 1, smoothness = 2), noise = TRUE) {
  scale <- param[1]
  range <- param[2]
  smoothness <- param[3]
  m <- length(t)
  S <- diag(x = Matern(0, scale = scale, range = range, smoothness = smoothness) + noise, nrow = m)
  i <- 1
#   while (i < m) {
#     S[cbind(1:(m - i), (1 + i):m)] <- S[cbind((1 + i):m, 1:(m - i))] <- Matern(abs(t[1:(m - i)] - t[(1 + i):m]), scale = scale, range = range, smoothness = smoothness)
#     i <- i + 1
#   }
  for (i in 1:(m - 1)) {
    S[i, (i + 1):m] <- S[(i + 1):m, i] <- Matern(abs(t[[i]] - t[(1 + i):m]), scale = scale, range = range, smoothness = smoothness)
  }
  return(S)
}

