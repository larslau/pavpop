# TODO:
# EXPLICITLY GENERATE INVERSE OF BROWNIAN COVARIANCES (+OTHERS?!)
# INCLUDE CODE FOR MATERN
# CHECK IF ATTRIBUTE CAN BE SET SMARTER
# CHECK IF SYMMETRIC MATRICES ARE MORE EFFICIENT

#' Zero covariance function
#'
#' Returns a matrix with zeros with rows and columns equal to the length of t.
#' @param t evaluation points.
#' @param param not used.
#' @keywords covariance
#' @export

zero_cov <- function(t, param = NULL) {
  return(matrix(0, length(t), length(t)))
}
attr(zero_cov, 'discrete') <- TRUE

#' Identity covariance function
#'
#' Returns a diagonal matrix with diagonal given by param.
#' @param t evaluation points.
#' @param param one-dimensional consisting of the diagonal entry.
#' @keywords covariance
#' @export

id_cov <- function(t, param = c(scale = 1)) {
  return(diag(param, length(t)))
}
attr(id_cov, 'discrete') <- TRUE


#' Matern covariance function
#'
#' Functional form of Matern covariance function. Code adapted from the fields package.
#' @param d evaluation points.
#' @param param parameter vector consisting of scale, range and smoothness.
#' @keywords covariance
#' @note Method taken from the \code{fields} package.
#' @export
#' @examples
#' Matern(seq(0, 1, length = 10), param = c(scale = 1, range = 0.5, smoothness = 2))
#' Matern(seq(0, 1, length = 10), param = c(scale = 1, range = 1, smoothness = 2))


Matern <- function(d, param = c(scale = 1, range = 1, smoothness = 2)) {
  scale <- param[1]
  range <- param[2]
  smoothness <- param[3]
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d <- d / range
  d[d == 0] <- 1e-10
  con <- (2^(smoothness - 1)) * gamma(smoothness)
  con <- 1 / con
  return(scale * con * (d^smoothness) * besselK(d, smoothness))
}
attr(Matern, 'stationary') <- TRUE

#' Brownian covariance functions
#'
#' Functional form of Brownian covariance functions. Brownian bridge is on [0, 1].
#' @param t two-dimensional vector of evaluation points.
#' @param param parameter vector consisting of scale parameter tau.
#' @param type type of covariance, either 'motion' or 'bridge'.
#' @keywords covariance
#' @export
#' @examples
#' Brownian(t = c(1, 1))
#' Brownian(t = c(1, 1), type = 'bridge')

Brownian <- function(t, param = c(tau = 1), type = 'motion') {
  tau <- param[[1]]
  is_bridge <- type == 'bridge'
  return(tau^2 * (min(t[1], t[2]) - is_bridge * t[1] * t[2]))
}
attr(Brownian, 'stationary') <- FALSE

#' Generate covariance matrix function from covariance function
#'
#' Function that takes a covariance function and returns a function that generates covariance matrices with the given covariance function
#' @param cov_fct covariance function.
#' @param noise logical. Should an identity matrix be added to the covariance matrix?
#' @param ... arguments passed to cov_fct.
#' @export


make_cov_fct <- function(cov_fct, noise = TRUE, param = NULL, inv_cov_fct = NULL, ...) {
  if (!is.null(attr(cov_fct, 'discrete'))) {
    if (attr(cov_fct, 'discrete')) {
      f <- cov_fct
    } else {
      stop('attribute \'discrete\' should be NULL for non-discrete covariances.')
    }
  } else {
    if (attr(cov_fct, 'stationary')) {
      # stationary covariance, fill rows and columns simultaneously
      f <- function (t, param) {
        m <- length(t)
        S <- diag(cov_fct(0, param, ...) + noise, m)
        for (i in 1:(m - 1)) {
          S[i, (i + 1):m] <- S[(i + 1):m, i] <- cov_fct(abs(t[i] - t[(1 + i):m]), param, ...)
        }
        return(S)
      }
    } else {
      # Non-stationary covariance, fill in all entries separately
      f <- function (t, param) {
        m <- length(t)
        S <- matrix(NA, m, m)
        for (i in 1:m) {
          for (j in i:m) {
            S[i, j] <- S[j, i] <- cov_fct(c(t[i], t[j]), param, ...)
          }
        }
        if (noise) diag(S) <- diag(S) + 1
        return(S)
      }
    }
  }
  if (is.null(param)) param <- formals(cov_fct)$param
  attr(f, 'param') <- param
  # Set solve method
  attr(f, 'inv_cov_fct') <- inv_cov_fct
  # Set the scale parameter
  attr(f, 'scale') <- ifelse(is.null(formals(cov_fct)$param$scale), NA, which(names(formals(cov_fct)$param) == 'scale') - 1)
  attr(f, 'cov_fct') <- cov_fct
  return(f)
}

#' Rectangular evaluation of covariance functions
#'
#' Generate rectangular evaluations of covariance functions that are typically used for prediction purposes.
#' @param t observation points.
#' @param t_p new "prediction" points.
#' @param cov_fct covariance function.
#' @param ... arguments passed to cov_fct.
#'
cov_rect <- function(t, t_p, cov_fct, param, ...) {
  m_p <- length(t_p)
  m <- length(t)
  S <- matrix(NA, m_p, m)

  if (attr(cov_fct, 'stationary')) {
    # stationary covariance, fill rows and columns simultaneously
    for (i in 1:m) {
      S[, i] <- cov_fct(abs(t[i] - t_p), param, ...)
    }
  } else {
    # Non-stationary covariance, fill in all entries separately
    for (i in 1:m_p) {
      for (j in i:m) {
        S[i, j] <- S[j, i] <- cov_fct(c(t_p[i], t[j]), param, ...)
      }
    }
  }
  return(S)
}
