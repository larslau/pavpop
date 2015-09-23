# TODO:
# EXPLICITLY GENERATE INVERSE OF BROWNIAN COVARIANCES (+OTHERS?!?)
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
#' @param param either one-dimensional or equal to the length of t consisting of the diagonal entries.
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
  return(tau^2 * (min(t) - is_bridge * prod(t)))
}
attr(Brownian, 'stationary') <- FALSE

#' Constant shift covariance
#'
#' Constant shift covariance. Note, the resulting covariance matrix is only positive definite if noise is added.
#' @param t two-dimensional vector of evaluation points.
#' @param param parameter vector consisting of scale parameter tau.
#' @keywords covariance
#' @export

const_cov <- function(t, param = c(scale = 1)) {
  return(matrix(param, length(t), length(t)))
}
attr(const_cov, 'discrete') <- TRUE
attr(const_cov, 'stationary') <- TRUE


#' Generate covariance matrix function from covariance function
#'
#' Function that takes a covariance function and returns a function that generates covariance matrices with the given covariance function
#' @param cov_fct covariance function.
#' @param noise logical. Should an identity matrix be added to the covariance matrix?
#' @param param standard values of parameters
#' @param ns list of arguments to make a stationary covariance locally adaptive.
#' @param ... arguments passed to cov_fct.
#' @export

#TODO: EXAMPLE

make_cov_fct <- function(cov_fct, noise = TRUE, param = NULL, inv_cov_fct = NULL, ns = NULL, ...) {
  if (!is.null(ns)) {
    if (is.null(ns$knots)) stop('ns must be a list with the argument knots.')
    if (ns$knots < 2) {
      warning('number of knots should at least be 2, ignoring argument.')
      ns <- NULL
    }
  }
  if (!is.null(attr(cov_fct, 'discrete'))) {
    if (attr(cov_fct, 'discrete')) {
      if (noise) {
        f <- function(t, param) cov_fct(t, param) + id_cov(t)
      } else {
        f <- cov_fct
      }

    } else {
      stop('attribute \'discrete\' should be NULL for non-discrete covariances.')
    }
  } else {
    if (attr(cov_fct, 'stationary')) {
      if (is.null(ns)) {
        # stationary covariance, fill rows and columns simultaneously
        f <- function (t, param) {
          m <- length(t)
          S <- diag(cov_fct(0, param, ...) + noise, m)
          if (m > 1) {
            for (i in 1:(m - 1)) {
              S[i, (i + 1):m] <- S[(i + 1):m, i] <- cov_fct(abs(t[i] - t[(1 + i):m]), param, ...)
            }
          }
          return(S)
        }
      } else {
        # ns = list(knots = 4)
        knots <- ns$knots
        f <- function (t, param) {
          ns <- spline(seq(0, 1, length = knots), c(param[1:(knots - 1)], 1), xout = t)$y
          m <- length(t)
          S <- diag(cov_fct(0, param[-(1:(knots - 1))], ...), m)
          if (m > 1) {
            for (i in 1:(m - 1)) {
              S[i, (i + 1):m] <- S[(i + 1):m, i] <- cov_fct(abs(t[i] - t[(1 + i):m]), param[-(1:(knots - 1))], ...)
            }
          }
          return(ns %*% t(ns) * S + diag(1, nrow = length(t)))
        }
      }
    } else {
      # Non-stationary covariance, fill in all entries separately
      if (is.null(ns)) {
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
      } else {
        knots <- ns$knots
        f <- function (t, param) {
          ns <- spline(seq(0, 1, length = knots), c(param[1:(knots - 1)], 1), xout = t)$y
          m <- length(t)
          S <- matrix(NA, m, m)
          for (i in 1:m) {
            for (j in i:m) {
              S[i, j] <- S[j, i] <- ns[i] * ns[j] * cov_fct(c(t[i], t[j]), param[-(1:(knots - 1))], ...)
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
    # Include covariance function evaluated with additional arguments
    eval_cov_fct <- function(t, param) cov_fct(t, param, ...)
    attributes(eval_cov_fct) <- attributes(cov_fct)
    attr(f, 'cov_fct') <- eval_cov_fct
    attr(f, 'noise') <- noise
    # If the covariance has been made non-stationary, modify
    if (!is.null(ns)) {
      attr(f, 'param') <- c(local_scale = rep(1, ns$knots - 1), param)
      ns_cov_fct <- function(t, param) {
        prod(spline(seq(0, 1, length = ns$knots), c(param[1:(ns$knots - 1)], 1), xout = t)$y) * cov_fct(abs(diff(t)), param[-(1:(ns$knots - 1))])
      }
      attr(ns_cov_fct, 'stationary') <- FALSE
      attr(f, 'cov_fct') <- ns_cov_fct
    }
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
        for (j in 1:m) {
          S[i, j] <- cov_fct(c(t_p[i], t[j]), param, ...)
        }
      }
    }
    return(S)
  }
