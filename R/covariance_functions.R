# TODO:
# EXPLICITLY GENERATE INVERSE OF BROWNIAN COVARIANCES (+OTHERS?!?)
# CHECK IF ATTRIBUTE CAN BE SET SMARTER
# CHECK IF SYMMETRIC MATRICES ARE MORE EFFICIENT

#' Discrete covariance-matrix functions
#'
#' Discrete covariance-matrix functions return a covariance matrix that does not depend
#' on evaluation points \code{t}. Possibilities are \code{zero_cov}, \code{id_cov},
#' \code{diag_cov} and \code{unstr_cov}. See details for further information.
#'
#' The discrete covariance-matrix functions return \eqn{m x m} covariance matrices
#' where \eqn{m} is the length of \code{t}. They are discrete in the sense that they
#' do not depend on the specific evaluation points \code{t}. These covariance matrix
#' functions are mainly used to model latent variables.
#'
#' \code{zero_cov} returns an \eqn{m x m} zero matrix.
#'
#' \code{id_cov} returns an \eqn{m x m} identity matrix.
#'
#' \code{const_cov} returns an \eqn{m x m} matrix consisting of all ones.
#' \strong{Note:} this covariance matrix has rank 1 and is thus generally
#' not positive definite.
#'
#' \code{diag_cov} returns an \eqn{m x m} diagonal matrix with \code{param} on the diagonal.
#'
#' \code{unstr_cov} returns an unstructured \eqn{m x m} covariance matrix with the diagonal
#' given by the first \eqn{m} elements in \code{param}, and the remaining filling
#' the upper and lower triangles. If the supplied parameters does not specify a positive
#' definite matrix, the function tries to return the nearest positive definite matrix (see
#' \code{\link[Matrix]{nearPD}}). This may cause \code{unstr_cov} to be slow if \eqn{m} is
#' not small.
#'
#' @param t evaluation points.
#' @param param parameters for the covariance matrix.
#' @return Covariance matrix of dimension \eqn{m x m} where \eqn{m} is the length
#' of \code{t}.
#' @keywords covariance
#' @seealso \code{link{make_cov_fct}}
#' @examples
#' # Evaluation points
#' t <- 0:1
#' # Generate zero, identity and constant covariance matrices
#' zero_cov(t)
#' id_cov(t)
#' const_cov(t)
#'
#' # Generate diagonal covariance
#' diag_cov(t, param = 1:3)
#'
#' # Generate unstructured covariance matrix
#' unstr_cov(t, param = c(1, 1, 0.5))
#'
#' # Generate unstructured covariance matrix with parameters
#' # that will not produce a positive matrix
#' (C <- unstr_cov(t, param = c(1, 1, 1.1)))
#' det(C)
#'
#' @name discrete_cov
NULL
#> NULL

#' @rdname discrete_cov
#' @export

zero_cov <- function(t, param = NULL) {
  if (!is.null(param)) warning('Parameters ignored.')
  return(matrix(0, length(t), length(t)))
}
attr(zero_cov, 'discrete') <- TRUE
attr(zero_cov, 'type') <- 'zero_cov'

#' @rdname discrete_cov
#' @export

id_cov <- function(t, param = c(scale = 1)) {
  return(diag(x = param, length(t)))
}
attr(id_cov, 'discrete') <- TRUE
attr(id_cov, 'type') <- 'id_cov'

#' @rdname discrete_cov
#' @export

const_cov <- function(t, param = c(scale = 1)) {
  return(matrix(param, length(t), length(t)))
}
attr(const_cov, 'discrete') <- TRUE
attr(const_cov, 'type') <- 'const_cov'

#' @rdname discrete_cov
#' @export

diag_cov <- function(t, param) {
  if (length(param) != length(t)) stop('Number of parameters must equal number of observation points.')
  return(diag(param, length(t)))
}
attr(diag_cov, 'discrete') <- TRUE
attr(diag_cov, 'type') <- 'diag_cov'

#' @rdname discrete_cov
#' @export

unstr_cov <- function(t, param) {
  m <- length(t)
  if (length(param) != m * (m + 1) / 2) stop('Unstructured covariances of dimension ', m, ' x ', m, ' must have ', m * (m + 1) / 2, ' parameters.')
  C <- matrix(NA, m, m)
  diag(C) <- param[1:m]
  C[lower.tri(C)] <- C[upper.tri(C)] <- param[-(1:m)]
  C <- as.matrix(Matrix::nearPD(C)$mat)
  return(C)
}
attr(unstr_cov, 'discrete') <- TRUE
attr(unstr_cov, 'type') <- 'unstr_cov'

#' Matern covariance function
#'
#' Functional form of Matern covariance function. Code adapted from the \code{fields} package.
#' @param d distance between two points.
#' @param param parameter vector consisting of scale, range and smoothness.
#' @keywords covariance
#' @note Method taken from the \code{fields} package.
#' @seealso \code{link{make_cov_fct}}
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
attr(Matern, 'type') <- 'Matern'

#' Brownian covariance functions
#'
#' Functional form of Brownian covariance functions. Brownian bridge is on [0, 1].
#' @param t two-dimensional vector of evaluation points.
#' @param param parameter vector consisting of scale parameter tau.
#' @param type type of covariance, either 'motion' or 'bridge'.
#' @param range interval of the process. If \code{type = 'motion'} only the first point
#' is used.
#' @keywords covariance
#' @seealso \code{link{make_cov_fct}}
#' @export
#' @examples
#' Brownian(t = c(1, 1))
#' Brownian(t = c(1, 1), type = 'bridge')

Brownian <- function(t, param = c(tau = 1), type = 'motion', range = c(0, 1)) {
  tau <- param[[1]]
  is_bridge <- type == 'bridge'
  result <- tau^2 * (min(t) - range[1])
  if (is_bridge) result <- result * (range[2] - max(t)) / diff(range)
  return(result)
}
attr(Brownian, 'stationary') <- FALSE
attr(Brownian, 'type') <- 'Brownian'

#' Generate covariance matrix function from covariance function
#'
#' Function that takes a covariance function
#' and returns a function that generates covariance matrices
#' according to the given covariance function.
#'
#' @param cov_fct covariance function.
#' @param noise logical. Should an identity matrix be added to the covariance matrix?
#' @param param standard values of parameters
#' @param ns list of arguments to make a stationary covariance locally adaptive.
#' The list has takes the following entries \code{knots} (mandatory) the number of parameters
#' to control the non-stationary, \code{fixed} (optional) vector of fixed knots with reference
#' weight of 1, if not supplied the last knot is chosen. See examples and Details.
#' @param ... arguments passed to cov_fct.
#' @export
#' @examples
#' # Generate observation points
#' t <- seq(0, 1, length = 3)
#'
#' # Generate covariances
#' matern_cov <- make_cov_fct(Matern, noise = FALSE)
#' bm_cov <- make_cov_fct(Brownian, noise = FALSE)
#' bb_cov <- make_cov_fct(Brownian, noise = FALSE, type = 'bridge')
#'
#' # Evaluate covariance matrices
#' matern_cov(t, param = c(1, 1, 1))
#' bm_cov(t, param = c(1, 1, 1))
#' bb_cov(t, param = c(1, 1, 1))
#'
#' # Plot covariance matrices
#' t <- seq(0, 1, length = 30)
#'
#' persp(t, t, matern_cov(t, param = c(1, 0.5, 1)),
#'       theta = -30, phi = 30, ticktype = 'd', zlab = '',
#'       col = 'lightblue', shade = 0.2,
#'       main = 'Matern covariance (1, 0.5, 1)', zlim = c(0, 1))
#'
#' persp(t, t, matern_cov(t, param = c(1, 0.5, 2)),
#'       theta = -30, phi = 30, ticktype = 'd', zlab = '',
#'       col = 'lightblue', shade = 0.2,
#'       main = 'Matern covariance (1, 0.5, 2)', zlim = c(0, 1))
#'
#' persp(t, t, matern_cov(t, param = c(1, 0.1, 2)),
#'       theta = -30, phi = 30, ticktype = 'd', zlab = '',
#'       col = 'lightblue', shade = 0.2,
#'       main = 'Matern covariance (1, 0.1, 2)', zlim = c(0, 1))
#'
#' persp(t, t, bm_cov(t, param = 1), theta = -30, phi = 30,
#' ticktype = 'd', zlab = '', col = 'lightblue', shade = 0.2,
#' main = 'Brownian motion covariance (tau = 1)', zlim = c(0, 1))
#'
#' persp(t, t, bb_cov(t, param = 2), theta = -30, phi = 30,
#' ticktype = 'd', zlab = '', col = 'lightblue', shade = 0.2,
#' main = 'Brownian bridge covariance (tau = 2)', zlim = c(0, 1))
#'
#' # Make covariance non-stationary
#' matern_cov_ns <- make_cov_fct(Matern, noise = FALSE,
#'                               ns = list(knots = 3, range = c(0, 1)))
#' # Use mid-point reference instead of last
#' matern_cov_ns_mid <- make_cov_fct(Matern, noise = FALSE,
#'                               ns = list(knots = 3, fixed = 2,
#'                                         range = c(0, 1)))
#'
#' # Original covariance
#' persp(t, t, matern_cov_ns(t, param = c(1, 1, 1, 0.3, 2)),
#'       theta = -30, phi = 30, ticktype = 'd', zlab = '',
#'       col = 'lightblue', shade = 0.2,
#'       main = 'Matern covariance (1, 0.1, 2)', zlim = c(0, 1))
#'
#' # Modified covariances
#' persp(t, t, matern_cov_ns(t, param = c(0.5, 0.7, 1, 0.3, 2)),
#'       theta = -30, phi = 30, ticktype = 'd', zlab = '',
#'       col = 'lightblue', shade = 0.2,
#'       main = 'Non-stationary Matern covariance (1, 0.1, 2)', zlim = c(0, 1.1))
#'
#' persp(t, t, matern_cov_ns(t, param = c(1, 0.7, 1, 0.3, 2)),
#'       theta = -30, phi = 30, ticktype = 'd', zlab = '',
#'       col = 'lightblue', shade = 0.2,
#'       main = 'Non-stationary Matern covariance (1, 0.1, 2)', zlim = c(0, 1.1))
#'
#' persp(t, t, matern_cov_ns_mid(t, param = c(1, 0.7, 1, 0.3, 2)),
#'       theta = -30, phi = 30, ticktype = 'd', zlab = '',
#'       col = 'lightblue', shade = 0.2,
#'       main = 'Non-stationary Matern covariance (1, 0.1, 2)', zlim = c(0, 1.1))

make_cov_fct <- function(cov_fct, noise = TRUE, param = NULL, inv_cov_fct = NULL, ns = NULL, ...) {
  # Check ns related things
  if (!is.null(ns)) {
    if (is.null(ns$knots)) stop('ns must be a list with the argument knots.')
    if (ns$knots < 2) {
      warning('number of knots should at least be 2, ignoring argument.')
      ns <- NULL
    }
    if (is.null(ns$range)) {
      warning('No range specified using c(0, 1).')
      ns$range <- c(0, 1)
    }
  }

  if (!is.null(attr(cov_fct, 'discrete'))) {
    # Discrete covariances
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
    # Non-discrete covariances
    if (attr(cov_fct, 'stationary')) {
      # Faster construction of stationary covariances
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
        # Stationary covariance made non-stationary
        knots <- ns$knots
        fixed <- ifelse(!is.null(ns$fixed), ns$fixed, max(knots))
        m_ns <- knots - length(fixed)
        f <- function (t, param) {
          ns_param <- rep(1, length = knots)
          ns_param[-fixed] <- param[1:m_ns]
          ns <- spline(seq(ns$range[1], ns$range[2], length = knots), ns_param, xout = t)$y
          m <- length(t)
          S <- diag(cov_fct(0, param[-(1:m_ns)], ...), m)
          if (m > 1) {
            for (i in 1:(m - 1)) {
              S[i, (i + 1):m] <- S[(i + 1):m, i] <- cov_fct(abs(t[i] - t[(1 + i):m]), param[-(1:m_ns)], ...)
            }
          }
          S <- ns %*% t(ns) * S
          if (noise) diag(S) <- diag(S) + 1
          return(S)
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
        fixed <- ifelse(!is.null(ns$fixed), ns$fixed, max(knots))
        m_ns <- knots - length(fixed)
        f <- function (t, param) {
          ns_param <- rep(1, length = knots)
          ns_param[-fixed] <- param[1:m_ns]
          ns <- spline(seq(ns$range[1], ns$range[2], length = knots), ns_param, xout = t)$y

          m <- length(t)
          S <- matrix(NA, m, m)
          for (i in 1:m) {
            for (j in i:m) {
              S[i, j] <- S[j, i] <- ns[i] * ns[j] * cov_fct(c(t[i], t[j]), param[-(1:m_ns)], ...)
            }
          }
          if (noise) diag(S) <- diag(S) + 1
          return(S)
        }
      }
    }
  }
  if (is.null(param)) param <- formals(cov_fct)$param
  attr(f, 'param') <- param

  # Set type
  attr(f, 'type') <- attr(cov_fct, 'type')
  if (attr(f, 'type') == 'Brownian') {
    args <- list(...)
    exist <- "type" %in% names(args)
    if (exist) {
      attr(f, 'type') <- attr(cov_fct, 'type') <- paste(attr(cov_fct, 'type'), args$type)
    } else {
      attr(f, 'type') <- attr(cov_fct, 'type') <- 'Brownian motion'
    }
  }
  # Set solve method
  attr(f, 'inv_cov_fct') <- inv_cov_fct
  # Set the scale parameter
  attr(f, 'scale') <- ifelse(any(names(formals(cov_fct)$param) == 'scale'), which(names(formals(cov_fct)$param) == 'scale') - 1, NA)
  # Include covariance function evaluated with additional arguments
  eval_cov_fct <- function(t, param) cov_fct(t, param, ...)
  attributes(eval_cov_fct) <- attributes(cov_fct)
  attr(f, 'cov_fct') <- eval_cov_fct
  attr(f, 'noise') <- noise
  # If the covariance has been made non-stationary, modify
  if (!is.null(ns)) {
    fixed <- ifelse(!is.null(ns$fixed), ns$fixed, max(knots))
    m_ns <- knots - length(fixed)

    attr(f, 'param') <- c(local_scale = rep(1, m_ns), param)
    ns_cov_fct <- function(t, param) {
      ns_param <- rep(1, length = knots)
      ns_param[-fixed] <- param[1:m_ns]
      prod(spline(seq(ns$range[1], ns$range[2], length = ns$knots), ns_param, xout = t)$y) *
        cov_fct(abs(diff(t)), param[-(1:m_ns)])
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
#' @seealso \code{link{predict_curve}}
#' @export

cov_rect <- function(t, t_p, cov_fct, param, ...) {
  if (!is.null(attr(cov_fct, 'discrete')) && attr(cov_fct, 'discrete')) stop('Rectangular covariance matrices are not defined for discrete covariances.')

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

#' Construct list of precision matrices
#'
#' Generates a list of inverse covariance (precision) matrices for the supplied list of observation points.
#' @param t list of observation points corresponding to functional samples.
#' @param cov_fct covariance function to use.
#' @param param parameters for covariance function.
#' @param amp_fct functional basis that amplitude variation should be expressed in. Default is \code{NULL}.
#' @export

fill_precision <- function(t, cov_fct, param, amp_fct = NULL) {
  m <- sapply(t, length)

  # Check if covariance function is specified
  if (!is.null(cov_fct)) {
    # Check if functional basis for the amplitude variation is supplied
    if (!is.null(amp_fct)) {
      df <- attr(amp_fct, 'df')
      SSinv <- chol2inv(chol(amp_cov(1:df, param)))

      # Fill precision matrices
      Sinv <- list()
      for (i in 1:n) {
        A <- amp_fct(t[[i]])
        Sinv[[i]] <- diag(1, m[i]) - A %*% solve(SSinv + t(A) %*% A, t(A))
      }
      return(Sinv)
    } else {
      # No functional basis specified, invert covariance matrices

      # Check if inverse method is given
      if (!is.null(attr(cov_fct, 'inv_cov_fct'))) {
        lapply(t, attr(cov_fct, 'inv_cov_fct'), param = param)
      } else {
        # No inverse method given, manually invert
        lapply(lapply(t, cov_fct, param = param), function(x)  chol2inv(chol(x)))
      }
    }
  } else {
    # No covariance specified, assuming iid. Gaussian noise
    lapply(m, Diagonal, x = 1)
  }
}
