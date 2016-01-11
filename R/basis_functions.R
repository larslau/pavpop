# TODO: LOOK INTO B-SPLINE BEHAVIOR AT BOUNDARY. PROBLEMATIC OR OKAY?

# TODO: Wavelet basis
# TODO: Fourier basis
# TODO: I-splines can be computed more efficiently using (33) in de Boor.

#' Generate basis function
#'
#' Method for generating a 'basis function' function, that outputs a functional basis.
#'
#' Basis types \code{'Fourier'} and \code{'wavelet'}
#' will become available soon. If needed please contact the developer.
#'
#' The control argument takes a list with the following entries
#' \code{order} and \code{constraints} and \code{sparse}
#' \describe{
#'   \item{\code{boundary}}{boundary knots for the basis spline.}
#'   \item{\code{order}}{order of the spline, if \code{NULL}, B-splines have order
#'   4 (cubic spline) and I-splines (\code{type = 'increasing'}) have order 3.}
#'   \item{\code{constraints}}{positivity constraints, if set to \code{'positive'},
#'   only positive weights are allowed}
#'   \item{\code{sparse}}{logical. Should sparse matrices be used?}
#' }
#'
#'
#' @param kts a sequence of increasing points specifying the placement of the knots.
#' @param df degrees of freedom of the spline basis. Knots are chosen equidistantly.
#' @param type the type of basis function you want. Currently supported choices are \code{'B-spline'},
#' \code{'increasing'} and \code{'intercept'}. See details for more information.
#' @param intercept logical. Should the basis include an intercept?
#' @param control list of control parameters. Most importantly is \code{boundary} which
#' contains boundary points for a B-spline basis. See details for more options.
#' @keywords spline
#' @export
#' @examples
#' # Basis function knots
#' kts <- seq(0, 1, length = 12)[2:11]
#'
#' # Construct B-spline basis function
#' basis_fct <- make_basis_fct(kts = kts, control = list(boundary = c(0, 1)))
#'
#' # Evaluation points
#' t <- seq(0, 1, length = 100)
#' A <- basis_fct(t)
#' plot(t, t, type = 'n', ylim = range(A))
#' for (i in 1:ncol(A)) lines(t, A[, i], col = rainbow(ncol(A))[i])
#'
#' # Evaluate derivatives
#' Ad <- basis_fct(t, TRUE)
#' plot(t, t, type = 'n', ylim = range(Ad))
#' for (i in 1:ncol(A)) lines(t, Ad[, i], col = rainbow(ncol(Ad))[i])
#'
#' # Construct I-spline
#' # Knots should contain the left and right endpoints
#' kts_inc <- seq(0, 1, length = 10)
#' basis_fct_inc <- make_basis_fct(kts = kts_inc, type = 'increasing')
#' A_inc <- basis_fct_inc(t)
#' plot(t, t, type = 'n', ylim = range(A_inc))
#' for (i in 1:ncol(A_inc)) lines(t, A_inc[, i], col = rainbow(ncol(A))[i])
#'
#' # Evaluate derivatives
#' Ad_inc <- basis_fct_inc(t, deriv = TRUE)
#' plot(t, t, type = 'n', ylim = range(Ad_inc))
#' for (i in 1:ncol(Ad_inc)) lines(t, Ad_inc[, i], col = rainbow(ncol(Ad))[i])
#'
#' # Simulate data
#' y <- t^2 * sin(8 * t) + t
#' plot(t, y, type = 'l', lwd = 2, lty = 2)
#'
#' # Add noise to data
#' y <- y + rnorm(length(y), sd = 0.1)
#' points(t, y, pch = 19, cex = 0.5)
#'
#' # Fit B-spline to data assuming iid. noise
#' weights <- spline_weights(y, t, basis_fct = basis_fct)
#' lines(t, A %*% weights, col = 'red', lwd = 2)
#'
#' # Fit increasing spline
#' pos_weights <- spline_weights(y, t, basis_fct = basis_fct_inc)
#' lines(t, A_inc %*% pos_weights, col = 'blue', lwd = 2)


make_basis_fct <- function(kts = NULL, df = NULL, type = 'B-spline', intercept = FALSE, control = list()) {
  # Match type
  types <- c('B-spline', 'increasing', 'intercept', 'Fourier', 'wavelet')
  type <- types[pmatch(type, types)]

  if (is.na(type)) stop('Invalid type of basis.')

  # Check that kts or df is supplied if type is different from 'intercept'
  if (type != 'intercept' & is.null(kts) & is.null(df)) stop('You must supply either list of knots or degrees of freedom.')


  # Control boundary
  if (!is.null(control$boundary)) {
    boundary = control$boundary
  } else {
    # If kts is supplied, boundary is given by extending the knot list
    # one "step-length" under the assumption of equidistant knots
    if (!is.null(kts)) {
      boundary <- range(kts) + diff(range(kts)) * c(-1 / (length(kts) - 1), 1 / (length(kts) - 1))
    } else if (type != 'intercept') {
      # No boundary supplied, no knots
      warning('No boundary knots or evaluation knots supplied, assuming observation points are in (0, 1)')
      boundary <- c(0, 1)
    }
  }

  # Control constaints
  if (is.null(control$constraints) & type != 'increasing') {
    constraints <- 'none'
  } else {
    constraints <- 'positive'
  }

  #
  # B-spline basis
  #
  if (type == 'B-spline') {
    # Extract order og B-spline basis
    if (is.null(control$order)) {
      order <- 4
    } else {
      order <- control$order
    }
    if (is.null(kts)){
      kts <- seq(boundary[1], boundary[2], length = df - order + (3L - intercept))
      kts <- head(tail(kts, -1), -1)
    }
    Aknots <- sort(c(rep(boundary, order), kts))

    # Should sparse matrices be used?
    if (is.null(control$sparse)) {
      sparse <- FALSE
    } else {
      sparse <- control$sparse
    }

    # Basis function to return
    b <- function(t, deriv = FALSE) {
      basis <- splines::splineDesign(Aknots, t, ord = order, derivs = deriv, outer.ok = TRUE, sparse = sparse)
      if (!intercept) basis <- basis[, -1]
      return(basis[])
    }
    attr(b, 'boundary') <- boundary
  }

  #
  # Increasing spline basis
  #
  if (type == 'increasing') {
    # Extract order og I-spline basis
    if (is.null(control$order)) {
      order <- 3
    } else {
      order <- control$order
    }

    # Set boundary and B-spline knots (for derivative)
    boundary <- range(kts)
    Aknots <- sort(c(rep(boundary, order - 1), kts))

    b <- function(t, deriv = FALSE) {
      if (!deriv) {
      basis <- ispline(t, knots = kts, d = order)
      } else {
        basis <- t(c(order / 1:(order - 1), rep(1, length(kts) - order), order / (order - 1):1) * t(splines::splineDesign(Aknots, t, ord = order, derivs = 0, outer.ok = TRUE) * (length(kts) - 1) / diff(range(kts))))
      }
      if (intercept) {
        basis <- cbind(!deriv, matrix(basis, nrow = length(t)))
      }
      dims <- dim(basis)
      basis <- as.numeric(basis)
      dim(basis) <- dims
      return(basis)
    }
  }

  #
  # Only intercept
  #
  if (type == 'intercept') {
    b <- function(t, deriv = FALSE) {
      t[] <- ifelse(deriv, 0, 1)
      dim(t) <- c(length(t), 1)
      return(t)
    }
    attr(b, 'df') <- 1
    attr(b, 'intercept') <- TRUE

    return(b)
  }

  attr(b, 'df') <- length(kts) + order - 1 - 1 * (type == 'increasing') + intercept
  attr(b, 'intercept') <- intercept
  attr(b, 'constraints') <- constraints
  return(b)
}


#' Generate increasing spline basis
#'
#' Method for generating a 'basis function' function
#' @param t evaluation points.
#' @param knots the internal breakpoints of the spline.
#' @param d order of the spline functions.
#' @keywords spline
#' @note Method is a corrected version of the increasing spline basis code in the \code{SVMMaj} package.
#' @export

ispline <- function (t, knots, d) {
  if (is.null(knots) || any(is.na(knots)) || any(diff(knots) == 0) || length(knots) <= 2)
    return(t)
  m <- length(knots)
  n <- length(t)
  interval <- findInterval(t, knots, all.inside = TRUE)
  M <- sapply(sequence(m - 1), `==`, interval)
  for (i in 2:(d + 1)) {
    tik <- c(knots[-1], rep(knots[m], i - 2))
    ti <- c(rep(knots[1], i - 2), knots[-m])
    M <- M %*% diag(1 / (tik - ti))
    Dx <- Dt <- array(0, dim = c(m + i - 3, m + i - 2))
    Dx[1L + 0L:(m + i - 4L) * (m + i - 2L)] <- -1
    Dx[1L:(m + i - 3L) * (m + i - 2L)] <- 1
    Dt[1L + 0L:(m + i - 4L) * (m + i - 2L)] <- tik
    Dt[1L:(m + i - 3L) * (m + i - 2L)] <- -ti
    M <- (M * t) %*% Dx + M %*% Dt
  }
  M <- M[, -1, drop = FALSE]
  S <- array(1, dim = rep(NCOL(M), 2))
  S[upper.tri(S)] <- 0
  I <- M %*% S
  I[t > max(knots), ] <- 1
  return(I)
}
