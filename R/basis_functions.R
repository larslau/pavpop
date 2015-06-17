#' Generate basis function
#'
#' Method for generating a 'basis function' function
#' @param kts a sequence of increasing points specifying the placement of the knots.
#' @param intercept logical. Should the basis include an intercept?
#' @param increasing logical. Should the basis be an I-spline?
#' @param order the order of the basis splines.
#' @param boundary boundary knots.
#' @keywords spline
#' @export


#TODO: EXPLICIT DERIVATIVES!!
make_basis_fct <- function(kts, intercept = FALSE, increasing = FALSE, order = 3, boundary = c(0, 1)) {
  epsilon <- 1e-5
  if (increasing) { # I-spline
    b <- function(t, deriv = FALSE) {
      basis <- ispline(t + deriv * epsilon, knots = kts, d = order)
      if (deriv) basis <- (basis - ispline(t - deriv * epsilon, knots = kts, d = order)) / (2 * epsilon)
      if (intercept) {
        basis <- cbind(!deriv, matrix(basis, ncol = length(kts)))
      }
      dims <- dim(basis)
      basis <- as.numeric(basis)
      dim(basis) <- dims
      return(basis)
    }
  } else { # B-spline
    #TODO: CAN BE EASILY HANDELED WITH splineDesign, includes derivs argument
    b <- function(t, deriv = FALSE) {
      basis <- bs(t + deriv * epsilon, knots = kts, degree = order, Boundary.knots = boundary, intercept = intercept)
      if (deriv) basis <- (basis - bs(t - deriv * epsilon, knots = kts, degree = order, Boundary.knots = boundary, intercept = intercept)) / (2 * epsilon)
      return(basis[])
    }
    
  }
  attr(b, 'df') <- length(kts) + order - 2 * increasing + intercept
  attr(b, 'intercept') <- intercept
  attr(b, 'increasing') <- increasing
  return(b)
}



#' Generate increasing spline basis
#'
#' Method for generating a 'basis function' function
#' @param x predictor variable.
#' @param knots the internal breakpoints of the spline.
#' @param d order of the spline functions.
#' @keywords spline
#' @note Method taken from the \code{SVMMaj} package.
#' @export

ispline <- function (x, knots, d) {
  if (is.null(knots) || any(is.na(knots)) || any(diff(knots) == 0) || length(knots) <= 2)
    return(x)
  m <- length(knots)
  n <- length(x)
  interval <- findInterval(x, knots, all.inside = TRUE)
  M <- sapply(sequence(m - 1), `==`, interval)
  for (i in 2:(d + 1)) {
    tik <- c(knots[-1], rep(knots[m], i - 2))
    ti <- c(rep(knots[1], i - 2), knots[-m])
    M <- M %*% diag(1/(tik - ti))
    Dx <- Dt <- array(0, dim = c(m + i - 3, m + i - 2))
    Dx[1L + 0L:(m + i - 4L) * (m + i - 2L)] <- -1
    Dx[1L:(m + i - 3L) * (m + i - 2L)] <- 1
    Dt[1L + 0L:(m + i - 4L) * (m + i - 2L)] <- tik
    Dt[1L:(m + i - 3L) * (m + i - 2L)] <- -ti
    M <- (M * x) %*% Dx + M %*% Dt
  }
  M <- M[, -1]
  S <- array(1, dim = rep(NCOL(M), 2))
  S[upper.tri(S)] <- 0
  I <- M %*% S
  return(I)
}
