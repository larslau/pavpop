#' Estimate spline weights
#'
#' Estimate spline weights for given warps and covariances
#' @param y list of \eqn{n} functional observations.
#' @param t list of time points corresponding to y.
#' @param w warping parameters.
#' @param tw anchor points for the warping parameters.
#' @param Ainv global precision matrix for amplitude variation.
#' @param kts B-spline knots.
#' @keywords warping
#' @export

spline_weights <- function(y, t, w, tw, Ainv, kts, intercept = FALSE) {
  n <- length(y)
  btime <- sapply(1:n, function(i) v(w[, i], t[[i]], tw))
  btime <- as.numeric(unlist(btime))

  basis <- bs(btime, knots = kts, Boundary.knots = c(0, 1), intercept = intercept)[]

  c <- as.numeric(solve(t(basis) %*% Ainv %*% basis) %*% t(basis) %*% Ainv %*% unlist(y))

  return(c)
}
