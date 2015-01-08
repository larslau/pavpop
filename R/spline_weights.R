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
#' @examples
#' t <- seq(0, 1, length = 10)
#' plot(v(0.4, t, 0.5), type = 'b', pch = 19)

spline_weights <- function(y, t, w, tw, Ainv, kts) {
  n <- length(y)
  btime <- sapply(1:n, function(i) v(w[, i], t[[i]], tw))
  btime <- as.numeric(unlist(btime))

  basis <- bs(btime, knots = kts, Boundary.knots = c(0, 1))[]

  c <- as.numeric(solve(t(basis) %*% Ainv %*% basis) %*% t(basis) %*% Ainv %*% unlist(y))

  return(c)
}
