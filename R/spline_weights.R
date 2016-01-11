#' Estimate spline weights
#'
#' Estimate spline weights for given warps and covariances. The method seamlessly
#' handles positivity constraints (specified in the basis function).
#' @param y list, matrix or data.frame containing \eqn{n} functional observations.
#' @param t (warped) time points corresponding to y in the same format.
#' @param Sinv list of precision matrices for amplitude variation. If \code{NULL}, iid.
#' Gaussian noise is assumed
#' @param basis_fct function for generating a basis.
#' @param weights weights for the individual observations. Mainly used for pavpop clustering.
#' @keywords spline basis
#' @seealso make_basis_fct
#' @export
#' @examples
#' # Evaluation points
#' t <- seq(0, 1, length = 100)
#'
#' # Simulate data
#' y <- t^2 * sin(8 * t) + t
#' plot(t, y, type = 'l', lwd = 2, lty = 2)
#'
#' # Add noise to data
#' y <- y + rnorm(length(y), sd = 0.1)
#' points(t, y, pch = 19, cex = 0.5)
#'
#' # Basis function knots
#' kts <- seq(0, 1, length = 12)[2:11]
#'
#' # Construct B-spline basis function
#' basis_fct <- make_basis_fct(kts = kts, control = list(boundary = c(0, 1)))
#'
#' # Fit B-spline to data assuming iid. noise
#' weights <- spline_weights(y, t, basis_fct = basis_fct)
#' lines(t, basis_fct(t) %*% weights, col = 'red', lwd = 2)

spline_weights <- function(y, t, Sinv = NULL, basis_fct, weights = NULL) {
  if (class(y) != 'list') {
    t <- as.matrix(t)
    y <- as.matrix(y)
    t <- lapply(1:ncol(t), function(i) t[, i])
    y <- lapply(1:ncol(y), function(i) y[, i])
  }

  n <- length(y)
  m <- sapply(y, length)

  nb <- attr(basis_fct, 'df')

  if (is.null(weights)) weights <- rep(1, n)

  Dmat <- matrix(0, nb, nb)
  dvec <- matrix(0, nb, 1)

  # If no precision matrix is supplied, use identity
  if (is.null(Sinv)) Sinv <- lapply(m, Matrix::Diagonal, x = 1)

  for (i in 1:n) {
    basis <- basis_fct(t[[i]])
    bSinv <- weights[i] * (t(basis) %*% Sinv[[i]])
    Dmat <- Dmat + bSinv %*% basis
    dvec <- dvec + bSinv %*% y[[i]]
  }

  if (attr(basis_fct, 'constraints') == 'positive') {
    qrDmat <- qr(Dmat)
    indices <- qrDmat$pivot[1:qrDmat$rank]

    c <- rep(0, ncol(basis))
    c[indices] <- quadprog::solve.QP(Dmat = Dmat[indices, indices],
                                     dvec = dvec[indices,],
                                     Amat = diag(nrow = length(indices)))$solution
  } else {
    # Faster to use qr?
    c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
  }

  return(c)
}
