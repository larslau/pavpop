#' Estimate spline weights
#'
#' Estimate spline weights for given warps and covariances
#' @param y list of \eqn{n} functional observations.
#' @param t list of time points corresponding to y.
#' @param warp_fct warping function.
#' @param w warping parameters.
#' @param Sinv global precision matrix for amplitude variation.
#' @param basis_fct function for generating a basis.
#' @param weights weights for the individual observations
#' @keywords spline basis
#' @export

# spline_weights <- function(y, t, w, tw, Ainv, basis_fct, weights = NULL, smooth_warp = FALSE) {
#   n <- length(y)
#   m <- sapply(y, length)
#
#   btime <- sapply(1:n, function(i) v(w[, i], t[[i]], tw, smooth = smooth_warp))
#   btime <- as.numeric(unlist(btime))
#
#   basis <- basis_fct(btime)
#   attr(basis, 'class') <- 'matrix'
#   if (!is.null(weights)) Ainv <- as.numeric(unlist(sapply(1:n, function(x) rep(weights[x], each = m[x])))) * Ainv
#   #TODO: OPTIMIZE, can be optimized by numerical procedures?
#   Dmat <- t(basis) %*% Ainv %*% basis
#   dvec <- t(basis) %*% Ainv %*% unlist(y)
#
#   if (attr(basis_fct, 'increasing')) {
#     intercept <- attr(basis_fct, 'intercept')
#     indices <- 1:ncol(basis)
#     for (i in 1:ncol(basis)) {
#       #TODO: DO A PROPER FIX!
#       if (length(unique(basis[, i])) <= 3) { # Numerical precision hack. Should be == 1 in an ideal world
#         indices <- indices[indices != i + intercept]
#       }
#     }
#     c <- rep(0, ncol(basis))
#     basis <- basis[, indices]
#     c[indices] <- solve.QP(Dmat = Dmat[indices, indices],
#                            dvec = dvec[indices,],
#                            Amat = diag(nrow = ncol(basis)))$solution
#   } else {
#     c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
#   }
#
#   return(c)
# }

spline_weights <- function(y, t, warp_fct, w, Sinv, basis_fct, weights = NULL) {
  n <- length(y)
  m <- sapply(y, length)
  nb <- attr(basis_fct, 'df')

  if (is.null(weights)) weights <- rep(1, n)

  Dmat <- matrix(0, nb, nb)
  dvec <- matrix(0, nb, 1)

  for (i in 1:n) {
    basis <- basis_fct(warp_fct(w[, i], t[[i]]))
    bSinv <- weights[i] * (t(basis) %*% Sinv[[i]])
    Dmat <- Dmat + bSinv %*% basis
    dvec <- dvec + bSinv %*% y[[i]]
  }

  if (attr(basis_fct, 'increasing')) {
    intercept <- attr(basis_fct, 'intercept')
    indices <- 1:nb
    qrDmat <- qr(Dmat)
    indices <- qrDmat$pivot[1:qrDmat$rank]

    c <- rep(0, ncol(basis))
    c[indices] <- solve.QP(Dmat = Dmat[indices, indices],
                           dvec = dvec[indices,],
                           Amat = diag(nrow = length(indices)))$solution
  } else {
    c <- as.numeric(MASS::ginv(as.matrix(Dmat)) %*% dvec)
  }

  return(c)
}
