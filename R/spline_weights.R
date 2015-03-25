#' Estimate spline weights
#'
#' Estimate spline weights for given warps and covariances
#' @param y list of \eqn{n} functional observations.
#' @param t list of time points corresponding to y.
#' @param w warping parameters.
#' @param tw anchor points for the warping parameters.
#' @param Ainv global precision matrix for amplitude variation.
#' @param kts B-spline knots.
#' @param smooth_warp logical. Should a smooth warping function be used?
#' @keywords warping
#' @export

spline_weights <- function(y, t, w, tw, Ainv, kts, intercept = FALSE, smooth_warp = FALSE, increasing = FALSE) {
  n <- length(y)
  btime <- sapply(1:n, function(i) v(w[, i], t[[i]], tw, smooth = smooth_warp))
  btime <- as.numeric(unlist(btime))

  if (!increasing) {
    basis <- bs(btime, knots = kts, Boundary.knots = c(0, 1), intercept = intercept)[]
    #TODO: Check! Can it be done smarter?!
    c <- as.numeric(MASS::ginv(as.matrix(t(basis) %*% Ainv %*% basis)) %*% t(basis) %*% Ainv %*% unlist(y))
  } else {
    basis <- t(Ispline(btime, 3, knots = kts))
    indices <- 1:(ncol(basis) + intercept)
    for (i in 1:ncol(basis)) {
      if (length(unique(basis[, i])) <= 3) { # Numerical precision hack. Should be == 1 in an ideal world
        indices <- indices[indices != i + intercept]
      }
    }
    if (intercept) basis <- cbind(1, basis)
    c <- rep(0, ncol(basis))
    basis <- basis[,indices]
    c[indices] <- solve.QP(Dmat = t(basis) %*% Ainv %*% basis,
                  dvec = Matrix::t(t(unlist(y)) %*% Ainv %*% basis),
                  Amat = diag(nrow = ncol(basis)))$solution
  }


  return(c)
}
