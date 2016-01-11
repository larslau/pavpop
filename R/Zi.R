#' Jacobian of warped B-spline function
#'
#' Computes the Jacobian of a warped B-spline function.
#' @param t (warped) evaluation points.
#' @param dwarp Jacobian of the warping function for the given evaluation points.
#' @param basis_fct basis function.
#' @param c spline weights.
#' @export

Zi <- function(t, dwarp, basis_fct, c) {
  basis <- basis_fct(t, deriv = TRUE)
  dwarp <- dwarp * ((basis %*% c)[ , 1])
  return(dwarp)
}

Zis <- function(t, dwarp, basis_fct, c) {
  n <- length(t)
  m <- sapply(t, length)
  mw <- ncol(dwarp[[1]])
  Zis <- list()
  for (i in 1:n) {
    Zis[[i]] <- matrix(Zi(t[[i]], dwarp[[i]], basis_fct, c), m[i], mw)
  }
  return(Zis)
}
