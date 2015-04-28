#' Jacobian of warped B-spline function
#'
#' Computes the Jacobian of a warped B-spline function.
#' @param t (warped) evaluation points.
#' @param dwarp Jacobian of the warping function for the given evaluation points.
#' @param c spline weights.
#' @param basis_fct basis function.
#' @export

Zi <- function(t, dwarp, c, basis_fct) {
  basis <- basis_fct(t, deriv = TRUE)
  dwarp <- dwarp * ((basis %*% c)[ , 1])
  return(dwarp)
}
