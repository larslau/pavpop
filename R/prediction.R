# TODO: Include possibility for scaling of results
# TODO: Automatic handeling of derivatives

#' Predict individual curves
#'
#' This function predicts individual curves from estimated parameters and warps.
#' @param t_p evaluation points for prediction.
#' @param t observation points.
#' @param y observation values.
#' @param c weights for basis functions.
#' @param basis_fct basis function to describe the mean function.
#' @param amp_cov_par amplitude covariance parameters.
#' @param amp_cov_fct amplitude covariance matrix function.
#' @param w warping parameters.
#' @param tw anchor points for the warping parameters.
#' @param deriv logical. Should the derivative of the curve be predicted?
#' @keywords curve prediction
#' @export

predict_curve <- function(t_p, t, y, c, basis_fct, warp_fct, amp_cov, amp_cov_par, w, deriv = FALSE) {
  S <- amp_cov(t, amp_cov_par)
  S_p <- cov_rect(t, t_p, attr(amp_cov, 'cov_fct'), amp_cov_par)
  Ainv <- chol2inv(chol(S))

  basis <- basis_fct(warp_fct(w, t_p))

  # adjust for estimation of derivatives
  if (deriv) {
    S_p <- (cov_rect(t, t_p + 1e-5, attr(amp_cov, 'cov_fct'), amp_cov_par) - cov_rect(t, t_p - 1e-5, attr(amp_cov, 'cov_fct'), amp_cov_par)) / (2e-5)
    basis <- (basis_fct(warp_fct(w, t_p + 1e-5)) - basis_fct(warp_fct(w, t_p - 1e-5))) / (2e-5)
  }

  r <- y - basis_fct(warp_fct(w, t)) %*% c
  theta_pers   <- basis %*% c
  pred         <- theta_pers + S_p %*% Ainv %*% r

  return(pred)
}
