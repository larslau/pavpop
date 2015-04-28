# TODO:
# Include smooth possibility
# Include possibility for scaling of results

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
#' @param warp_cov_par warp covariance parameters.
#' @param warp_cov_fct warp covariance matrix function.
#' @param w warping parameters.
#' @param tw anchor points for the warping parameters.
#' @param deriv logical. Should the derivative of the curve be predicted?
#' @param type type of covariance, either 'motion' or 'bridge'.
#' @keywords curve prediction
#' @export

predict_curve <- function(t_p, t, y, c, basis_fct, amp_cov_fct, amp_cov_par, w, tw, deriv = FALSE) {
  S <- amp_cov_fct(t, amp_cov_par)
  S_p <- cov_rect(t, t_p, attr(amp_cov_fct, 'cov_fct'), amp_cov_par)
  Ainv <- chol2inv(chol(S))

  basis <- basis_fct(v(w, t_p, tw, smooth = TRUE))
  # adjust for estimation of derivatives
  if (deriv) {
    S_p <- (cov_rect(t, t_p + 1e-5, attr(amp_cov_fct, 'cov_fct'), amp_cov_par) - cov_rect(t, t_p - 1e-5, attr(amp_cov_fct, 'cov_fct'), amp_cov_par)) / (2e-5)
    basis <- (basis_fct(v(w, t_p + 1e-5, tw, smooth = TRUE)) - basis_fct(v(w, t_p - 1e-5, tw, smooth = TRUE))) / (2e-5)
  }

  theta_pers   <- basis %*% c
  r            <- y - basis_fct(v(w, t, tw, smooth = TRUE)) %*% c
  pred         <- theta_pers + S_p %*% Ainv %*% r

  return(pred)
}
