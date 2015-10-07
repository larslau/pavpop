# TODO: Include possibility for scaling of results
# TODO: Automatic handeling of derivatives
# TODO: Make faster for non-stationary Matern covariance (smarter rectangular covariance function!)

#' Predict individual curves
#'
#' This function predicts individual curves from estimated parameters and warps and returns verbose information about the prediction
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
#' @return list with predicted curve, warp parameters, and additional information
#' @keywords curve prediction
#' @export

predict_curve_verbose <- function(t_p, t, y, c, basis_fct, warp_fct, amp_cov, amp_cov_par, w = NULL, warp_cov = NULL, warp_cov_par = NULL, deriv = FALSE) {
  S <- amp_cov(t, amp_cov_par)
  S_p <- cov_rect(t, t_p, attr(amp_cov, 'cov_fct'), amp_cov_par)
  Ainv <- chol2inv(chol(S))

  if (is.null(w)) {
    if (is.null(warp_cov) | is.null(warp_cov_par)) {
      stop('Warp covariance function and warp covariance parameters must be supplied if no warping parameters are given.')
    }
    tw <- attr(warp_fct, 'tw')
    mw <- attr(warp_fct, 'mw')
    Cinv <- solve(warp_cov(tw, warp_cov_par))
    w <- optim(par = rep(0, mw), fn = posterior, warp_fct = warp_fct, t = t, y = y, c = c, Sinv = Ainv, Cinv = Cinv, basis_fct = basis_fct)$par
    w <- make_homeo(w, tw)
  }

  basis <- basis_fct(warp_fct(w, t_p))

  # adjust for estimation of derivatives
  if (deriv) {
    S_p <- (cov_rect(t, t_p + 1e-5, attr(amp_cov, 'cov_fct'), amp_cov_par) - cov_rect(t, t_p - 1e-5, attr(amp_cov, 'cov_fct'), amp_cov_par)) / (2e-5)
    basis <- (basis_fct(warp_fct(w, t_p + 1e-5)) - basis_fct(warp_fct(w, t_p - 1e-5))) / (2e-5)
  }

  r <- y - basis_fct(warp_fct(w, t)) %*% c
  theta_pers   <- basis %*% c
  pred         <- theta_pers + S_p %*% Ainv %*% r
  if (is.null(warp_cov_par)) {
    post <- NULL
  } else {
    tw <- attr(warp_fct, 'tw')
    mw <- attr(warp_fct, 'mw')
    Cinv <- solve(warp_cov(tw, warp_cov_par))
    post <- posterior(w = w, t = t, y = y, c = c, Cinv = Cinv, basis_fct = basis_fct, warp_fct = warp_fct, Sinv = Ainv)
  }
  return(list(pred = pred, w = w, r = r, post = post))
}

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
#' @return predicted curve
#' @keywords curve prediction
#' @export

predict_curve <- function(t_p, t, y, c, basis_fct, warp_fct, amp_cov, amp_cov_par, w = NULL, warp_cov = NULL, warp_cov_par = NULL, deriv = FALSE) {
  return(predict_curve_verbose(t_p, t, y, c, basis_fct, warp_fct, amp_cov, amp_cov_par, w, warp_cov, warp_cov_par, deriv)$pred)
}
