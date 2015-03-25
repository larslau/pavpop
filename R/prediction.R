#' Predict individual curves
#'
#' This function generates predicts individual curves from
#' @param t evaluation points.
#' @param tau scale parameter.
#' @param type type of covariance, either 'motion' or 'bridge'.
#' @keywords curve prediction
#' @export
#' @examples
#'
#'

predf <- function(t_p, t, y, amp_cov_fct, warp_cov_fct, res, deriv = FALSE) {
  S <- amp_cov_fct(t, res$amp_cov_par)
  Ainv <- chol2inv(chol(S))

  C <- warp_cov_fct(tw, res$warp_cov_par)
  Cinv <- solve(C)

  At_p <- amp_cov_fct(t, c(est_amp_par, smoothness), t_p = t_p)
  if (deriv) {
    At_p       <- (amp_cov_fct(t, c(est_amp_par, smoothness), t_p = t_p + 1e-5) - amp_cov_fct(t, c(est_amp_par, smoothness), t_p = t_p - 1e-5)) / age_scale(2e-5)

  }
  Ainv         <- chol2inv(chol(amp_cov_fct(t, c(est_amp_par, smoothness))))

  At_p_d       <- (amp_cov_rect(t_p + 1e-5, t, c(est_amp_par, smoothness)) - amp_cov_rect(t_p - 1e-5, t, c(est_amp_par, smoothness))) / age_scale(2e-5)
  Ainv         <- chol2inv(chol(amp_cov_fct(t, c(est_amp_par, smoothness))))

  w <- optim(par = array(0, dim = c(mw)), fn = posterior, gr = NULL, method = 'Nelder-Mead', t = t, y = y, tw = tw, c = res$c, Ainv = Ainv, Cinv = Cinv, kts = kts, intercept = TRUE, smooth_warp = TRUE, increasing = TRUE)$par
  homeomorphisms <- 'soft'
  if (homeomorphisms == 'soft') w <- make_homeo(w, tw)

  theta_pers   <- basis_fct(v(w, t_p, tw, smooth = TRUE)) %*% res$c
  theta_pers_d <- basis_d_fct(v(w, t_p, tw, smooth = TRUE)) %*% res$c[-1]
  r            <- y - basis_fct(v(w, t, tw, smooth = TRUE)) %*% res$c
  pred         <- theta_pers + At_p %*% Ainv %*% r
  velo_pred    <- theta_pers_d + At_p_d %*% Ainv %*% r

  return(list(w = w, pred = pred, velo_pred = velo_pred, r = r, theta_pers = theta_pers, theta_pers_d = theta_pers_d, Ainv = Ainv, At_p = At_p, At_p_d = At_p_d))
}
