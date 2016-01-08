# TODO: HANDLE FORCE INCREASING (PREDICTION)
# TODO: HANDLE WARPED_AMP (PREDICTION)
# TODO: HANDLE HOMEOMORPHIC CONSTRAINTS!
# TODO: Make faster for non-stationary Matern covariance (smarter rectangular covariance function!)

#' Predict individual curves
#'
#' This function predicts individual curves given the estimated parameters and predicted warps.
#' @param t_p evaluation points for prediction.
#' @param t observation points.
#' @param y observation values.
#' @param basis_fct basis function to describe the mean function.
#' @param c weights for basis functions.
#' @param warp_fct warping function.
#' @param w warping parameters (optional). If parameters are not supplied, they are
#' predicted from the parameter estimates in the model.
#' @param amp_cov amplitude covariance matrix function.
#' @param amp_cov_par amplitude covariance parameters.
#' @param amp_fct functional basis that amplitude variation should be expressed in.
#' If used, it is assumed that iid. Gaussian noise is present in the data.
#' @param warp_cov warp covariance matrix function. If \code{w} is supplied,
#' the argument is ignored.
#' @param warp_cov_par warp covariance parameters. If \code{w} is supplied,
#' the argument is ignored.
#' @param deriv logical. Should the derivative of the curve be predicted?
#' @param verbose_output logical. Should output also include predicted warps, residuals
#' and the value of the posterior w | y.
#' @return predicted curve
#' @keywords curve prediction
#' @export

predict_curve <- function(t_p, t, y, basis_fct, c, warp_fct, w = NULL, amp_cov = NULL, amp_cov_par = NULL, amp_fct = NULL, warp_cov = NULL, warp_cov_par = NULL, deriv = FALSE, verbose_output = FALSE) {
  m <- length(t)
  m_p <- length(t_p)
  # Construct amplitude covariance/precision corresponding to prediction/observation points
  if (is.null(amp_cov)) {
    Sinv <- diag(m)
    S_p <- matrix(0, m_p, m)
  } else {
    if (!is.null(amp_fct)) {
      df <- attr(amp_fct, 'df')
      S <- amp_cov(1:df, amp_cov_par)
      SSinv <- chol2inv(chol(S))
      A <- amp_fct(t)
      A_p <- amp_fct(t_p, deriv = deriv)
      Sinv <- diag(1, m) - A %*% solve(SSinv + t(A) %*% A, t(A))
      S_p <- A_p %*% S %*% t(A)
    } else {
      if (!is.null(attr(amp_cov, 'inv_cov_fct'))) {
        Sinv <- attr(amp_cov, 'inv_cov_fct')(t, amp_cov_par)
      } else {
        # No inverse method given, manually invert
        Sinv <- chol2inv(chol(amp_cov(t, amp_cov_par)))
      }
      if (deriv) {
        S_p <- (cov_rect(t, t_p + 1e-5, attr(amp_cov, 'cov_fct'), amp_cov_par)
                - cov_rect(t, t_p - 1e-5, attr(amp_cov, 'cov_fct'), amp_cov_par)) / (2e-5)
      } else {
        S_p <- cov_rect(t, t_p, attr(amp_cov, 'cov_fct'), amp_cov_par)
      }
    }
  }

  # Construct warp precision
  if (!is.null(warp_cov)) {
    mw <- attr(warp_fct, 'mw')
    tw <- attr(warp_fct, 'tw')[2:(mw + 1)]
    Cinv <- solve(warp_cov(tw, warp_cov_par))
  }

  # Check if predicted warp parameters are given
  if (is.null(w)) {
    if (is.null(warp_cov) | is.null(warp_cov_par)) {
      stop('Warp covariance function and warp covariance parameters must be supplied if no warping parameters are given.')
    }

    w <- optim(par = rep(0, mw), fn = posterior, method = 'CG', warp_fct = warp_fct, t = t, y = y, c = c, Sinv = Sinv, Cinv = Cinv, basis_fct = basis_fct)$par

    # w <- make_homeo(w, tw)
  }

  # Construct basis
  basis <- basis_fct(warp_fct(w, t_p), deriv)
  if (deriv) basis <- basis * (warp_fct(w, t_p + 1e-5) - warp_fct(w, t_p - 1e-5)) / 2e-5

  r <- y - basis_fct(warp_fct(w, t)) %*% c
  theta_pers   <- basis %*% c
  pred         <- theta_pers + S_p %*% Sinv %*% r

#   if (force_increasing) {
#     ispline <- make_basis_fct(kts = kts, intercept = TRUE, control = list(type = 'increasing'))
#     c <- spline_weights(list(pred), list(t_p), warp_fct, array(0,dim=c(mw,1)), list(diag(length(t_p))), ispline)
#     pred <- ispline(t_p) %*% c
#   }

  if (is.null(warp_cov_par)) {
    post <- NULL
  } else {
    post <- posterior(w = w, t = t, y = y, c = c, Cinv = Cinv, basis_fct = basis_fct, warp_fct = warp_fct, Sinv = Sinv)
  }
  result <- list(pred = pred, w = w, r = r, post = post)
  if (!verbose_output) result <- result$pred

  return(result)
}
