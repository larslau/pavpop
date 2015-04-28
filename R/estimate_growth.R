#' Estimate parameters and predict warps for growth curve data
#'
#' This function does likelihood estimation in the model \deqn{y_i(t)=\theta(v(t, w_i))+x_i(t)+\epsilon_i(t)} based on iterative local linearization of the model around predictions of the random warping parameters \eqn{w_i}.
#' @param y list of \eqn{n} functional observations. Missing values are allowed.
#' @param t list of time points corresponding to y. Should be scaled to have outer endpoints at 0 and 1.
#' @param kts a sequence of increasing points specifying the placement of the knots for the basis.
#' @param tw anchor points for the warping parameters.
#' @param iter two-dimensional numeric consisting of number of outer and inner iterations.
#' @keywords likelihood estimation
#' @seealso \code{\link{estimate_generic}} for a more flexible framework.
#' @export


estimate_growth <- function(y, t, kts, tw, iter = c(5, 5)) {
  nouter <- iter[1] + 1
  ninner <- iter[2]

  basis_fct <- make_basis_fct(kts = kts, intercept = TRUE, increasing = TRUE,
                              order = 3, boundary = c(0, 1))

  # Set up covariance functions
  warp_cov_par <- c(tau = 1)
  warp_cov_fct <- make_cov_fct(Brownian, noise = FALSE, type = 'bridge')

  amp_cov_par <- c(scale = 50, range = 1, smoothness = 1.5)
  amp_cov_fct <- make_cov_fct(Matern, noise = TRUE)

  # Bounds of parameters
  lower <- c(1e-2, 1e-2, 0.5001, 1e-2)
  upper <- c(200, 1, 3, 1)

  res <- estimate_generic(y, t, basis_fct, amp_cov_par, amp_cov_fct,
                          warp_cov_par, warp_cov_fct, tw, iter = iter,
                          smooth_warp = TRUE, homeomorphism = 'soft',
                          like_optim_control = list(lower = lower,
                                                    upper = upper))

  return(res)
}
