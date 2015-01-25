# TODO: PYRAMID SCHEME! (2 types)
# TODO: AUTOMATICALLY INITIALIZE!
# TODO: MAKE HARD HOMEOMORPHIC CONSTRAINTS POSSIBLE
# TODO LATER: MAKE LARGE WARPS POSSIBLE (AND CHEAP)
# TODO LATER: ALLOW LARGE m (USING SPARSE CHOL)
# TODO LATER: ALLOW SHIFT WARPS


#' Estimate parameters and predict warps for curve data
#'
#' This function does likelihood estimation in the model \deqn{y_i(t)=\theta(v(t, w_i))+x_i(t)+\epsilon_i(t)} based on iterative local linearization of the model around predictions of the random warping parameters \eqn{w_i}.
#' @param y list of \eqn{n} functional observations. Missing values are allowed.
#' @param t list of time points corresponding to y. Should be scaled to have outer endpoints at 0 and 1.
#' @param amp_cov_par amplitude covariance parameters.
#' @param amp_cov_fct amplitude covariance matrix function.
#' @param warp_cov_par warp covariance parameters.
#' @param warp_cov_fct warp covariance matrix function.
#' @param kts anchor points for the B-spline basis used to model \eqn{\theta}.
#' @param tw anchor points for the warping parameters.
#' @param iter two-dimensional numeric consisting of number of outer and inner iterations.
#' @param use_warp_gradient logical. Should warp prediction use gradient based optimization?
#' @param homeomorphisms should warps be constrained to be homeomorphisms? Options are: \code{'no'}, \code{'soft'} or \code{'hard'}. 'soft' will project the prediction onto the space of homeomorphisms after each prediction. 'hard' will do the optimiziation in the constrained space (not implemented yet!).
#' @param like_optim_control list of control options for likelihood optimization. Parameters are given as \code{c(amp_cov_par, warp_cov_par)} and options include lower, upper, method, ndev (see \code{\link[stats::optim]{optim}}).
#' @keywords likelihood estimation
#' @export


estimate_growth <- function(y, t, amp_cov_par, amp_cov_fct, warp_cov_par, warp_cov_fct, kts, tw, intercept = FALSE, iter = c(5, 5), use_warp_gradient = FALSE, homeomorphisms = 'no', like_optim_control = list()) {
  nouter <- iter[1] + 1
  ninner <- iter[2]

  # Set size parameters
  n <- length(y)
  m <- sapply(y, length)
  mw <- length(tw)
  n_par_amp <- length(amp_cov_par)
  n_par_warp <- length(warp_cov_par)

  # Check for same data structures of y and t
  if (length(t) != n) stop("y and t must have same length.")
  if (!all(sapply(t, length) == m)) stop("Observations in y and t must have same length.")

  # Remove missing values
  for (i in 1:n) {
    missing_indices <- is.na(y[[i]])
    y[[i]] <- y[[i]][!missing_indices]
    t[[i]] <- t[[i]][!missing_indices]
  }
  # Update m with cleaned data
  m <- sapply(y, length)

  # Initialize warp parameters
  w <- array(0, dim = c(mw, n))

  # Build amplitude covariances and inverse covariances
  S <- Ainv <- list()
  for (i in 1:n) {
    S[[i]] <- amp_cov_fct(t[[i]], amp_cov_par)
    Ainv[[i]] <- chol2inv(chol(S[[i]]))
  }

  # Build warp covariance and inverse
  C <- warp_cov_fct(tw, warp_cov_par)
  Cinv <- solve(C)

  # Estimate spline weights
  c <- spline_weights(y, t, w, tw, bdiag(Ainv), kts, intercept = intercept, smooth_warp = TRUE, increasing = TRUE)

  for (iouter in 1:nouter) {
    # Outer loop
    if (iouter != nouter) cat(iouter, ':\t')

    for (iinner in 1:ninner) {
      # Inner loop
      if (iouter != nouter) cat(iinner, '\t')

      # Predict warping parameters for all functional samples
      if (homeomorphisms == 'hard') {
        #TODO: constrainOptim
        stop("Hard homeomorphic constrained optimization for warps is not implemented.")
      } else {
        for (i in 1:n) {
          gr <- NULL
          warp_optim_method <- 'Nelder-Mead'
          if (use_warp_gradient) {
            stop("Exact gradient based optimization is not available at the moment.")
            #TODO: construct dwarp explicitly using splinefun
#             gr <- function(w, t, y, tw, c, Ainv, Cinv, kts, intercept) posterior_grad(w, dwarp[[i]], t, y, tw, c, Ainv, Cinv, kts, intercept, )
            warp_optim_method <- 'BFGS'
          }
          w[, i] <- optim(par = w[,i], fn = posterior, gr = gr, method = warp_optim_method, t = t[[i]], y = y[[i]], tw = tw, c = c, Ainv = Ainv[[i]], Cinv = Cinv, kts = kts, intercept = intercept, smooth_warp = TRUE, increasing = TRUE)$par
          if (homeomorphisms == 'soft') w[, i] <- make_homeo(w[, i], tw)
        }
      }

      # Update spline weights
      c <- spline_weights(y, t, w, tw, bdiag(Ainv), kts, intercept = intercept, smooth_warp = TRUE, increasing = TRUE)
    }

    # Likelihood estimation of parameters (outer loop)

    # Construct residual vector for given warp prediction
    Zis <- list()
    r <- y
    for (i in 1:n) {
      twarped <- v(w[, i], t[[i]], tw, smooth = TRUE)
      dwarp <- matrix(0, m[i], mw)
      for (k in 1:mw) {
        h_tmp <- rep(0, mw)
        h_tmp[k] <- 1e-4
        dwarp[,k] <- (v(w[, i] + h_tmp, t[[i]], tw, smooth = TRUE) - v(w[, i] - h_tmp, t[[i]], tw, smooth = TRUE)) / (2e-4)
      }
      Zis[[i]] <- Zi(twarped, dwarp, c, kts, increasing = TRUE)
      basis <- t(Ispline(twarped, 3, kts))
      if (intercept) basis <- cbind(1, basis)

      r[[i]] <- as.numeric(r[[i]] - basis %*% c + Zis[[i]] %*% w[,i])
    }

    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      # Estimate parameters using locally linearized likelihood
      lower  <- if (is.null(like_optim_control$lower)) rep(0, n_par_amp + n_par_warp) else like_optim_control$lower
      upper  <- if (is.null(like_optim_control$upper)) rep(Inf, n_par_amp + n_par_warp) else like_optim_control$upper
      method <- if (is.null(like_optim_control$method)) "L-BFGS-B" else like_optim_control$method
      ndeps <- if (is.null(like_optim_control$ndeps)) rep(1e-3, n_par_amp + n_par_warp) else like_optim_control$ndeps

      param <- optim(c(amp_cov_par, warp_cov_par), like, n_par = c(n_par_amp, n_par_warp), r = r, Zis = Zis, amp_cov_fct = amp_cov_fct, warp_cov_fct = warp_cov_fct, t = t, tw = tw, method = method, lower = lower, upper = upper, control = list(ndeps = ndeps, maxit = 10))$par

      amp_cov_par <- param[1:n_par_amp]
      warp_cov_par <- param[(n_par_amp + 1):length(param)]

      # Update covariances
      for (i in 1:n) {
        S[[i]] <- amp_cov_fct(t[[i]], amp_cov_par)
        Ainv[[i]] <- chol2inv(chol(S[[i]]))
      }

      C <- warp_cov_fct(tw, warp_cov_par)
      Cinv <- solve(C)

    } else {
      # Estimate of sigma if final iteration is reached
      sigma <- sqrt(sigmasq(c(amp_cov_par, warp_cov_par), c(n_par_amp, n_par_warp), r, Zis, amp_cov_fct, warp_cov_fct, t, tw))
    }

    if (iouter != nouter) cat(':\t', param, '\n')
  }
  return(list(c = c, w = w, amp_cov_par = amp_cov_par, warp_cov_par = warp_cov_par, sigma = sigma))
}
