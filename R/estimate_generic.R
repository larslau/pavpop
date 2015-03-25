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
#' @param kts anchor points for the spline basis used to model \eqn{\theta}.
#' @param tw anchor points for the warping parameters.
#' @param intercept should the basis include an intercept (default is \code{FALSE})
#' @param iter two-dimensional numeric consisting of number of outer and inner iterations.
#' @param use_warp_gradient logical. Should warp prediction use the exact gradient for based optimization?
#' @param homeomorphisms should warps be constrained to be homeomorphisms? Options are: \code{'no'}, \code{'soft'} or \code{'hard'}. 'soft' will project the prediction onto the space of homeomorphisms after each prediction. 'hard' will do the optimiziation in the constrained space (not implemented yet!).
#' @param like_optim_control list of control options for likelihood optimization. Parameters are given as \code{c(amp_cov_par, warp_cov_par)} and options include lower, upper, method, ndev (see \code{\link[stats::optim]{optim}}).
#' @keywords likelihood estimation
#' @export
#' @examples
#' # Load female growth data from the Berkeley growth study
#' t_orig <- fda::growth$age
#' y <- fda::growth$hgtf
#'
#' m <- nrow(y)
#' n <- ncol(y)
#'
#' # Construct velocities from finite differences
#' y <- (y[c(2:m, m), ] - y) / (t_orig[c(2:m, m)] - t_orig)
#' y[is.na(y)] <- 0
#' theta <- rowMeans(y)
#'
#' t_range <- c(0, 19)
#' t <- replicate(n, t_orig / t_range[2], simplify = FALSE)
#' y <- lapply(1:n, function (x) y[, x])
#'
#' # Set up covariance functions
#' warp_cov_par <- c(tau = 1)
#' warp_cov_fct <- function(t, param) Brownian_cov(t, param, type = 'bridge')
#'
#' amp_cov_par <- c(scale = 1, range = 1)
#' amp_cov_fct <- function(t, param) Matern_cov(t, c(param, 2), noise = TRUE)
#'
#' # Set up parametrization
#' tw <- seq(0, 1, length = 7)[2:6]
#' kts <- seq(0, 1, length = 12)[2:11]
#'
#' # Estimate in the model
#' res <- estimate_generic(y, t, amp_cov_par, amp_cov_fct, warp_cov_par, warp_cov_fct, kts, tw, iter = c(3, 3), homeomorphism = 'soft', like_optim_control = list(lower = rep(1e-2, 3), upper = c(100, 1, 1)))
#'
#' # Display data
#' par(mfrow = c(1, 2))
#'
#' plot(t_orig, theta, ylim = range(sapply(y, range)), type = 'l', lwd = 2, lty = 2, main = 'Original growth velocities', xlab = 'Age', ylab = 'Growth velocity')
#' for (i in 1:n) lines(t_orig, y[[i]], lwd = 0.2)
#'
#' basis <- bs(x = t[[1]], knots = kts, Boundary.knots = c(0, 1))
#' plot(t_orig, basis %*% res$c, ylim = range(sapply(y, range)), type = 'l', lwd = 2, lty = 2, main = 'Warped growth velocities', xlab = 'Biological age', ylab = 'Growth velocity')
#' for (i in 1:n) lines(t_range[2] * v(res$w[,i], t[[i]], tw), y[[i]], lwd = 0.2)
#'
#' # Predicted warping functions
#' par(mfrow = c(1, 1))
#' plot(t_orig, t_orig, type = 'l', lwd = 2, lty = 2, main = 'Warping functions', xlab = 'Age', ylab = 'Biological age')
#' for (i in 1:n) lines(t_orig, t_range[2] * v(res$w[,i], t[[i]], tw), lwd = 0.2)


estimate_generic <- function(y, t, amp_cov_par, amp_cov_fct, warp_cov_par, warp_cov_fct, kts, tw, intercept = FALSE, iter = c(5, 5), use_warp_gradient = FALSE, homeomorphisms = 'no', like_optim_control = list()) {
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
  c <- spline_weights(y, t, w, tw, bdiag(Ainv), kts, intercept = intercept)

  # Construct warp derivative
  dwarp <- list()
  for (i in 1:n) {
    dwarp[[i]] <- as(dv(t[[i]], tw), "dgCMatrix")
  }


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
            gr <- function(w, t, y, tw, c, Ainv, Cinv, kts, intercept) posterior_grad(w, dwarp[[i]], t, y, tw, c, Ainv, Cinv, kts, intercept)
            warp_optim_method <- 'BFGS'
          }
          w[, i] <- optim(par = w[,i], fn = posterior, gr = gr, method = warp_optim_method, t = t[[i]], y = y[[i]], tw = tw, c = c, Ainv = Ainv[[i]], Cinv = Cinv, kts = kts, intercept = intercept)$par
          if (homeomorphisms == 'soft') w[, i] <- make_homeo(w[, i], tw)
        }
      }

      # Update spline weights
      c <- spline_weights(y, t, w, tw, bdiag(Ainv), kts, intercept = intercept)
    }

    # Likelihood estimation of parameters (outer loop)

    # Construct residual vector for given warp prediction
    Zis <- list()
    r <- y
    for (i in 1:n) {
      twarped <- v(w[, i], t[[i]], tw)
      Zis[[i]] <- Zi(twarped, dwarp[[i]], c, kts)
      r[[i]] <- as.numeric(r[[i]] - bs(twarped, knots = kts, Boundary.knots = c(0, 1), intercept = intercept) %*% c + Zis[[i]] %*% w[,i])
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
