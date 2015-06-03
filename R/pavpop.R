# TODO: Piecewise linear basis
# TODO: AUTOMATICALLY REMOVE SIGMA SCALE PARAMETER IF INCLUDED
# TODO: TEST PARALLEL LIKELIHOOD FOR MANY SAMPLES
# TODO: MAKE SINGLE CLUSTER LIKELIHOOD (EASY!)
# TODO: MAKE DETAILS SECTION IN DOCUMENTATION
# TODO: UPDATE DOCUMENTATION
# TODO: INCORPORATE SCALING OF t VALUES NATURALLY
# TODO: PYRAMID SCHEME! (2 types) (ALSO ndeps)
# TODO: AUTOMATICALLY INITIALIZE!
# TODO: MAKE HARD HOMEOMORPHIC CONSTRAINTS POSSIBLE
# TODO LATER: MAKE LARGE WARPS POSSIBLE (AND CHEAP)
# TODO LATER: ALLOW LARGE m (USING SPARSE CHOL)
# TODO: FREE WARP COVARIANCE, ALREADY POSSIBLE?!
# TODO: REML
# TODO: ALLOW SOLUTION METHODS IN COVARIANCES AND BASIS
# TODO: METHOD BASED ON fPCA
# TODO: SPARSE OVERCOMPLETE BASIS


#' Estimate parameters and predict warps for curve data
#'
#' This function does likelihood estimation in the model \deqn{y_i(t)=\theta(v(t, w_i))+x_i(t)+\epsilon_i(t)} based on iterative local linearization of the model around predictions of the random warping parameters \eqn{w_i}.
#' @param y list of \eqn{n} functional observations. Missing values are allowed.
#' @param t list of time points corresponding to y. Should be scaled to have outer endpoints at 0 and 1.
#' @param basis_fct basis function to describe the mean function.
#' @param amp_cov_par amplitude covariance parameters.
#' @param amp_cov_fct amplitude covariance matrix function.
#' @param warp_cov_par warp covariance parameters.
#' @param warp_cov_fct warp covariance matrix function.
#' @param tw anchor points for the warping parameters.
#' @param iter two-dimensional numeric consisting of number of outer and inner iterations.
#' @param use_warp_gradient logical. Should warp prediction use the exact gradient for based optimization?
#' @param smooth_warp logical. Should warping functions be based on a monotonic cubic spline?
#' @param homeomorphisms should warps be constrained to be homeomorphisms? Options are: \code{'no'}, \code{'soft'} or \code{'hard'}. 'soft' will project the prediction onto the space of homeomorphisms after each prediction. 'hard' will do the optimiziation in the constrained space (not implemented yet!).
#' @param n_cores number of cores to use.
#' @param like_optim_control list of control options for likelihood optimization. Parameters are given as \code{c(amp_cov_par, warp_cov_par)} and options include lower, upper, method, ndev (see \code{\link[stats::optim]{optim}}).
#' @keywords likelihood estimation
#' @export
#' @examples
#' # Load female growth data from the Berkeley growth study
#' t_orig <- fda::growth$age
#' y <- fda::growth$hgtf
#' m <- nrow(y)
#' n <- ncol(y)
#'
#' # Specify age rage for normalization, endpoints are 0 and 1 in normalized data
#' t_range <- c(0, 20)
#' t <- replicate(n, t_orig / t_range[2], simplify = FALSE)
#' y <- lapply(1:n, function (x) y[, x])
#'
#' # Set up basis function
#' kts <- seq(0, 1, length = 16)[1:15]
#' basis_fct <- make_basis_fct(kts = kts, intercept = TRUE, increasing = TRUE,
#'                             order = 3, boundary = c(0, 1))
#'
#' # Set up warp function
#' tw <- seq(0, 1, length = 5)[2:4]
#' warp_fct <- make_warp_fct('smooth', tw)
#'
#' # Set up covariance functions
#' warp_cov_par <- c(tau = 1)
#' warp_cov <- make_cov_fct(Brownian, noise = FALSE, param = warp_cov_par, type = 'bridge')
#'
#' amp_cov_par <- c(scale = 200, range = 1, smoothness = 2)
#' amp_cov <- make_cov_fct(Matern, noise = TRUE, param = amp_cov_par)
#'
#'
#' # Estimate in the model
#'
#' # Bounds of parameters
#' # NOTE: Prediction of velocities is only meaningful
#' #       when the smoothness parameter is > 0.5
#' lower <- c(1e-2, 1e-2, 0.5001, 1e-2)
#' upper <- c(1000, 1, 3, 1)
#'
#' res <- pavpop(y, t, basis_fct, warp_fct, amp_cov, warp_cov, iter = c(3, 3),
#'               homeomorphism = 'soft', like_optim_control = list(lower = lower,
#'                                                                 upper = upper))
#' # Plot results
#' t_p <- seq(range(t)[1], range(t)[2], length = 100)
#' t_p_orig <- t_p * 20
#'
#' # Functional fixed effect
#' theta <- basis_fct(t_p) %*% res$c
#'
#' # Display data with predictions
#' plot(t_p_orig, theta, ylim = range(y), type = 'n', main = 'Original heights and predicted',
#'      xlab = 'Age', ylab = 'Height')
#' for (i in 1:n) {
#'   points(t_orig, y[[i]], pch = 19, cex = 0.3, col = rainbow(n)[i])
#'   lines(t_p_orig, predict_curve(t_p, t[[i]], y[[i]], res$c, basis_fct, warp_fct, amp_cov, res$amp_cov_par,
#'                                 res$w[,i], deriv = FALSE),
#'         lwd = 0.5, col = rainbow(n)[i])
#'
#' }
#' lines(t_p_orig, theta, ylim = range(y), lwd = 2, lty = 2,)
#'
#' # Compute and display growth velocities
#' plot(t_p_orig, t_p_orig, ylim = c(0, 23), type = 'n', main = 'Predicted growth velocities',
#'      xlab = 'Age', ylab = 'Growth velocity')
#' for (i in 1:n) {
#'   lines(t_p_orig, predict_curve(t_p, t[[i]], y[[i]], res$c, basis_fct, warp_fct, amp_cov, res$amp_cov_par,
#'                                 res$w[,i], deriv = TRUE) / t_range[2],
#'         lwd = 0.5, col = rainbow(n)[i])
#' }
#'
#'
#' # Display predicted warping functions
#' plot(t_orig, t_orig, type = 'l', lwd = 2, lty = 2, main = 'Warping functions', xlab = 'Age', ylab = 'Biological age')
#' for (i in 1:n) lines(t_orig, t_range[2] * warp_fct(res$w[,i], t[[i]]), lwd = 0.2)
#'

pavpop <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, iter = c(5, 5), parallel = list(n_cores = 1, parallel_likelihood = FALSE), use_warp_gradient = FALSE, homeomorphisms = 'no', like_optim_control = list()) {
  nouter <- iter[1] + 1
  if (is.null(amp_cov) & is.null(warp_cov)) nouter <- 1
  ninner <- iter[2]
  halt_iteration <- FALSE
  # Set size parameters
  n <- length(y)
  m <- sapply(y, length)

  # Warp parameters
  tw <- attr(warp_fct, 'tw')
  mw <- attr(warp_fct, 'mw')
  if (all(is.na(tw))) tw <- rep(tw, mw)
  warp_type <- attr(warp_fct, 'type')
  if (warp_type != 'piecewise linear' & warp_type != 'smooth') homeomorphisms <- 'no'

  # Unknown parameters
  amp_cov_par <- eval(attr(amp_cov, 'param'))
  warp_cov_par <- eval(attr(warp_cov, 'param'))
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

  # Initialize cluster
  registerDoParallel(cores = parallel$n_cores)

  # Initialize warp parameters
  w <- array(attr(warp_fct, 'init'), dim = c(mw, n))

  # Build amplitude covariances and inverse covariances
  inv_amp_cov <- attr(amp_cov, 'inv_cov_fct')
  inv_amp <- !is.null(attr(amp_cov, 'inv_cov_fct'))

  S <- Sinv <- list()
  for (i in 1:n) {
    # Check if an amplitude covariance is defined
    if (!is.null(amp_cov)) {
      S[[i]] <- amp_cov(t[[i]], amp_cov_par)
      if (inv_amp) {
        Sinv[[i]] <- inv_amp_cov(t[[i]], amp_cov_par)
      } else {
        Sinv[[i]] <- chol2inv(chol(S[[i]]))
      }
    } else {
      S[[i]] <- Sinv[[i]] <- Diagonal(m[i], x = 1)
    }
  }

  # Build warp covariance and inverse
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, warp_cov_par)
    Cinv <- solve(C)
  } else {
    C <- Cinv <- matrix(0, mw, mw)
  }

  # Estimate spline weights
  c <- spline_weights(y, t, warp_fct, w, Sinv, basis_fct)

  # Construct warp derivative
  dwarp <- list()
  if (warp_type != 'smooth') {
    for (i in 1:n) {
      dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
      if (warp_type == 'piecewise linear') dwarp[[i]] <- as(dwarp[[i]], "dgCMatrix")
    }
  }

  # Initialize best parameters
  like_best <- Inf
  w_best <- w
  c_best <- c
  amp_cov_par_best <- amp_cov_par
  warp_cov_par_best <- warp_cov_par

  cat('Outer\t:\tInner \t:\tEstimates\n')
  for (iouter in 1:nouter) {
    if (halt_iteration & iouter != nouter) next
    # Outer loop
    if (iouter != nouter) cat(iouter, '\t:\t')
    for (iinner in 1:ninner) {
      # Inner loop
      if (iouter != nouter | nouter == 1) cat(iinner, '\t')

      # Predict warping parameters for all functional samples
      warp_change <- c(0, 0)
      if (homeomorphisms == 'hard') {
        #TODO: constrainOptim
        stop("Hard homeomorphic constrained optimization for warps is not implemented.")
      } else {
        # Parallel prediction of warping parameters
        w_res <- foreach(i = 1:n) %dopar% {
          gr <- NULL
          warp_optim_method <- 'Nelder-Mead'
          ww <- optim(par = w[, i], fn = posterior, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = c, Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
          if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
          return(ww)
        }

        for (i in 1:n) {
          warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
          warp_change[2] <- max(warp_change[2], abs(w[, i] -  w_res[[i]]))
        }

        for (i in 1:n) w[, i] <- w_res[[i]]
      }

      # Update spline weights
      c <- spline_weights(y, t, warp_fct, w, Sinv, basis_fct)
      if (warp_change[2] < 1e-2 / sqrt(mw)) break #TODO: Consider other criteria
    }

    # Likelihood estimation of parameters (outer loop)

    # Construct residual vector for given warp prediction
    Zis <- list()
    r <- y
    for (i in 1:n) {
      twarped <- warp_fct(w[, i], t[[i]])
      if (!is.null(warp_cov)) {
        if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
        Zis[[i]] <- matrix(Zi(twarped, dwarp[[i]], basis_fct, c), m[i], mw)
      } else {
        Zis[[i]] <- Matrix(0, m[i], mw)
      }
      r[[i]] <- as.numeric(r[[i]] - basis_fct(twarped) %*% c + Zis[[i]] %*% w[, i])
    }

    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      # Likelihood function
      like_fct <- function(par) {
        like(par, n_par = c(n_par_amp, n_par_warp), r = r, Zis = Zis, amp_cov = amp_cov, warp_cov = warp_cov, t = t, tw = tw)
      }

      # Likelihood gradient
      like_gr <- NULL
      if (parallel$parallel_likelihood) {
        # Construct parallel gradient
        like_gr <- function(par) {
          epsilon <- 1e-5
          rep(1:length(par), each = 2)
          res <- foreach(ip = 1:length(par), .combine = 'c') %:%
            foreach(sign = c(1, -1), .combine= '-') %dopar% {
              h <- rep(0, length(par))
              h[ip] <- sign * epsilon
              return(like_fct(par + h) / (2 * epsilon))
            }
          return(res)
        }
      } else {
        # Construct parallel gradient
        like_gr <- function(par) {
          epsilon <- 1e-5
          rep(1:length(par), each = 2)
          res <- rep(0, length(par))
          for (ip in  1:length(par)) {
            for (sign in c(1, -1)) {
              h <- rep(0, length(par))
              h[ip] <- sign * epsilon
              res[ip] <- res[ip] + sign * like_fct(par + h) / (2 * epsilon)
            }
          }
          return(res)
        }
      }

      # Estimate parameters using locally linearized likelihood
      lower  <- if (is.null(like_optim_control$lower)) rep(1e-3, n_par_amp + n_par_warp) else like_optim_control$lower
      upper  <- if (is.null(like_optim_control$upper)) rep(Inf, n_par_amp + n_par_warp) else like_optim_control$upper
      method <- if (is.null(like_optim_control$method)) "L-BFGS-B" else like_optim_control$method
      ndeps <- if (is.null(like_optim_control$ndeps)) rep(1e-3, n_par_amp + n_par_warp) else like_optim_control$ndeps

      like_optim <- optim(c(amp_cov_par, warp_cov_par), like_fct, gr = like_gr, method = method, lower = lower, upper = upper, control = list(ndeps = ndeps, maxit = 20))
      param <- like_optim$par

      if (!is.null(amp_cov)) amp_cov_par <- param[1:n_par_amp]
      if (!is.null(warp_cov)) warp_cov_par <- param[(n_par_amp + 1):length(param)]

      if (like_optim$value <= like_best) {
        # Save parameters
        like_best <- like_optim$value
        w_best <- w
        c_best <- c
        amp_cov_par_best <- amp_cov_par
        warp_cov_par_best <- warp_cov_par

        # Update covariances
        S <- Sinv <- list()
        for (i in 1:n) {
          # Check if an amplitude covariance is defined
          if (!is.null(amp_cov)) {
            S[[i]] <- amp_cov(t[[i]], amp_cov_par)
            if (inv_amp) {
              Sinv[[i]] <- inv_amp_cov(t[[i]], amp_cov_par)
            } else {
              Sinv[[i]] <- chol2inv(chol(S[[i]]))
            }
          } else {
            S[[i]] <- Sinv[[i]] <- Diagonal(m[i], x = 1)
          }
        }

        if (!is.null(warp_cov)) {
          C <- warp_cov(tw, warp_cov_par)
          Cinv <- solve(C)
        } else {
          C <- Cinv <- matrix(0, mw, mw)
        }
        cat(':\t', param, '\n')
        cat('Linearized likelihood:\t', like_best, '\n')
      } else {
        cat(':\tLikelihood not improved, returning best likelihood estimates.\n')
        halt_iteration <- TRUE
      }
    } else {
      # TODO: Should in principle be done before warps are updated in the final iteration!
      # Estimate of sigma if final iteration is reached
      if (nouter == 1) {
        w_best <- w
        c_best <- c
      }
      sigma <- sqrt(sigmasq(c(amp_cov_par, warp_cov_par), c(n_par_amp, n_par_warp), r, Zis, amp_cov, warp_cov, t, tw))
    }
  }
  return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, sigma = sigma))
}
