# TODO: MAKE SINGLE CLUSTER LIKELIHOOD (EASY!)
# TODO: MAKE DETAILS SECTION IN DOCUMENTATION
# TODO: UPDATE DOCUMENTATION
# TODO: INCORPORATE SCALING OF t VALUES NATURALLY
# TODO: SUPPLY WARPING FUNCTION
# TODO: PYRAMID SCHEME! (2 types) (ALSO ndeps)
# TODO: AUTOMATICALLY INITIALIZE!
# TODO: MAKE HARD HOMEOMORPHIC CONSTRAINTS POSSIBLE
# TODO LATER: MAKE LARGE WARPS POSSIBLE (AND CHEAP)
# TODO LATER: ALLOW LARGE m (USING SPARSE CHOL)
# TODO LATER: ALLOW SHIFT WARPS
# TOTO: CHOOSE BEST LINEARIZED LIKELIHOOD

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
#' @param intercept should the basis include an intercept (default is \code{FALSE})
#' @param iter two-dimensional numeric consisting of number of outer and inner iterations.
#' @param use_warp_gradient logical. Should warp prediction use the exact gradient for based optimization?
#' @param smooth_warp logical. Should warping functions be based on a monotonic cubic spline?
#' @param homeomorphisms should warps be constrained to be homeomorphisms? Options are: \code{'no'}, \code{'soft'} or \code{'hard'}. 'soft' will project the prediction onto the space of homeomorphisms after each prediction. 'hard' will do the optimiziation in the constrained space (not implemented yet!).
#' @param cluster_options a list of control parameters if clustering should be used. See 'Details'.
#' @param like_optim_control list of control options for likelihood optimization. Parameters are given as \code{c(amp_cov_par, warp_cov_par)} and options include lower, upper, method, ndev (see \code{\link[stats::optim]{optim}}).
#' @keywords likelihood estimation
#' @export
#' @examples
#' Load female growth data from the Berkeley growth study
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
#' # Set up covariance functions
#' warp_cov_par <- c(tau = 1)
#' warp_cov_fct <- make_cov_fct(Brownian, noise = FALSE, type = 'bridge')
#'
#' amp_cov_par <- c(scale = 50, range = 1, smoothness = 1.5)
#' amp_cov_fct <- make_cov_fct(Matern, noise = TRUE)
#'
#' # Set up parametrization
#' tw <- seq(0, 1, length = 5)[2:4]
#'
#' # Estimate in the model
#'
#' # Bounds of parameters
#' # NOTE: Prediction of velocities is only meaningful
#' #       when the smoothness parameter is > 0.5
#' lower <- c(1e-2, 1e-2, 0.5001, 1e-2)
#' upper <- c(200, 1, 3, 1)
#'
#' res <- estimate_generic(y, t, basis_fct, amp_cov_par, amp_cov_fct,
#'                         warp_cov_par, warp_cov_fct, tw, iter = c(3, 3),
#'                         smooth_warp = TRUE, homeomorphism = 'soft',
#'                         like_optim_control = list(lower = lower,
#'                                                   upper = upper))
#'
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
#'   lines(t_p_orig, predict_curve(t_p, t[[i]], y[[i]], res$c, basis_fct, amp_cov_fct, res$amp_cov_par,
#'                                 res$w[,i], tw),
#'         lwd = 0.5, col = rainbow(n)[i])
#' }
#' lines(t_p_orig, theta, ylim = range(y), lwd = 2, lty = 2,)
#'
#' # Compute and display growth velocities
#' plot(t_p_orig, t_p_orig, ylim = c(0, 23), type = 'n', main = 'Predicted growth velocities',
#'      xlab = 'Age', ylab = 'Growth velocity')
#' for (i in 1:n) {
#'   lines(t_p_orig, predict_curve(t_p, t[[i]], y[[i]], res$c, basis_fct, amp_cov_fct, res$amp_cov_par,
#'                                 res$w[,i], tw, deriv = TRUE) / t_range[2],
#'         lwd = 0.5, col = rainbow(n)[i])
#' }
#'
#'
#' # Display predicted warping functions
#' plot(t_orig, t_orig, type = 'l', lwd = 2, lty = 2, main = 'Warping functions', xlab = 'Age', ylab = 'Biological age')
#' for (i in 1:n) lines(t_orig, t_range[2] * v(res$w[,i], t[[i]], tw, smooth = TRUE), lwd = 0.2)
#'



estimate_generic <- function(y, t, basis_fct, amp_cov_par, amp_cov_fct, warp_cov_par, warp_cov_fct, tw, iter = c(5, 5), use_warp_gradient = FALSE, smooth_warp = FALSE, homeomorphisms = 'no', cluster_options = list(), like_optim_control = list()) {
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

  # Cluster options

  # Number of mean functions
  if (is.null(cluster_options$clusters)) {
    k <- 1
  } else {
    k <- cluster_options$clusters
  }

  # Number of outer clustering iterations
  if (is.null(cluster_options$ncluster)) {
    ncluster <- ifelse(k == 1, 1, 5)
  } else {
    ncluster <- cluster_options$ncluster
  }
  # Weights for functional samples
  return_ind_like <- TRUE
  if (is.null(cluster_options$observation_weights))  {
    observation_weights <- array(1, c(n, k))
    return_ind_like <- FALSE
  } else {
    #TODO: Check for proper size
    observation_weights <- cluster_options$observation_weights
  }

  #TODO: option for fixed, supplied clusters

  # Initialize warp parameters
  w <- array(0, dim = c(mw, n, k))

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
  c <- list()
  for (j in 1:k) {
    c[[j]] <- spline_weights(y, t, w[,, j], tw, bdiag(Ainv), basis_fct, observation_weights[, j], smooth_warp = smooth_warp)
  }

  # Construct warp derivative
  dwarp <- list()
  for (i in 1:n) {
    dwarp[[i]] <- as(dv(t[[i]], tw), "dgCMatrix")
  }

  #TODO: DO WE NEED AN OUTER CLUSTER LOOP AND MOVE UPDATE STEP TO THAT?
  for (icluster in 1:ncluster) {
    for (iouter in 1:nouter) {
      # Outer loop
      if (iouter != nouter) cat(ifelse(ncluster != 1, paste(icluster, '\t'), ''), iouter, ':\t')

      for (iinner in 1:ninner) {
        # Inner loop
        if (iouter != nouter) cat(iinner, '\t')

        # Predict warping parameters for all functional samples
        if (homeomorphisms == 'hard') {
          #TODO: constrainOptim
          stop("Hard homeomorphic constrained optimization for warps is not implemented.")
        } else {
          for (i in 1:n) {
            for (j in 1:k) {
              gr <- NULL
              warp_optim_method <- 'Nelder-Mead'
              if (use_warp_gradient) {
                #TODO: ADD SMOOTH WARP OPTION?!
                gr <- function(w, t, y, tw, c, Ainv, Cinv, basis_fct) posterior_grad(w, dwarp[[i]], t, y, tw, c, Ainv, Cinv, basis_fct)
                warp_optim_method <- 'BFGS'
              }
              w[, i, j] <- optim(par = w[, i, j], fn = posterior, gr = gr, method = warp_optim_method, t = t[[i]], y = y[[i]], tw = tw, c = c[[j]], Ainv = Ainv[[i]], Cinv = Cinv, basis_fct = basis_fct, smooth_warp = smooth_warp)$par
              if (homeomorphisms == 'soft') w[, i, j] <- make_homeo(w[, i, j], tw)
            }
          }
        }

        # Update spline weights
        for (j in 1:k) {
          c[[j]] <- spline_weights(y, t, w[,, j], tw, bdiag(Ainv), basis_fct,observation_weights[, j], smooth_warp = smooth_warp)
        }
      }

      # Likelihood estimation of parameters (outer loop)

      # Construct residual vector for given warp prediction
      # TODO: use warp derivative directly from warping function
      Zis <- r <- list()
      for (j in 1:k) {
        Zis[[j]] <- list()
        r[[j]] <- list()
        for (i in 1:n) {
          twarped <- v(w[, i, j], t[[i]], tw, smooth = smooth_warp)
          dwarp <- matrix(0, m[i], mw)
          for (kk in 1:mw) {
            h_tmp <- rep(0, mw)
            h_tmp[kk] <- 1e-5
            dwarp[, kk] <- (v(w[, i, j] + h_tmp, t[[i]], tw, smooth = smooth_warp) - v(w[, i, j] - h_tmp, t[[i]], tw, smooth = smooth_warp)) / (2e-5)
          }
          Zis[[j]][[i]] <- Zi(twarped, dwarp, c[[j]], basis_fct)
          basis <- basis_fct(twarped)

          r[[j]][[i]] <- as.numeric(y[[i]] - basis %*% c[[j]] + Zis[[j]][[i]] %*% w[, i, j])
          i <- i + 1
        }
      }

      # Check wheter the final outer loop has been reached
      if (!(iouter == nouter & icluster == ncluster)) {
        # Estimate parameters using locally linearized likelihood
        lower  <- if (is.null(like_optim_control$lower)) rep(0, n_par_amp + n_par_warp) else like_optim_control$lower
        upper  <- if (is.null(like_optim_control$upper)) rep(Inf, n_par_amp + n_par_warp) else like_optim_control$upper
        method <- if (is.null(like_optim_control$method)) "L-BFGS-B" else like_optim_control$method
        ndeps <- if (is.null(like_optim_control$ndeps)) rep(1e-2, n_par_amp + n_par_warp) else like_optim_control$ndeps
        maxit <- if (is.null(like_optim_control$maxit)) 20 else like_optim_control$maxit

        like_fct <- function(par) {
          res <- 0
          for (j in 1:k) {
            res <- like(par, n_par = c(n_par_amp, n_par_warp), r = r[[j]], Zis = Zis[[j]], amp_cov_fct = amp_cov_fct, warp_cov_fct = warp_cov_fct, t = t, tw = tw, observation_weights = observation_weights[, j])
          }
          return(res)
        }

        param <- optim(c(amp_cov_par, warp_cov_par), like_fct, method = method, lower = lower, upper = upper, control = list(ndeps = ndeps, maxit = maxit))$par

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
        sigma_fct <- function(par) {
          res <- 0
          for (j in 1:k) {
            res <- sigmasq(par, n_par = c(n_par_amp, n_par_warp), r = r[[j]], Zis = Zis[[j]], amp_cov_fct = amp_cov_fct, warp_cov_fct = warp_cov_fct, t = t, tw = tw, observation_weights = observation_weights[, j])
          }
          return(res)
        }
        sigma <- sqrt(sigma_fct(c(amp_cov_par, warp_cov_par)))
      }
      if (iouter != nouter) cat(':\t', param, '\n')
    }
    ind_like <- array(NA, dim = c(n, k))
    for (j in 1:k) {
      for (i in 1:n) ind_like[i, j] <- like_ind(c(amp_cov_par, warp_cov_par), n_par = c(n_par_amp, n_par_warp), r = r[[j]][[i]], Zi = Zis[[j]][[i]], amp_cov_fct = amp_cov_fct, warp_cov_fct = warp_cov_fct, t = t[[i]], tw = tw)
    }
    # If multiple clusters are supplied, update weights, M step in EM algorithm
    if (k > 1) {
      pi <- colMeans(weights)
      for (j in 1:k) {
        weights[, j] <- pi[j] * ind_like[, j] / (rowSums(t(pi * t(ind_like))))
      }
    }

  }
  if (k == 1) c <- c[[1]]
  result <- list(c = c, w = w[,,1:k], amp_cov_par = amp_cov_par, warp_cov_par = warp_cov_par, sigma = sigma)
  if (return_ind_like) { #TODO: rename?
    result$ind_like <- ind_like
    result$pi <- pi
  }
  return(result)
}
