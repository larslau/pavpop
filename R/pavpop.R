#' pavpop: A package for simultaneously analyzing functional data
#' containing systematic phase and amplitude variation.
#'
#' @docType package
#' @author Lars Lau Raket \email{larslau@@math.ku.dk}
#' @references L.L. Raket, S. Sommer, and B. Markussen, “A nonlinear mixed-effects model for simultaneous smoothing and registration of functional data,” Pattern Recognition Letters, vol. 38, pp. 1-7, 2014.
#'
#' @name pavpop-package
NULL
#> NULL



# TODO: MLE.Matern
# TODO: Initialize warp parameter variance in first iteration
# TODO: WARP OPTIM METHOD!
# TODO: DTW
# TODO: SITAR MODEL
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
# TODO: REML
# TODO: METHOD BASED ON fPCA
# TODO: SPARSE SOLUTION FOR OVERCOMPLETE BASIS

#' Estimate parameters and predict warps for curve data
#'
#' This function does likelihood estimation in the model \deqn{y_i(t)=\theta(v(t, w_i))+x_i(t)+\epsilon_i(t)} based on iterative local linearization of the model around predictions of the random warping parameters \eqn{w_i}.
#' @param y list of \eqn{n} functional observations. Missing values are allowed.
#' @param t list of time points corresponding to y. Should be scaled to have outer endpoints at 0 and 1.
#' @param basis_fct basis function to describe the mean function.
#' @param warp_fct warp function that models the warping procedure.
#' @param amp_cov amplitude covariance matrix function. If `NULL`, the amplitude is assumed to only contain iid Gaussian noise.
#' @param warp_cov warp covariance matrix function. If `NULL` the warps are treated as fixed parameters.
#' @param amp_fct functional basis that amplitude variation should be expressed in. If used, pavpop will automatically assume that iid. Gaussian noise is present in the data. Default is \code{NULL}.
#' @param warped_amp logical. Does the amplitude correlation follow observed or warped time?
#' @param iter two-dimensional numeric consisting of number of outer and inner iterations.
#' @param parallel list containing elements \code{n_cores}, the number of cores to use when
#' predicting warps and \code{parallel_likelihood}, a logical indicating whether likelihood
#' optimization should use a parallel evaluation of the gradient. In order for the parallel
#' likelihood to give a speed-up, large data sizes are required.
#' @param use_warp_gradient logical. Should warp prediction use the exact gradient for based optimization?
#' @param warp_optim_method optimization method for predicting warp. Defaults to 'CG' which gives robust results. 'Nelder-Mead' is often much faster, but not as reliable.
#' @param homeomorphisms should warps be constrained to be homeomorphisms? Options are: \code{'no'}, \code{'soft'} or \code{'hard'}. 'soft' will project the prediction onto the space of homeomorphisms after each prediction. 'hard' will do the optimiziation in the constrained space (not implemented yet!).
#' @param like_optim_control list of control options for likelihood optimization. Parameters are given as \code{c(amp_cov_par, warp_cov_par)} and options include lower, upper, method, ndev (see \code{\link[stats::optim]{optim}}).
#' @keywords likelihood estimation
#' @seealso See the vignettes for many more examples.
#' @export
#' @examples
#' # Load male growth data from the Berkeley growth study
#' t <- fda::growth$age
#' y <- fda::growth$hgtm
#' m <- nrow(y)
#' n <- ncol(y)
#'
#' # Specify age rage for controlling boundary points
#' t_range <- c(0, 20)
#' t <- replicate(n, t, simplify = FALSE)
#' y <- lapply(1:n, function (x) y[, x])
#'
#' # Set up basis function
#' kts <- seq(t_range[1], t_range[2], length = 15)
#' basis_fct <- make_basis_fct(kts = kts, type = 'increasing', intercept = TRUE,
#'                             control = list(boundary = t_range))
#'
#' # Set up warp function
#' tw <- seq(t_range[1], t_range[2], length = 6)
#' warp_fct <- make_warp_fct('smooth', tw, control = list(wright = 'extrapolate'))
#' mw <- attr(warp_fct, 'mw')
#'
#' # Set up covariance functions
#' warp_cov_par <- c(tau = 10)
#' warp_cov <- make_cov_fct(Brownian, noise = FALSE, param = warp_cov_par, type = 'motion',
#'                          range = t_range)
#'
#' amp_cov_par <- c(scale = 200, range = 10, smoothness = 2)
#' amp_cov <- make_cov_fct(Matern, noise = TRUE, param = amp_cov_par)
#'
#'
#' # Estimate in the model
#'
#' # Bounds of parameters
#' # NOTE: Prediction of velocities is only meaningful
#' #       when the smoothness parameter is > 0.5
#' lower <- c(1e-2, 1e-2, 0.5001, 1e-2)
#' upper <- c(1000, Inf, Inf, Inf)
#'
#' res <- pavpop(y, t, basis_fct, warp_fct, amp_cov, warp_cov, homeomorphisms = 'soft',
#'               like_optim_control = list(lower = lower, upper = upper))
#' #
#' # Plot results
#' #
#'
#' t_p <- seq(range(t)[1], range(t)[2], length = 100)
#'
#' # Functional fixed effect
#' theta <- basis_fct(t_p) %*% res$c
#'
#' # Display data with predictions
#' plot(t_p, theta, ylim = range(y), type = 'n', main = 'Original heights and predicted',
#'      xlab = 'Age (years)', ylab = 'Height (cm)')
#' for (i in 1:n) {
#'   points(t[[i]], y[[i]], pch = 19, cex = 0.3, col = rainbow(n)[i])
#'   lines(t_p, predict_curve(t_p, t[[i]], y[[i]], basis_fct, res$c, warp_fct, res$w[, i],
#'                            amp_cov, res$amp_cov_par),
#'         lwd = 0.5, col = rainbow(n)[i])
#'
#' }
#' lines(t_p, theta, ylim = range(y), lwd = 2, lty = 2)
#'
#' # Compute and display growth velocities
#' plot(t_p, t_p, ylim = c(0, 23), type = 'n', main = 'Predicted growth velocities',
#'      xlab = 'Age (years)', ylab = 'Growth velocity (cm/year)')
#' for (i in 1:n) {
#'   lines(t_p, predict_curve(t_p, t[[i]], y[[i]], basis_fct, res$c, warp_fct, res$w[, i],
#'                            amp_cov, res$amp_cov_par, deriv = TRUE),
#'         lwd = 0.5, col = rainbow(n)[i])
#' }
#'
#'
#' # Display predicted warping functions
#' plot(t_p, t_p, type = 'l', lwd = 2, lty = 2, main = 'Warping functions',
#'      xlab = 'Age (years)', ylab = 'Biological age (years)')
#' for (i in 1:n) lines(t[[i]], warp_fct(res$w[,i], t[[i]]), lwd = 0.4, col = rainbow(n)[i])

pavpop <- function(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, amp_fct = NULL, warped_amp = FALSE, iter = c(5, 5), parallel = list(n_cores = 1, parallel_likelihood = FALSE), use_warp_gradient = FALSE, warp_optim_method = 'CG', homeomorphisms = 'no', like_optim_control = list()) {

  if (!is.null(amp_fct) & is.null(amp_cov))
    warning('Amplitude variation basis ignored when no amplitude covariance is specified.')

  warp_centering <- TRUE
  # If the warps are not regularized, we should be careful when centering.
  if (is.null(warp_cov)) warp_centering <- FALSE

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
  if (all(is.na(tw))) tw <- rep(tw, mw + 2)
  tw_int <- tw[2:(mw + 1)]

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
  # Stored warped time
  t_warped <- t

  # Update m with cleaned data
  m <- sapply(y, length)

  # Initialize cluster
  registerDoParallel(cores = parallel$n_cores)

  # Initialize warp parameters
  w <- array(0, dim = c(mw, n))

  # Compute amplitude precision matrices
  Sinv <- fill_precision(t, amp_cov, amp_cov_par, amp_fct)


  # Build warp covariance and inverse
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw_int, warp_cov_par)
    Cinv <- solve(C)
  } else {
    C <- Cinv <- matrix(0, mw, mw)
  }

  # Estimate spline weights
  c <- spline_weights(y, t, Sinv, basis_fct)

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
          ww <- optim(par = w[, i], fn = posterior, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = c, Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
          if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
          return(ww)
        }

        for (i in 1:n) {
          warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
          warp_change[2] <- max(warp_change[2], abs(w[, i] -  w_res[[i]]))
        }

        # Fix warps that are considerably beyond basis-function endpoints
        # Pre-compute warped time
        t_warped <- lapply(1:n, function(i) warp_fct(w_res[[i]], t[[i]]))
        max_outside <- 0.75
        outside_frac <- sapply(t_warped, function(t) mean(attr(basis_fct, 'boundary')[1] < t & t < attr(basis_fct, 'boundary')[2]))
        for (i in which(outside_frac < max_outside)) w_res[[i]][] <- 0

        for (i in 1:n) w[, i] <- w_res[[i]]

        # Center warps
        if (warp_centering & iinner != ninner) w <- w - rowMeans(w)
      }
      # Update precisions, if
      if (warped_amp) Sinv <- fill_precision(t_warped, amp_cov, amp_cov_par, amp_fct)

      # Update spline weights
      c <- spline_weights(y, t_warped, Sinv, basis_fct)
      if (warp_change[2] < 1e-2 / sqrt(mw)) break #TODO: Consider other criteria
    }

    #
    # Likelihood estimation of parameters (outer loop)
    #

    # Construct residual vector for given warp prediction
    if (warp_type == 'smooth') {
      for (i in 1:n) dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
    }
    if (!is.null(warp_cov)) {
      Zis <- Zis(t_warped, dwarp, basis_fct, c)
    } else {
      Zis <- lapply(m, function(m) Matrix(0, m, mw))
    }
    r <- y
    for (i in 1:n) {
      r[[i]] <- as.numeric(r[[i]] - basis_fct(t_warped[[i]]) %*% c + Zis[[i]] %*% w[, i])
    }

    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      t_like <- t
      if (warped_amp) t_like <- t_warped
      # Likelihood function
      like_fct <- function(par) {
        if (!is.null(amp_fct)) {
          like_amp(par, n_par = c(n_par_amp, n_par_warp), r = r, Zis = Zis, amp_cov = amp_cov, warp_cov = warp_cov, amp_fct = amp_fct, t = t_like, tw = tw_int)
        } else {
          like(par, n_par = c(n_par_amp, n_par_warp), r = r, Zis = Zis, amp_cov = amp_cov, warp_cov = warp_cov, t = t_like, tw = tw_int)
        }
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

      if (!is.null(amp_cov)) {
        # Handle type 'unstr_cov' which can be used in combination with an
        # amplitude function. In this case, the found variance parameters m
        # may not correspond to the nearest positive definite matrix that
        # was actually used in computations
        if (!is.null(amp_fct) & attr(amp_cov, 'type') == 'unstr_cov') {
          amp_cov_mat <- warp_cov(1:n_par_amp, param[1:n_par_amp])
          amp_cov_par <- c(diag(amp_cov_mat), warp_cov_mat[upper.tri(amp_cov_mat)])
        } else {
          amp_cov_par <- param[1:n_par_amp]
        }
      }
      if (!is.null(warp_cov)) {
        # Handle type 'unstr_cov' where found parameters may not correspond to the nearest
        # positive definite matrix that was actually used in computations
        if (attr(warp_cov, 'type') == 'unstr_cov') {
          warp_cov_mat <- warp_cov(tw_int, param[(n_par_amp + 1):length(param)])
          warp_cov_par <- c(diag(warp_cov_mat), warp_cov_mat[upper.tri(warp_cov_mat)])
        } else {
          warp_cov_par <- param[(n_par_amp + 1):length(param)]
        }
      }


      if (like_optim$value <= like_best) {
        # Save parameters
        like_best <- like_optim$value
        w_best <- w
        c_best <- c
        amp_cov_par_best <- amp_cov_par
        warp_cov_par_best <- warp_cov_par

        # Update covariances
        Sinv <- fill_precision(t, amp_cov, amp_cov_par, amp_fct)

        if (!is.null(warp_cov)) {
          C <- warp_cov(tw_int, warp_cov_par)
          Cinv <- solve(C)
        }

        cat(':\t', param, '\n')
        cat('Linearized likelihood:\t', like_best, '\n')
      } else {
        cat('.\nLikelihood not improved, returning best likelihood estimates.\n')
        halt_iteration <- TRUE
      }
    } else {
      # TODO: Should in principle be done before warps are updated in the final iteration!
      # Estimate of sigma if final iteration is reached
      if (nouter == 1) {
        w_best <- w
        c_best <- c
      }

      t_like <- t
      if (warped_amp) t_like <- t_warped

      if (!is.null(amp_fct)) {
        sigma <- sqrt(sigmasq_amp(c(amp_cov_par, warp_cov_par), c(n_par_amp, n_par_warp), r, Zis, amp_cov, warp_cov, amp_fct, t_like, tw_int))
      } else {
        sigma <- sqrt(sigmasq(c(amp_cov_par, warp_cov_par), c(n_par_amp, n_par_warp), r, Zis, amp_cov, warp_cov, t_like, tw_int))
      }
    }
  }
  return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, sigma = sigma, like = like_best))
}
