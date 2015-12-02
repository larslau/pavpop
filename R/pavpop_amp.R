# TODO: Check, bruger jeg SS og S?
# TODO: LINK TO pavpop HELP PAGE

#' Estimate parameters and predict warps for curve data
#'
#' This function does likelihood estimation in the model \deqn{y_i(t)=\theta(v(t, w_i))+x_i(t)+\epsilon_i(t)} based on iterative local linearization of the model around predictions of the random warping parameters \eqn{w_i}.
#' @param y list of \eqn{n} functional observations. Missing values are allowed.
#' @param t list of time points corresponding to y. Should be scaled to have outer endpoints at 0 and 1.
#' @param basis_fct basis function to describe the mean function.
#' @param amp_cov amplitude covariance matrix function. If `NULL`, the amplitude is assumed to only contain iid Gaussian noise.
#' @param warp_cov warp covariance matrix function. If `NULL` the warps are treated as fixed parameters.
#' @param warped_amp logical. Does the amplitude correlation follow observed or warped time?
#' @param iter two-dimensional numeric consisting of number of outer and inner iterations.
#' @param use_warp_gradient logical. Should warp prediction use the exact gradient for based optimization?
#' @param smooth_warp logical. Should warping functions be based on a monotonic cubic spline?
#' @param homeomorphisms should warps be constrained to be homeomorphisms? Options are: \code{'no'}, \code{'soft'} or \code{'hard'}. 'soft' will project the prediction onto the space of homeomorphisms after each prediction. 'hard' will do the optimiziation in the constrained space (not implemented yet!).
#' @param n_cores number of cores to use.
#' @param basis_amp logical. Should the amplitude variation be expressed in terms of the basis functions?
#' @param like_optim_control list of control options for likelihood optimization. Parameters are given as \code{c(amp_cov_par, warp_cov_par)} and options include lower, upper, method, ndev (see \code{\link[stats::optim]{optim}}).
#' @keywords likelihood estimation
#' @export

pavpop_amp <- function(y, t, basis_fct, amp_fct, warp_fct, amp_cov = NULL, warp_cov = NULL, warped_amp = FALSE, iter = c(5, 5), parallel = list(n_cores = 1, parallel_likelihood = FALSE), use_warp_gradient = FALSE, warp_optim_method = 'CG', homeomorphisms = 'no', like_optim_control = list()) {
  # Check compatability between warp_cov and basis_amp
  if (!(is.null(warp_cov) || !is.null(attr(attr(warp_cov, 'cov_fct'), 'discrete'))))
    stop('A discrete covariance structure must be used when amplitude variation
         is modeled using a functional basis.')

  # If amp_cov contains noise term, make new covariance without noise
  if (attr(amp_cov, 'noise')) {
    warning('Noise term in amplitude covariance removed.
If required, manually construct covariance function with noise term.')
    amp_cov <- make_cov_fct(attr(amp_cov, 'cov_fct'), noise = FALSE, param = attr(amp_cov, 'param'))
  }


  # SET AS PARAMETER, MAY NOT ALWAYS BE IDEAL
  # OR BETTER, MAKE CHECK TO SEE IF IT IMPROVES
  warp_centering <- TRUE

  nouter <- iter[1] + 1
  if (is.null(amp_cov) & is.null(warp_cov)) nouter <- 1
  ninner <- iter[2]
  halt_iteration <- FALSE
  # Set size parameters
  n <- length(y)
  m <- sapply(y, length)
  df <- attr(amp_fct, 'df')

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
  # Stored warped time
  t_warped <- t

  # Update m with cleaned data
  m <- sapply(y, length)

  # Initialize cluster
  registerDoParallel(cores = parallel$n_cores)

  # Initialize warp parameters
  w <- array(attr(warp_fct, 'init'), dim = c(mw, n))

  # Build amplitude covariances and inverse covariances
  inv_amp_cov <- attr(amp_cov, 'inv_cov_fct')
  inv_amp <- !is.null(attr(amp_cov, 'inv_cov_fct'))

  if (!is.null(amp_cov)) {
    SS <- amp_cov(1:df, amp_cov_par)
    SSinv <- chol2inv(chol(SS))
  }

  S <- Sinv <- list()
  for (i in 1:n) {
    # Check if an amplitude covariance is defined
    if (!is.null(amp_cov)) {
      A <- amp_fct(t[[i]])
      S[[i]] <- A %*% SS %*% t(A) + diag(1, m[i])
      Sinv[[i]] <- diag(1, m[i]) - A %*% solve(SSinv + t(A) %*% A, t(A))
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
          ww <- optim(par = w[, i], fn = posterior, gr = gr, method = warp_optim_method, warp_fct = warp_fct, t = t[[i]], y = y[[i]], c = c, Sinv = Sinv[[i]], Cinv = Cinv, basis_fct = basis_fct)$par
          if (homeomorphisms == 'soft') ww <- make_homeo(ww, tw)
          return(ww)
        }

        for (i in 1:n) {
          warp_change[1] <- warp_change[1] + sum((w[, i] - w_res[[i]])^2)
          warp_change[2] <- max(warp_change[2], abs(w[, i] -  w_res[[i]]))
        }

        for (i in 1:n) w[, i] <- w_res[[i]]
        if (warp_centering & iinner != ninner) w <- w - rowMeans(w)
      }

      # Update spline weights
      c <- spline_weights(y, t, warp_fct, w, Sinv, basis_fct)
      if (warp_change[2] < 1e-2 / sqrt(mw)) break #TODO: Consider other criteria
    }

    # Likelihood estimation of parameters (outer loop)

    # Pre-compute warped time
    for (i in 1:n) {
      t_warped[[i]] <- warp_fct(w[, i], t[[i]])
    }

    # If the amplitude variation is assumed to be varying in warped time, the amplitude covariances are updated
    if (warped_amp) {
      warning('Not currently implemented!')
      if (!is.null(amp_cov)) {
        A <- amp_fct(t[[i]])
        S[[i]] <- A %*% SS %*% t(A) + diag(1, m[i])
        Sinv[[i]] <- diag(1, m[i]) - A %*% solve(SSinv + t(A) %*% A, t(A))
      } else {
        S[[i]] <- Sinv[[i]] <- Diagonal(m[i], x = 1)
      }
    }

    # Construct residual vector for given warp prediction
    Zis <- list()
    r <- y
    for (i in 1:n) {
      if (!is.null(warp_cov)) {
        if (warp_type == 'smooth') dwarp[[i]] <- warp_fct(w[, i], t[[i]], w_grad = TRUE)
        Zis[[i]] <- matrix(Zi(t_warped[[i]], dwarp[[i]], basis_fct, c), m[i], mw)
      } else {
        Zis[[i]] <- Matrix(0, m[i], mw)
      }
      r[[i]] <- as.numeric(r[[i]] - basis_fct(t_warped[[i]]) %*% c + Zis[[i]] %*% w[, i])
    }

    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      t_like <- t
      if (warped_amp) t_like <- t_warped
      # Likelihood function
      like_fct <- function(par) {
        like_amp(par, n_par = c(n_par_amp, n_par_warp), r = r, Zis = Zis, amp_cov = amp_cov, warp_cov = warp_cov, amp_fct = amp_fct, warped_amp = warped_amp, t = t_like, tw = tw)
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
        if (!is.null(amp_cov)) {
          SS <- amp_cov(1:df, amp_cov_par)
          SSinv <- chol2inv(chol(SS))
        }

        for (i in 1:n) {
          twarped <- t[[i]]
          if (warped_amp) twarped <- t_warped[[i]]
          # Check if an amplitude covariance is defined
          if (!is.null(amp_cov)) {
            A <- amp_fct(t[[i]])
            S[[i]] <- A %*% SS %*% t(A) + diag(1, m[i])
            Sinv[[i]] <- diag(1, m[i]) - A %*% solve(SSinv + t(A) %*% A, t(A))
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
      sigma <- sqrt(sigmasq_amp(c(amp_cov_par, warp_cov_par), c(n_par_amp, n_par_warp), r, Zis, amp_cov, warp_cov, amp_fct, warped_amp, t, tw))
    }
  }
  return(list(c = c_best, w = w_best, amp_cov_par = amp_cov_par_best, warp_cov_par = warp_cov_par_best, sigma = sigma, like = like_best))
}
