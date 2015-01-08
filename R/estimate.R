#' Estimate parameters and predict warps
#'
#' This function does likelihood estimation in the model \deqn{y_i(t)=\theta(v(t, \boldsymbol{w}_i))+x_i(t)+\varepsilon_i(t)} based on iterative local linearization of the model around predictions of the random warping parameters \eqn{\boldsymbol{w}_i}.
#' @param y list of \eqn{n} functional observations. Missing values are allowed.
#' @param t list of time points corresponding to y.
#' @param tw anchor points for the warping parameters.
#' @param tau initialization of warp variance.
#' @param scale initialization of scale of serially correlated variation.
#' @param range initialization of range of serially correlated variation.
#' @keywords likelihood estimation
#' @export
#' @examples
#' #TODO

# TODO: BOUNDARY KNOTS
# TODO: WARPS WITH/WITHOUT GRADIENT
# TODO: PYRAMID SCHEME!
# TODO: AUTOMATICALLY INITIALIZE!
# TODO: MAKE ARBITRARY COVARIANCES POSSIBLE
# TODO: MAKE LARGE WARPS POSSIBLE (AND CHEAP)
# TODO: MAKE HOMEOMORPHIC CONSTRAINTS POSSIBLE
# TODO: ALLOW LARGE m (USING SPARSE CHOL)
# TODO: NOISE/NO NOISE
# TODO: CONTROL OPTIM PARAMETERS

estimate <- function(y, t, kts, tw, tau, scale, range, smoothness = 2, nouter = 5, ninner = 5) {
  # Set size parameters
  n <- length(y)
  m <- sapply(y, length)
  mw <- length(tw)

  # Check for same data structures of y and t
  if (length(t) != n) stop("y and t must have same length.")
  if (sapply(t, length)!= m) stop("Observations in y and t must have same length.")

  # Remove missing values
  for (i in 1:n) {
    missing_indices <- is.na(y[[i]])
    y[[i]] <- y[[i]][missing_indices]
    t[[i]] <- t[[i]][missing_indices]
  }
  # Update m with cleaned data
  m <- sapply(y, length)


  # Initialize warp parameters
  w <- array(0, dim = c(mw, n))

  # Build amplitude covariances and inverse covariances
  S <- Ainv <- list()
  for (i in 1:n) {
    S[[i]] <- Matern_cov(t[[i]], scale, range, smoothness)
    Ainv[[i]] <- chol2inv(chol(diag(nrow = m[i]) + S[[i]]))
  }

  # Build warp covariance and inverse
  C <- Brownian_cov(t = tw, type = 'bridge')
  Cinv <- solve(C)

  # Estimate spline weights
  c <- spline_weights(y, t, w, tw, bdiag(Ainv), kts)

  # Construct warp derivative
  dwarp <- list()
  for (i in 1:n) {
    dwarp[[i]] <- as(Matrix(dv(t[[i]], tw)), "dgCMatrix")
  }


  for (iouter in 1:nouter) {
    # Outer loop
    cat(iouter, ':\t')

    for (iinner in 1:ninner) {
      # Inner loop
      cat(iinner, '\t')

      # Predict warping parameters for all functional samples
      for (i in 1:n) {
        w[,i] <- optim(par = w[,i], fn = posterior, t = t[[i]], y = y[[i]], tw = tw, c = c, Ainv = Ainv[[i]], tau = tau, Cinv = Cinv, kts = kts)$par
      }

      # Update spline weights
      c <- spline_weights(y, t, w, tw, bdiag(Ainv), kts)
    }

    cat('\n')

    # Likelihood estimation of parameters (outer loop)

    # Construct residual vector for given warp prediction
    Zis <- list()
    r <- y
    for (i in 1:n) {
      twarped <- v(w[, i], t[[i]], tw)
      Zis[[i]] <- Zi(twarped, dwarp[[i]], c, kts)
      r[[i]] <- r[[i]] - bs(twarped, knots = kts) %*% c + Zis[[i]] %*% w[,i]
    }

    # Check wheter the final outer loop has been reached
    if (iouter != nouter) {
      # Estimate parameters using locally linearized likelihood
      param <- optim(c(tau, scale, range), like, r = r, Zis = Zis, Cinv = Cinv, smoothness = smoothness, method = "L-BFGS-B", lower = c(1e-3, 1e-3, 1e-07), upper = c(1e6, 1e6, 5), control = list(maxit = 10, trace = 0, parscale = c(tau, scale, range)))$par

      tau <- param[1]
      scale <- param[2]
      range <- param[3]

      # Update covariances
      for (i in 1:n) {
        S[[i]] <- Matern_cov(t[[i]], scale, range, smoothness)
        Ainv[[i]] <- chol2inv(chol(diag(nrow = m[i]) + S[[i]]))
      }
    } else {
      # Estimate of sigma if final iteration is reached
      sigma <- sqrt(sigmasq(c(tau, scale, range), r, Zis, Cinv, smoothness))
    }
  }

  return(list(c = c, w = w, tau = tau, scale = scale, range = range, sigma = sigma))
}
