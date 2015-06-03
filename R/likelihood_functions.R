#TODO: CHECK OBSERVATION WEIGHTS = 1 CALCULATIONS

#' Locally linearized likelihood function
#'
#' Computes the linearized likelihood given the residual around a given warp and the corresponding Jacobians.
#' @param param variance parameters.
#' @param n_par vector consisting of number of variance parameters for each covariance function.
#' @param r residual.
#' @param Zis list of Jacobians in the warps of the mean function around the given warp.
#' @param amp_cov function for generating amplitude covariance matrix.
#' @param warp_cov function for generating warp covariance function
#' @param t array of time variables corresponding to r.
#' @param tw anchor points for warp variables.
#' @param observation_weights vector of weights for the individual functional samples to be applied to the likelihood. This is useful for clustering analysis.
#' @keywords likelihood
#' @keywords linearization
#' @export
#' @importFrom Matrix t

like <- function(param, n_par, r, Zis, amp_cov, warp_cov, t, tw) {
  amp_cov_par <- param[1:n_par[1]]
  warp_cov_par <- param[(n_par[1] + 1):length(param)]
  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, warp_cov_par)
    Cinv <- chol2inv(chol(C))
  } else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }

  n <- length(r)
  m <- sapply(r, length)

  sq <- logdet <- 0
  for (i in 1:n) {
    if (!is.null(amp_cov)) {
      S <- amp_cov(t[[i]], amp_cov_par)
      U <- chol(S)
    } else {
      # TODO: SPARSE MATRIX COULD MAKE IT FASTER, THEN 'diag' cannot be used for logdet
      S <- U <- diag(1, m[i])
    }
    rr <- r[[i]]
    ZZ <- Zis[[i]]

    if (!is.null(warp_cov)) {
      A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
      LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
      x <- t(A) %*% rr
    } else {
      LR <- x <- 0
    }
    sq <- sq + (sum(backsolve(U, rr, transpose = TRUE)^2)
                - t(x) %*% LR %*% x)
    logdet_tmp <- 0
    if (!is.null(warp_cov)) logdet_tmp <- determinant(LR)$modulus[1]
    logdet <- logdet - (logdet_tmp - 2 * sum(log(diag(U))))
  }
  logdet <- logdet - n * determinant(Cinv)$modulus[1]

  sigmahat <- as.numeric(sq / sum(m))
  res <- sum(m) * log(sigmahat) + logdet
  return(res)
}


#' Locally linearized likelihood function cluster
#'
#' Computes the linearized likelihood given the residual around a given warp and the corresponding Jacobians.
#' @param param variance parameters.
#' @param n_par vector consisting of number of variance parameters for each covariance function.
#' @param r residual.
#' @param Zis list of Jacobians in the warps of the mean function around the given warp.
#' @param amp_cov function for generating amplitude covariance matrix.
#' @param warp_cov function for generating warp covariance function
#' @param t array of time variables corresponding to r.
#' @param tw anchor points for warp variables.
#' @param observation_weights vector of weights for the individual functional samples to be applied to the likelihood. This is useful for clustering analysis.
#' @keywords likelihood
#' @keywords linearization
#' @export
#' @importFrom Matrix t

like_clust <- function(param, n_par, r, Zis, amp_cov, warp_cov, t, tw, observation_weights) {
  if (is.null(observation_weights)) stop('Observation weights must be supplied.')
  amp_cov_par <- param[1:n_par[1]]
  warp_cov_par <- param[(n_par[1] + 1):length(param)]

  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, warp_cov_par)
    Cinv <- chol2inv(chol(C))
  } else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }

  n <- length(t)
  m <- sapply(t, length)
  k <- ncol(observation_weights)

  sq <- logdet <- 0

  for (i in 1:n) {
    # Construct covariance matrix
    if (!is.null(amp_cov)) {
      S <- amp_cov(t[[i]], amp_cov_par)
      U <- chol(S)
    } else {
      # TODO: SPARSE MATRIX COULD MAKE IT FASTER, THEN 'diag' cannot be used for logdet
      S <- U <- diag(1, m[i])
    }

    logdet_tmp <- 0

    for (j in 1:k) {
      rr <- r[[i]][[j]]
      ZZ <- Zis[[i]][[j]]

      if (!is.null(warp_cov)) {
        A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
        LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
        x <- t(A) %*% rr
        logdet_tmp <- logdet_tmp + determinant(LR)$modulus[1] * observation_weights[i, j]
      } else {
        LR <- x <- 0
      }
      sq <- sq + (sum(backsolve(U, rr, transpose = TRUE)^2)
                  - t(x) %*% LR %*% x) * observation_weights[i, j]
    }
    logdet <- logdet - (logdet_tmp - 2 * sum(log(diag(U))))
  }

  logdet <- logdet - n * determinant(Cinv)$modulus[1]

  sigmahat <- as.numeric(sq / sum(m))
  res <- sum(m) * log(sigmahat) + logdet
  return(res)
}

#' Maximum likelihood estimate of residual variance
#'
#' Computes the linearized likelihood given the residual around a given warp and the corresponding Jacobians.
#' @param param variance parameters.
#' @param n_par vector consisting of number of variance parameters for each covariance function.
#' @param r residual.
#' @param Zis list of Jacobians in the warps of the mean function around the given warp.
#' @param amp_cov function for generating amplitude covariance matrix.
#' @param warp_cov function for generating warp covariance function
#' @param t array of time variables corresponding to r.
#' @param tw anchor points for warp variables.
#' @param observation_weights vector of weights for the individual functional samples to be applied to the likelihood. This is useful for clustering analysis.
#' @keywords likelihood
#' @keywords linearization
#' @export
#' @importFrom Matrix t


sigmasq_clust <- function(param, n_par, r, Zis, amp_cov, warp_cov, t, tw, observation_weights) {
  if (is.null(observation_weights)) stop('Observation weights must be supplied.')
  amp_cov_par <- param[1:n_par[1]]
  warp_cov_par <- param[(n_par[1] + 1):length(param)]

  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, warp_cov_par)
    Cinv <- chol2inv(chol(C))
  } else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }

  n <- length(t)
  m <- sapply(t, length)
  k <- ncol(observation_weights)

  sq <- logdet <- 0

  for (i in 1:n) {
    # Construct covariance matrix
    if (!is.null(amp_cov)) {
      S <- amp_cov(t[[i]], amp_cov_par)
      U <- chol(S)
    } else {
      # TODO: SPARSE MATRIX COULD MAKE IT FASTER, THEN 'diag' cannot be used for logdet
      S <- U <- diag(1, m[i])
    }
    for (j in 1:k) {
      rr <- r[[i]][[j]]
      ZZ <- Zis[[i]][[j]]

      if (!is.null(warp_cov)) {
        A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
        LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
        x <- t(A) %*% rr
      } else {
        LR <- x <- 0
      }
      sq <- sq + (sum(backsolve(U, rr, transpose = TRUE)^2)
                  - t(x) %*% LR %*% x) * observation_weights[i, j]
    }
  }
  sigmahat <- as.numeric(sq / sum(m))

  return(sigmahat)
}

#' Individual locally linearized likelihood function
#'
#' Computes the linearized likelihood for a single functional sample given the residual around a given warp and the corresponding Jacobians.
#' @param param variance parameters.
#' @param sigma_sq maximum likelihood estiamte for sigma_sq.
#' @param n_par vector consisting of number of variance parameters for each covariance function.
#' @param r residual.
#' @param Zi Jacobian in the warps of the mean function around the given warp.
#' @param amp_cov function for generating amplitude covariance matrix.
#' @param warp_cov function for generating warp covariance function
#' @param t time variable corresponding to r.
#' @param tw anchor points for warp variables.
#' @keywords likelihood
#' @keywords linearization

ind_like <- function(param, sigma_sq, n_par, r, Zi, amp_cov, warp_cov, t, tw) {
  amp_cov_par <- param[1:n_par[1]]
  warp_cov_par <- param[(n_par[1] + 1):length(param)]

  C <- warp_cov(tw, warp_cov_par)
  Cinv <- chol2inv(chol(C))
  m <- length(r)

  sq <- logdet <- 0

  S <- amp_cov(t, amp_cov_par)
  U <- chol(S)
  A <- backsolve(U, backsolve(U, Zi, transpose = TRUE))
  LR <- chol2inv(chol(Cinv + Matrix::t(Zi) %*% A))
  x <- t(A) %*% r
  sq <- sq + (sum(backsolve(U, r, transpose = TRUE)^2)
              - t(x) %*% LR %*% x)

  logdet <- logdet - (determinant(LR)$modulus[1]
                      - 2 * sum(log(diag(U)))) - determinant(Cinv)$modulus[1]

  res <- m / 2 * log(sigma_sq) + 1 / 2 * logdet + 1 / (2 * sigma_sq) * sq
  res <- exp(-res)
  return(res)
}


#' Noise scale estimate from locally linearized likelihood function
#'
#' Computes the noise scale estimate from the linearized likelihood given the residual around a predicted warp and the corresponding Jacobians.
#' @inheritParams like
#' @keywords likelihood
#' @keywords linearization
#' @export


sigmasq <- function(param, n_par, r, Zis, amp_cov, warp_cov, t, tw) {
  amp_cov_par <- param[1:n_par[1]]
  warp_cov_par <- param[(n_par[1] + 1):length(param)]

  if (!is.null(warp_cov)) {
    C <- warp_cov(tw, warp_cov_par)
    Cinv <- chol2inv(chol(C))
  } else {
    C <- Cinv <- matrix(0, length(tw), length(tw))
  }

  n <- length(r)
  m <- sapply(r, length)

  sq <- 0
  for (i in 1:n) {
    if (!is.null(amp_cov)) {
      S <- amp_cov(t[[i]], amp_cov_par)
      U <- chol(S)
    } else {
      S <- U <- Diagonal(m[i], 1)
    }
    rr <- r[[i]]
    ZZ <- Zis[[i]]
    if (!is.null(warp_cov)) {
      A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
      LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
      x <- t(A) %*% rr
    } else {
      LR <- x <- 0
    }
    sq <- sq + (sum(backsolve(U, rr, transpose = TRUE)^2)
                - t(x) %*% LR %*% x)
  }
  sigmahat <- as.numeric(sq / sum(m))
  return(sigmahat)
}

#' Posterior of the data given the random warping parameters
#'
#' This function calculates the posterior of the data given the random warping parameters
#' @param w warping parameters.
#' @param t evaluation points.
#' @param y values at evaluation points.
#' @param basis_fct basis function to describe the mean function.
#' @param c spline coefficients.
#' @param Ainv precision matrix for amplitude variation.
#' @param Cinv precision matrix for the warping parameters.
#' @keywords warping
#' @keywords posterior
#' @export

posterior <- function(w, warp_fct, t, y, basis_fct, c, Sinv, Cinv) {
  vt <- warp_fct(w, t)
  basis <- basis_fct(vt)
  r <- y - basis %*% c
  return((t(r) %*% Sinv %*% r + t(w) %*% Cinv %*% w)[1])
}

#' Posterior of the data given the random warping parameters
#'
#' This function calculates the posterior of the data given the random warping parameters
#' @param w warping parameters.
#' @param dwarp Jaobian of the vector of observed warped points in the warping parameters.
#' @param t evaluation points.
#' @param y values at evaluation points.
#' @param tw anchor points for the warping parameters.
#' @param c B-spline coefficients.
#' @param Ainv precision matrix for amplitude variation.
#' @param Cinv precision matrix for the warping parameters.
#' @param basis_fct basis function to describe the mean function.
#' @param smooth_warp logical. Should warping functions be based on a monotonic cubic spline?
#' @keywords warping
#' @keywords posterior
#' @export

posterior_grad <- function(w, dwarp, t, y, tw, c, Ainv, Cinv, basis_fct) {
  vt <- v(w, t, tw, smooth = smooth_warp, smooth_warp = FALSE)
  vt[vt < 0] <- 0
  vt[vt > 1] <- 1
  r <- basis_fct(vt) %*% c - y
  theta_d <- basis_fct(vt, deriv = TRUE) %*% c
  grad <- 2 * t(r) %*% Ainv %*% (dwarp * theta_d[, 1])
  return(as.numeric(grad + 2 * w %*% Cinv))
}

# Minimize posterior of the data given the random warping parameters using a pyramidal primal-dual approach
# TODO: UPDATE, DON'T USE IN CURRENT FORM!!
# TODO: UPDATE DOCUMENTATION!!
# TODO: DOES NOT WORK WITH SPLINE DERIVATIVE AND INTERCEPT
#
# predict_warp_pyramid <- function(w, y, Ainv, t, tw, kts, warp_cov_par, warp_cov_fct, plevels = 5, beta0 = 1, start = 1, homeomorphic = TRUE, iter = c(10, 5), intercept = FALSE) {
#   n_outer <- iter[1]
#   n_inner <- iter[2]
#   n <- length(y)
#   m <- sapply(y, length)
#   nw <- nrow(w)
#   w_pred <- w
#
#   # Normalize y
# #   y_range <- range(sapply(y, range))
# #   y_mean <- mean(unlist(y))
# #   y <- lapply(y, function (x) (x - y_mean) / (diff(y_range)))
#
#   # Construct warp levels
#   nw_levels <- round(seq(start, nw, length = min(nw + 1 - start, plevels)))
#   twp <- tw
#
#   for (plvl in 1:plevels) {
#     w_pred_prev <- w_pred
#
#     twp_prev <- twp
#     nwp <- nw_levels[plvl]
#
#     twp <- seq(0, 1, length = nwp + 2)[2:(nwp + 1)]
#     w_pred <- array(0, dim = c(nwp, n))
#     for (i in 1:n) w_pred[, i] <- approx(c(0, twp_prev, 1), c(0, w_pred_prev[, i], 0), xout = twp)$y
#
#     # Build warp covariance
#     Cp <- warp_cov_fct(twp, warp_cov_par)
#     Cp_inv <- solve(Cp) #TODO: Optimize
#
#     # Construct warp derivative
#     dwarp <- list()
#     for (i in 1:n) {
#       dwarp[[i]] <- dv(t[[i]], twp) #as(Matrix(dv(t[[i]], twp)), "dgCMatrix")
#     }
#
#     beta <- beta0
#     kappa <- 1.5
#     for (outer in 1:n_outer) {
#       reg_Matrix <- solve(diag(1, nwp, nwp) + 1 / beta * Cp_inv)
#       c <- spline_weights(y, t, w_pred, twp, bdiag(Ainv), kts, intercept = intercept)
#       for (inner in 1:n_inner) {
#         for (i in 1:n) { #TODO: CHECK THIS CODE!!!
#           t_warped <- v(w_pred[, i], t[[i]], twp)
#           # The basis construction is also done in spline_weights
#           wbasis <- bs(t_warped, knots = kts, Boundary.knots = c(0, 1), intercept = intercept)
#           wdbasis <- bsd(t_warped, knots = kts, Boundary.knots = c(0, 1))
#           thetahat_warped <- wbasis %*% c
#           thetahat_warped_d <- wdbasis %*% c
#           b <- array(NA, dim = c(m[i], nwp))
#           for (j in 1:nwp) b[, j] <- as.numeric(-dwarp[[i]][, j] * thetahat_warped_d)
#           a <- y[[i]] - thetahat_warped - b %*% w_pred[, i]
#           w_pred[, i] <- solve(beta * diag(1, nwp) + t(b) %*% Ainv[[i]] %*% b,
#                                beta * w_pred[, i] - t(b) %*% Ainv[[i]] %*% a) #TODO: OPTIMIZE
#           w_pred[, i] <- reg_Matrix %*% w_pred[, i]
#           if (inner == n_inner & homeomorphic) w_pred[, i] <- make_homeo(w_pred[, i], twp)
#         }
#       }
#       beta <- beta * kappa
#     }
#   }
#   return(w_pred)
# }
