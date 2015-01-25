#' Locally linearized likelihood function
#'
#' Computes the linearized likelihood given the residual around a given warp and the corresponding Jacobians.
#' @param param variance parameters.
#' @param n_par vector consisting of number of variance parameters for each covariance function.
#' @param r residual.
#' @param Zis list of Jacobians in the warps of the mean function around the given warp.
#' @param amp_cov_fct function for generating amplitude covariance matrix.
#' @param warp_cov_fct function for generating warp covariance function
#' @param t array of time variables corresponding to r.
#' @param tw anchor points for warp variables.
#' @keywords likelihood
#' @keywords linearization
#' @export
#' @importFrom Matrix t

#TODO: Allow for no warping/amplitude covariance

like <- function(param, n_par, r, Zis, amp_cov_fct, warp_cov_fct, t, tw) {
  amp_cov_par <- param[1:n_par[1]]
  warp_cov_par <- param[(n_par[1] + 1):length(param)]

  C <- warp_cov_fct(tw, warp_cov_par)
  Cinv <- chol2inv(chol(C))

  n <- length(r)
  m <- sapply(r, length)

  sq <- logdet <- 0
  for (i in 1:n) {
    S <- amp_cov_fct(t[[i]], amp_cov_par)
    rr <- r[[i]]
    U <- chol(S)
    ZZ <- Zis[[i]]
    A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
    LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
    x <- t(A) %*% rr
    sq <- sq + sum(backsolve(U, rr, transpose = TRUE)^2) - t(x) %*% LR %*% x
    logdet <- logdet - determinant(LR)$modulus[1] + 2 * sum(log(diag(U)))
  }
  logdet <- logdet - n * determinant(Cinv)$modulus[1]

  sigmahat <- 1/sum(m) * as.numeric(sq)
  res <- sum(m) * log(sigmahat) + logdet

  return(res)
}

#' Noise scale estimate from locally linearized likelihood function
#'
#' Computes the noise scale estimate from the linearized likelihood given the residual around a predicted warp and the corresponding Jacobians.
#' @inheritParams like
#' @keywords likelihood
#' @keywords linearization
#' @export


sigmasq <- function(param, n_par, r, Zis, amp_cov_fct, warp_cov_fct, t, tw) {
  amp_cov_par <- param[1:n_par[1]]
  warp_cov_par <- param[(n_par[1] + 1):length(param)]

  C <- warp_cov_fct(tw, warp_cov_par)
  Cinv <- chol2inv(chol(C))

  n <- length(r)
  m <- sapply(r, length)

  sq <- 0
  for (i in 1:n) {
    S <- amp_cov_fct(t[[i]], amp_cov_par)
    rr <- r[[i]]
    U <- chol(S)
    ZZ <- Zis[[i]]
    A <- backsolve(U, backsolve(U, ZZ, transpose = TRUE))
    LR <- chol2inv(chol(Cinv + Matrix::t(ZZ) %*% A))
    x <- t(A) %*% rr
    sq <- sq + sum(backsolve(U, rr, transpose = TRUE)^2) - t(x) %*% LR %*% x
  }

  sigmahat <- 1/sum(m) * as.numeric(sq)

  return(sigmahat)
}

#' Posterior of the data given the random warping parameters
#'
#' This function calculates the posterior of the data given the random warping parameters
#' @param w warping parameters.
#' @param t evaluation points.
#' @param tw anchor points for the warping parameters.
#' @param c B-spline coefficients.
#' @param Ainv precision matrix for amplitude variation.
#' @param Cinv precision matrix for the warping parameters.
#' @param kts anchor points for the B-spline basis used to model the functional parameter \eqn{\theta}.
#' @keywords warping
#' @keywords posterior
#' @export

posterior <- function(w, t, y, tw, c, Ainv, Cinv, kts, intercept = FALSE, smooth_warp = FALSE, increasing = FALSE) {
  vt <- v(w, t, tw, smooth = smooth_warp)
  vt[vt < 0] <- 0
  vt[vt > 1] <- 1
  if (!increasing) {
    basis <- bs(vt, knots = kts, Boundary.knots = c(0, 1), intercept = intercept)
  } else {
    basis <- t(Ispline(vt, 3, kts))
    if (intercept) basis <- cbind(1, basis)
  }
  r <- y - basis %*% c
  return((t(r) %*% Ainv %*% r + t(w) %*% Cinv %*% w)[1])
}

#' Posterior of the data given the random warping parameters
#'
#' This function calculates the posterior of the data given the random warping parameters
#' @param dwarp Jaobian of the vector of observed warped points in the warping parameters.
#' @inheritParams like
#' @keywords warping
#' @keywords posterior
#' @export

posterior_grad <- function(w, dwarp, t, y, tw, c, Ainv, Cinv, kts, intercept = FALSE) {
  vt <- v(w, t, tw)
  vt[vt < 0] <- 0
  vt[vt > 1] <- 1
  r <- bs(vt, knots = kts, Boundary.knots = c(0, 1), intercept = intercept) %*% c - y
  theta_d <- bsd(vt, knots = kts, Boundary.knots = c(0, 1)) %*% c
  grad <- 2 * t(r) %*% Ainv %*% (dwarp * theta_d[, 1])
  return(as.numeric(grad + 2 * w %*% Cinv))
}

#' Minimize posterior of the data given the random warping parameters using a pyramidal primal-dual approach
#'
#' This function minimizes the posterior of the data given the random warping parameters using a primal-dual approach in a pyramidal setup
#' @param w warping parameters.
#' @param t evaluation points.
#' @param tw anchor points for the warping parameters.
#' @param c B-spline coefficients.
#' @param Ainv precision matrix for amplitude variation.
#' @param Cinv precision matrix for the warping parameters.
#' @param Boundary.knots boundary knots for the B-spline basis.
#' @keywords warping
#' @keywords posterior
#' @export

# TODO: UPDATE, DON'T USE IN CURRENT FORM!!
# TODO: UPDATE DOCUMENTATION!!
# TODO: DOES NOT WORK WITH SPLINE DERIVATIVE AND INTERCEPT

predict_warp_pyramid <- function(w, y, Ainv, t, tw, kts, warp_cov_par, warp_cov_fct, plevels = 5, beta0 = 1, start = 1, homeomorphic = TRUE, iter = c(10, 5), intercept = FALSE) {
  n_outer <- iter[1]
  n_inner <- iter[2]
  n <- length(y)
  m <- sapply(y, length)
  nw <- nrow(w)
  w_pred <- w

  # Normalize y
#   y_range <- range(sapply(y, range))
#   y_mean <- mean(unlist(y))
#   y <- lapply(y, function (x) (x - y_mean) / (diff(y_range)))

  # Construct warp levels
  nw_levels <- round(seq(start, nw, length = min(nw + 1 - start, plevels)))
  twp <- tw

  for (plvl in 1:plevels) {
    w_pred_prev <- w_pred

    twp_prev <- twp
    nwp <- nw_levels[plvl]

    twp <- seq(0, 1, length = nwp + 2)[2:(nwp + 1)]
    w_pred <- array(0, dim = c(nwp, n))
    for (i in 1:n) w_pred[, i] <- approx(c(0, twp_prev, 1), c(0, w_pred_prev[, i], 0), xout = twp)$y

    # Build warp covariance
    Cp <- warp_cov_fct(twp, warp_cov_par)
    Cp_inv <- solve(Cp) #TODO: Optimize

    # Construct warp derivative
    dwarp <- list()
    for (i in 1:n) {
      dwarp[[i]] <- dv(t[[i]], twp) #as(Matrix(dv(t[[i]], twp)), "dgCMatrix")
    }

    beta <- beta0
    kappa <- 1.5
    for (outer in 1:n_outer) {
      reg_Matrix <- solve(diag(1, nwp, nwp) + 1 / beta * Cp_inv)
      c <- spline_weights(y, t, w_pred, twp, bdiag(Ainv), kts, intercept = intercept)
      for (inner in 1:n_inner) {
        for (i in 1:n) { #TODO: CHECK THIS CODE!!!
          t_warped <- v(w_pred[, i], t[[i]], twp)
          # The basis construction is also done in spline_weights
          wbasis <- bs(t_warped, knots = kts, Boundary.knots = c(0, 1), intercept = intercept)
          wdbasis <- bsd(t_warped, knots = kts, Boundary.knots = c(0, 1))
          thetahat_warped <- wbasis %*% c
          thetahat_warped_d <- wdbasis %*% c
          b <- array(NA, dim = c(m[i], nwp))
          for (j in 1:nwp) b[, j] <- as.numeric(-dwarp[[i]][, j] * thetahat_warped_d)
          a <- y[[i]] - thetahat_warped - b %*% w_pred[, i]
          w_pred[, i] <- solve(beta * diag(1, nwp) + t(b) %*% Ainv[[i]] %*% b,
                               beta * w_pred[, i] - t(b) %*% Ainv[[i]] %*% a) #TODO: OPTIMIZE
          w_pred[, i] <- reg_Matrix %*% w_pred[, i]
          if (inner == n_inner & homeomorphic) w_pred[, i] <- make_homeo(w_pred[, i], twp)
        }
      }
      beta <- beta * kappa
    }
  }
  return(w_pred)
}
