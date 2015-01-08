#TODO: CONSISTENTLY USE c(0, 1) as boundary knots!!


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

#TODO: Change n_par
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




#TODO: DELETE!!

#' Stable locally linearized likelihood function
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

#TODO: Change n_par

like_stable <- function(param, n_par, r, Zis, amp_cov_fct, warp_cov_fct, t, tw) {
  amp_cov_par <- param[1:n_par[1]]
  warp_cov_par <- param[(n_par[1] + 1):length(param)]

  C <- warp_cov_fct(tw, warp_cov_par)
  Cinv <- chol2inv(chol(C))

  n <- length(r)
  m <- sapply(r, length)

  sq <- logdet <- 0
  for (i in 1:n) {
    ZZ <- Zis[[i]]
    S <- amp_cov_fct(t[[i]], amp_cov_par)
    V <- as.matrix(S + ZZ %*% Cinv %*% Matrix::t(ZZ))
    rr <- r[[i]]
    sq <- sq + as.numeric(t(rr) %*% solve(V, rr))
    logdet <- logdet + determinant(V)$modulus
  }
  sigmahat <- 1/sum(m) * as.numeric(sq)
  res <- sum(m) * log(sigmahat) + logdet

  return(res)
}



#' Noise scale estimate from locally linearized likelihood function
#'
#' Computes the noise scale estimate from the linearized likelihood given the residual around a predicted warp and the corresponding Jacobians.
#' @param param variance parameters.
#' @param r residual.
#' @param Zis list of Jacobians in the warps of the mean function around the given warp.
#' @param Cinv inverse warp covariance.
#' @keywords likelihood
#' @keywords linearization
#' @export
#' @examples
#' #TODO


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
#' @param Boundary.knots boundary knots for the B-spline basis.
#' @keywords warping
#' @keywords posterior
#' @export
#' @examples
#' TODO

#TODO: ROXYGEN USE PARAMETERS FROM LIKE

posterior <- function(w, t, y, tw, c, Ainv, Cinv, kts, Boundary.knots = c(0, 1)) {
  r <- y - bs(v(w, t, tw), knots = kts, Boundary.knots = Boundary.knots) %*% c
  return((t(r) %*% Ainv %*% r + t(w) %*% Cinv %*% w)[1])
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
#' @param Boundary.knots boundary knots for the B-spline basis.
#' @keywords warping
#' @keywords posterior
#' @export
#' @examples
#' TODO

#TODO: ROXYGEN USE PARAMETERS FROM LIKE

posterior_grad <- function(w, dwarp, t, y, tw, c, Ainv, Cinv, kts) {
  vt <- v(w, t, tw)
  r <- bs(vt, knots = kts, Boundary.knots = c(0, 1)) %*% c - y
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
#' @examples
#' TODO

# TODO: UPDATE, DON'T USE NOW!

predict_warp_pyramid <- function(w, y, Ainv, t, tw, kts, warp_cov_par, warp_cov_fct, plevels = 5, beta0 = 1, start = 1, homeomorphic = TRUE, iter = c(10, 5)) {
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
      c <- spline_weights(y, t, w_pred, twp, bdiag(Ainv), kts)
      for (inner in 1:n_inner) {
        for (i in 1:n) { #TODO: CHECK THIS CODE!!!
          t_warped <- v(w_pred[, i], t[[i]], twp)
          # The basis construction is also done in spline_weights
          wbasis <- bs(t_warped, knots = kts, Boundary.knots = c(0, 1))
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
