#' Generate a warping function
#'
#' This function generates a certain kind of warping function that can be used by pavpop.
#' @param type type of warping function.
#' @param tw anchor points for the warping parameters. Needed for types 'piecewise-linear' and 'smooth'.
#' @keywords warping
#' @export

make_warp_fct <- function(type = c('shift', 'linear', 'piecewise-linear', 'smooth'), tw = NULL) {
  # Match type argument
  types <- c('shift', 'linear', 'piecewise-linear', 'smooth')
  type <- types[pmatch(type, types)]
  if (is.null(tw)) tw <- NA

  if (type == 'shift') {
    v <- function(w, t, w_grad = FALSE) {
      if (!w_grad) {
        return(w + t)
      } else {
        return(rep(1, length = length(t)))
      }
    }
    attr(v, 'initialize') <- 0
    attr(v, 'mw') <- 1
  } else if (type == 'linear') {
    v <- function(w, t, w_grad = FALSE) {
      if (!w_grad) {
        return(w[1] + w[2] * t)
      } else {
        dv <- matrix(1, length(t), 2)
        dv[, 2] <- t
        return(dv)
      }
    }
    attr(v, 'initialize') <- c(0, 1)
    attr(v, 'mw') <- 2
  } else if (type == 'piecewise-linear') {
    if (any(is.na(tw))) stop('all anchor points tw should be supplied for type \'piecewise-linear\'')
    if (min(tw) < 0 | max(tw) > 1) stop('anchor points tw should be within the interval (0, 1)')
    mw <- length(tw)

    v <- function(w, t, w_grad = FALSE) {
      if (!w_grad) {
        vt <- t + approx(c(0, tw, 1), c(0, w, 0), xout = t, rule = 2)$y
        vt[vt < 0] <- 0
        vt[vt > 1] <- 1
        return(vt)
      } else {
        # Derivative of warp function
        # Note: does not depend on w
        dv <- apply(cbind(c(0, tw[-mw]), tw, c(tw[-1], 1)), 1, function(x) {
          a <- rep(0, length(t))
          a[t > x[1] & t < x[2]] <- ((t - x[1])/(x[2] - x[1]))[t > x[1] & t < x[2]]
          a[t >= x[2] & t < x[3]] <- (1 - ((t - x[2])/(x[3] - x[2])))[t >= x[2] & t < x[3]]
          return(a)
        })
        return(dv)
      }
    }
    attr(v, 'initialize') <- 0
    attr(v, 'mw') <- mw
  } else if (type == 'smooth') {
    if (any(is.na(tw))) stop('all anchor points tw should be supplied for type \'smooth\'')
    if (min(tw) < 0 | max(tw) > 1) stop('anchor points tw should be within the interval (0, 1)')
    mw <- length(tw)
    # Hyman spline warping function
    v <- function(w, t, w_grad = FALSE) {
      x <- c(0, tw, 1)
      y <- c(0, tw + w, 1)
      if (!all(diff(y) > 0)) {
        y <- c(0, tw + make_homeo(w, tw, epsilon = 0.1), 1)
      }
      if (!w_grad){
        return(spline(x, y, xout = t, method = 'hyman')$y)
      } else {
        # Derivative of warp function
        # Finite difference, could we do better?
        epsilon <- 1e-5
        m <- length(t)
        dv <- matrix(0, m, mw)
        for (j in 1:mw) {
          h_tmp <- rep(0, mw + 2)
          h_tmp[j + 1] <- epsilon
          dv[, j] <- (spline(x, y + h_tmp, xout = t, method = 'hyman')$y
                      - spline(x, y - h_tmp, xout = t, method = 'hyman')$y) / (2 * epsilon)
        }
        return(dv)
      }
    }
    attr(v, 'initialize') <- 0
    attr(v, 'mw') <- mw
  }
  attr(v, 'tw') <- tw
  attr(v, 'type') <- type
  return(v)
}


#' Project warping parameters onto the set of values that produce homeomorphic warping functions
#'
#' Projects a set of warping parameters onto the set of values that produce homeomorphic warping functions using constrOptim.
#' @param w warping parameters.
#' @param tw anchor points for the warping parameters.
#' @param epsilon the smallest allowed slope between two points at unit spacing.
#' @keywords warping
#' @export

make_homeo <- function(w, tw, epsilon = 0.1) {
  nw <- length(w)
  ui <- diag(1, nrow = nw)
  ui <- rbind(ui, rep(0, nw))
  ui[cbind(2:(nw + 1), 1:nw)] <- -1

  ci <- (epsilon - 1) * diff(c(0, tw, 1))
  res <- constrOptim(theta = rep(0, nw), f = function(x) mean((x - w)^2),
                     grad = function(x) 2 / nw * (x - w), ui = ui, ci = ci, method = "BFGS")

  return(res$par)
}
