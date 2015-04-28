#' Evaluate warping function given weights
#'
#' This function evaluates the warping function corresponding to a specified set of warping parameters at a set of evaluation points using linear interpolation.
#' @param w warping parameters.
#' @param t evaluation points.
#' @param tw anchor points for the warping parameters.
#' @param smooth logical. Should warping function be based on a monotonic cubic spline?
#' @keywords warping
#' @export
#' @examples
#' t <- seq(0, 1, length = 10)
#' plot(t, v(0.4, t, 0.5), type = 'l', pch = 19)

v <- function(w, t, tw, smooth = FALSE) {
  if (smooth) {
    # Make possible argument
    x <- c(0, tw, 1)
    y <- c(0, tw + w, 1)
    if (!all(diff(y) > 0)) {
      y <- c(0, tw + make_homeo(w, tw, epsilon = 0.2), 1)
    }
    return(spline(x, y, xout = t, method = 'hyman')$y)
  } else {
    return(t + approx(c(0, tw, 1), c(0, w, 0), xout = t, rule = 2)$y)
  }
}

#' Derivative of piecewise linear warping function warping
#'
#' Evaluate the derivative of a piecewise linear warping function
#' @param t evaluation points.
#' @param tw anchor points for the warping parameters.
#' @keywords warp derivative
#' @export

dv <- function(t, tw) {
  mw <- length(tw)
  # Derivative of warp function
  apply(cbind(c(0, tw[-mw]), tw, c(tw[-1], 1)), 1, function(x) {
    a <- rep(0, length(t))
    a[t > x[1] & t < x[2]] <- ((t - x[1])/(x[2] - x[1]))[t > x[1] & t < x[2]]
    a[t >= x[2] & t < x[3]] <- (1 - ((t - x[2])/(x[3] - x[2])))[t >= x[2] & t < x[3]]
    return(a)
  })
}

#' Compute the inverse warping function weights
#'
#' This function evaluates the inverse warping function corresponding to a specified set of warping parameters at a set of evaluation points using linear interpolation.
#' @param w warping parameters.
#' @param t evaluation points.
#' @param tw anchor points for the warping parameters.
#' @keywords warping
#' @export
#' @examples
#' t <- seq(0, 1, length = 10)
#' plot(vinv(0.4, t, 0.5), type = 'l', pch = 19)

vinv <- function(w, t, tw) {
  approx(c(0, tw + w, 1), c(0, tw, 1), xout = t, rule = 2)$y
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
  ci <- rep(-1 / (nw + 1), nw + 1)
  epsilon <- epsilon / (nw + 1)
  res <- constrOptim(theta = rep(0, nw), f = function(x) mean((x - w)^2),
                     grad = function(x) 2 * (x - w), ui = ui, ci = (ci + epsilon), method = "BFGS")

  return(res$par)
}
