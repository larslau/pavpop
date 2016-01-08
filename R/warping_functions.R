#' Generate a warping function
#'
#' This function generates different types of warping functions that can be used with
#' pavpop.
#'
#' The possible types of warping functions are
#' \describe{
#'   \item{\code{shift}}{a horizontal shift.}
#'   \item{\code{linear}}{a linear streach of the curves.}
#'   \item{\code{pieceiwse-linear}}{a piecewise linear deformation of the curves.
#'                                  Anchor points are given by \code{tw}.}
#'   \item{\code{smooth}}{a smooth nonlinear deformation of the curves,
#'                                  produced by an increasing spline (Hyman filtered).
#'                                  Anchor points are given by \code{tw}.}
#' }
#'
#'
#'
#' @param type type of warping function. See 'Details'.
#' @param tw anchor points for the warping parameters. Needed for types 'piecewise-linear' and 'smooth'. \bold{Note:} boundary points must be included.
#' @param control list of control arguments. Entries \code{wleft} and \code{wright} control
#' whether the warping function is fixed at the boundary (\code{'fixed'}) or wheter
#' extrapolation from the last interior knot should be done (\code{'extrapolate'}).
#' @keywords warping
#' @export

make_warp_fct <- function(type = c('shift', 'linear', 'piecewise-linear', 'smooth'), tw = NULL, control = list(wleft = 'fixed', wright = 'fixed')) {
  # Match type argument
  if (length(type) > 1) type <- type[1]
  types <- c('shift', 'linear', 'piecewise-linear', 'smooth')
  type <- types[pmatch(type, types)]
  if (is.null(tw)) tw <- NA

  # If boundary control is missing, assume fixed
  if (is.null(control$wleft)) control$wleft <- 'fixed'
  if (is.null(control$wright)) control$wright <- 'fixed'


  if (type == 'shift') {
    v <- function(w, t, w_grad = FALSE) {
      if (!w_grad) {
        return(w + t)
      } else {
        return(matrix(1, length(t), 1))
      }
    }
    attr(v, 'mw') <- 1
  } else if (type == 'linear') {
    v <- function(w, t, w_grad = FALSE) {
      if (!w_grad) {
        return(w[1] + (w[2] + 1) * t)
      } else {
        dv <- matrix(1, length(t), 2)
        dv[, 2] <- t
        return(dv)
      }
    }
    attr(v, 'mw') <- 2
  } else if (type == 'piecewise-linear') {
    if (any(is.na(tw))) stop('all anchor points tw should be supplied for type \'piecewise-linear\'')
    mw <- length(tw)
    v <- function(w, t, w_grad = FALSE) {
      if (!w_grad) {
        vt <- t + approx(tw, c(ifelse(control$wleft == 'fixed', 0, w[1]),
                               w,
                               ifelse(control$wright == 'fixed', 0, w[length(w)])),
                         xout = t, rule = 2)$y
        return(vt)
      } else {
        # Derivative of warp function
        # Note: does not depend on w
        dv <- apply(cbind(tw[1:(mw - 2)], tw[2:(mw - 1)], tw[3:mw]), 1, function(x) {
          a <- rep(0, length(t))
          a[t > x[1] & t < x[2]] <- ((t - x[1]) / (x[2] - x[1]))[t > x[1] & t < x[2]]
          a[t >= x[2] & t < x[3]] <- (1 - ((t - x[2]) / (x[3] - x[2])))[t >= x[2] & t < x[3]]
          return(a)
        })
        if (control$wleft != 'fixed') dv[t < tw[2], 1] <- 1
        if (control$wright != 'fixed') dv[t > tw[mw - 1], mw - 2] <- 1
        return(dv)
      }
    }
    attr(v, 'mw') <- mw - 2
  } else if (type == 'smooth') {
    if (any(is.na(tw))) stop('all anchor points tw should be supplied for type \'smooth\'')
    mw <- length(tw)
    # Hyman spline warping function
    v <- function(w, t, w_grad = FALSE) {
      y <- tw + c(ifelse(control$wleft == 'fixed', 0, w[1]), w, ifelse(control$wright == 'fixed', 0, w[length(w)]))
      if (!all(diff(y) > 0)) {
        w <- make_homeo(w, tw, epsilon = 0.1)
        y <- tw + c(ifelse(control$wleft == 'fixed', 0, w[1]), w, ifelse(control$wright == 'fixed', 0, w[length(w)]))
      }
      if (!w_grad){
        return(spline(tw, y, xout = t, method = 'hyman')$y)
      } else {
        # Derivative of warp function
        # Finite difference, could we do better?
        epsilon <- 1e-5
        m <- length(t)
        dv <- matrix(0, m, mw - 2)
        for (j in 1:(mw - 2)) {
          h_tmp <- rep(0, mw)
          if (j == 1 & control$wleft != 'fixed') h_tmp[j] <- epsilon
          if (j == (mw - 2) & control$wright != 'fixed') h_tmp[j + 2] <- epsilon
          h_tmp[j + 1] <- epsilon
          dv[, j] <- (spline(tw, y + h_tmp, xout = t, method = 'hyman')$y
                      - spline(tw, y - h_tmp, xout = t, method = 'hyman')$y) / (2 * epsilon)
        }
        return(dv)
      }
    }
    attr(v, 'mw') <- mw - 2
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
  mw <- length(w)
  ui <- diag(1, nrow = mw)
  ui <- rbind(ui, rep(0, mw))
  ui[cbind(2:(mw + 1), 1:mw)] <- -1

  ci <- (epsilon - 1) * diff(tw)
  res <- constrOptim_inf(theta = rep(0, mw), f = function(x) mean((x - w)^2),
                         grad = function(x) 2 / mw * (x - w), ui = ui, ci = ci, method = "BFGS")

  return(res$par)
}
