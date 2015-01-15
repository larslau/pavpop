#' Jacobian of warped B-spline function
#'
#' Computes the Jacobian of a warped B-spline function.
#' @param t (warped) evaluation points.
#' @param dvt Jacobian of the warping function for the given evaluation points.
#' @param kts B-spline knots.
#' @export

Zi <- function(t, dwarp, c, kts) {
  basis <- bsd(t, knots = kts, Boundary.knots = c(0, 1))
  if (ncol(basis) + 1 == length(c)) c <- c[-1]
  dwarp <- dwarp * ((basis %*% c)[, 1])
  return(dwarp)
}


# TODO: CLEAN bsd AND DOCUMENT!

#' Spline basis derivative
#'
#' Spline basis derivative.
#' @param x the predictor variable.  Missing values are allowed.
#' @param df
#' @param knots
#' @param degree
#' @param intercept
#' @param Boundary.knots
#' @export

bsd <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
    Boundary.knots = range(x)) {
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax))
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
    } else outside <- FALSE
    ord <- 1L + (degree <- as.integer(degree))
    if (ord <= 1)
        stop("'degree' must be integer >= 1")
    if (!is.null(df) && is.null(knots)) {
        nIknots <- df - ord + (1L - intercept)
        if (nIknots < 0L) {
            nIknots <- 0L
            warning(gettextf("'df' was too small; have used %d", ord - (1L -
                intercept)), domain = NA)
        }
        knots <- if (nIknots > 0L) {
            knots <- seq.int(from = 0, to = 1, length.out = nIknots + 2L)[-c(1L,
                nIknots + 2L)]
            stats::quantile(x[!outside], knots)
        }
    }
    Aknots <- sort(c(rep(Boundary.knots, ord), knots))
    if (any(outside)) {
        warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
        derivs <- 0:degree
        scalef <- gamma(1L:ord)
        basis <- array(0, c(length(x), length(Aknots) - degree - 1L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree, "^"))
            tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, derivs)
            basis[ol, ] <- xl %*% (tt/scalef)
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree, "^"))
            tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, derivs)
            basis[or, ] <- xr %*% (tt/scalef)
        }
        if (any(inside <- !outside))
            basis[inside, ] <- splineDesign(Aknots, x[inside], ord, derivs = rep(1,
                length = length(x[inside])))
    } else basis <- splineDesign(Aknots, x, ord, derivs = rep(1, length = length(x)))
    if (!intercept)
        basis <- basis[, -1L, drop = FALSE]
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots,
        Boundary.knots = Boundary.knots, intercept = intercept)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("bs", "basis", "matrix")
    basis
}
