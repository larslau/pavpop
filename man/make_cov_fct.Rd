% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariance_functions.R
\name{make_cov_fct}
\alias{make_cov_fct}
\title{Generate covariance matrix function from covariance function}
\usage{
make_cov_fct(cov_fct, noise = TRUE, param = NULL, ns = NULL, ...)
}
\arguments{
\item{cov_fct}{covariance function.}

\item{noise}{logical. Should an identity matrix be added to the covariance matrix?}

\item{param}{standard values of parameters}

\item{ns}{list of arguments to make a stationary covariance locally adaptive.
The list has takes the following entries \code{knots} (mandatory) the number of parameters
to control the non-stationary, \code{fixed} (optional) vector of fixed knots with reference
weight of 1, if not supplied the last knot is chosen. See examples and Details.}

\item{...}{arguments passed to cov_fct.}
}
\description{
Function that takes a covariance function
and returns a function that generates covariance matrices
according to the given covariance function.
}
\examples{
# Generate observation points
t <- seq(0, 1, length = 3)

# Generate covariances
matern_cov <- make_cov_fct(Matern, noise = FALSE)
bm_cov <- make_cov_fct(Brownian, noise = FALSE)
bb_cov <- make_cov_fct(Brownian, noise = FALSE, type = 'bridge')

# Evaluate covariance matrices
matern_cov(t, param = c(1, 1, 1))
bm_cov(t, param = c(1, 1, 1))
bb_cov(t, param = c(1, 1, 1))

# Plot covariance matrices
t <- seq(0, 1, length = 30)

persp(t, t, matern_cov(t, param = c(1, 0.5, 1)),
      theta = -30, phi = 30, ticktype = 'd', zlab = '',
      col = 'lightblue', shade = 0.2,
      main = 'Matern covariance (1, 0.5, 1)', zlim = c(0, 1))

persp(t, t, matern_cov(t, param = c(1, 0.5, 2)),
      theta = -30, phi = 30, ticktype = 'd', zlab = '',
      col = 'lightblue', shade = 0.2,
      main = 'Matern covariance (1, 0.5, 2)', zlim = c(0, 1))

persp(t, t, matern_cov(t, param = c(1, 0.1, 2)),
      theta = -30, phi = 30, ticktype = 'd', zlab = '',
      col = 'lightblue', shade = 0.2,
      main = 'Matern covariance (1, 0.1, 2)', zlim = c(0, 1))

persp(t, t, bm_cov(t, param = 1), theta = -30, phi = 30,
ticktype = 'd', zlab = '', col = 'lightblue', shade = 0.2,
main = 'Brownian motion covariance (tau = 1)', zlim = c(0, 1))

persp(t, t, bb_cov(t, param = 2), theta = -30, phi = 30,
ticktype = 'd', zlab = '', col = 'lightblue', shade = 0.2,
main = 'Brownian bridge covariance (tau = 2)', zlim = c(0, 1))

# Make covariance non-stationary
matern_cov_ns <- make_cov_fct(Matern, noise = FALSE,
                              ns = list(knots = 3, range = c(0, 1)))
# Use mid-point reference instead of last
matern_cov_ns_mid <- make_cov_fct(Matern, noise = FALSE,
                              ns = list(knots = 3, fixed = 2,
                                        range = c(0, 1)))

# Original covariance
persp(t, t, matern_cov_ns(t, param = c(1, 1, 1, 0.3, 2)),
      theta = -30, phi = 30, ticktype = 'd', zlab = '',
      col = 'lightblue', shade = 0.2,
      main = 'Matern covariance (1, 0.1, 2)', zlim = c(0, 1))

# Modified covariances
persp(t, t, matern_cov_ns(t, param = c(0.5, 0.7, 1, 0.3, 2)),
      theta = -30, phi = 30, ticktype = 'd', zlab = '',
      col = 'lightblue', shade = 0.2,
      main = 'Non-stationary Matern covariance (1, 0.1, 2)', zlim = c(0, 1.1))

persp(t, t, matern_cov_ns(t, param = c(1, 0.7, 1, 0.3, 2)),
      theta = -30, phi = 30, ticktype = 'd', zlab = '',
      col = 'lightblue', shade = 0.2,
      main = 'Non-stationary Matern covariance (1, 0.1, 2)', zlim = c(0, 1.1))

persp(t, t, matern_cov_ns_mid(t, param = c(1, 0.7, 1, 0.3, 2)),
      theta = -30, phi = 30, ticktype = 'd', zlab = '',
      col = 'lightblue', shade = 0.2,
      main = 'Non-stationary Matern covariance (1, 0.1, 2)', zlim = c(0, 1.1))
}

