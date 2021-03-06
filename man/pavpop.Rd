% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pavpop.R
\name{pavpop}
\alias{pavpop}
\title{Estimate parameters and predict warps for curve data}
\usage{
pavpop(y, t, basis_fct, warp_fct, amp_cov = NULL, warp_cov = NULL,
  amp_fct = NULL, warped_amp = FALSE, iter = c(5, 5),
  parallel = list(n_cores = 1, parallel_likelihood = FALSE),
  use_warp_gradient = FALSE, warp_optim_method = "CG",
  homeomorphisms = "no", like_optim_control = list())
}
\arguments{
\item{y}{list of \eqn{n} functional observations. Missing values are allowed.}

\item{t}{list of time points corresponding to y. Should be scaled to have outer endpoints at 0 and 1.}

\item{basis_fct}{basis function to describe the mean function.}

\item{warp_fct}{warp function that models the warping procedure.}

\item{amp_cov}{amplitude covariance matrix function. If `NULL`, the amplitude is assumed to only contain iid Gaussian noise.}

\item{warp_cov}{warp covariance matrix function. If `NULL` the warps are treated as fixed parameters.}

\item{amp_fct}{functional basis that amplitude variation should be expressed in. If used, pavpop will automatically assume that iid. Gaussian noise is present in the data. Default is \code{NULL}.}

\item{warped_amp}{logical. Does the amplitude correlation follow observed or warped time?}

\item{iter}{two-dimensional numeric consisting of number of outer and inner iterations.}

\item{parallel}{list containing elements \code{n_cores}, the number of cores to use when
predicting warps and \code{parallel_likelihood}, a logical indicating whether likelihood
optimization should use a parallel evaluation of the gradient. In order for the parallel
likelihood to give a speed-up, large data sizes are required.}

\item{use_warp_gradient}{logical. Should warp prediction use the exact gradient for based optimization?}

\item{warp_optim_method}{optimization method for predicting warp. Defaults to 'CG' which gives robust results. 'Nelder-Mead' is often much faster, but not as reliable.}

\item{homeomorphisms}{should warps be constrained to be homeomorphisms? Options are: \code{'no'}, \code{'soft'} or \code{'hard'}. 'soft' will project the prediction onto the space of homeomorphisms after each prediction. 'hard' will do the optimiziation in the constrained space (not implemented yet!).}

\item{like_optim_control}{list of control options for likelihood optimization. Parameters are given as \code{c(amp_cov_par, warp_cov_par)} and options include lower, upper, method, ndev (see \code{\link[stats]{optim}}).}
}
\description{
This function does likelihood estimation in the model \deqn{y_i(t)=\theta(v(t, w_i))+x_i(t)+\epsilon_i(t)} based on iterative local linearization of the model around predictions of the random warping parameters \eqn{w_i}.
}
\examples{
# Load male growth data from the Berkeley growth study
t <- fda::growth$age
y <- fda::growth$hgtm
m <- nrow(y)
n <- ncol(y)

# Specify age rage for controlling boundary points
t_range <- c(0, 20)
t <- replicate(n, t, simplify = FALSE)
y <- lapply(1:n, function (x) y[, x])

# Set up basis function
kts <- seq(t_range[1], t_range[2], length = 15)
basis_fct <- make_basis_fct(kts = kts, type = 'increasing', intercept = TRUE,
                            control = list(boundary = t_range))

# Set up warp function
tw <- seq(t_range[1], t_range[2], length = 6)
warp_fct <- make_warp_fct('smooth', tw, control = list(wright = 'extrapolate'))
mw <- attr(warp_fct, 'mw')

# Set up covariance functions
warp_cov_par <- c(tau = 10)
warp_cov <- make_cov_fct(Brownian, noise = FALSE, param = warp_cov_par, type = 'motion',
                         range = t_range)

amp_cov_par <- c(scale = 200, range = 10, smoothness = 2)
amp_cov <- make_cov_fct(Matern, noise = TRUE, param = amp_cov_par)


# Estimate in the model

# Bounds of parameters
# NOTE: Prediction of velocities is only meaningful
#       when the smoothness parameter is > 0.5
lower <- c(1e-2, 1e-2, 0.5001, 1e-2)
upper <- c(1000, Inf, Inf, Inf)

res <- pavpop(y, t, basis_fct, warp_fct, amp_cov, warp_cov, homeomorphisms = 'soft',
              like_optim_control = list(lower = lower, upper = upper))
#
# Plot results
#

t_p <- seq(range(t)[1], range(t)[2], length = 100)

# Functional fixed effect
theta <- basis_fct(t_p) \%*\% res$c

# Display data with predictions
plot(t_p, theta, ylim = range(y), type = 'n', main = 'Original heights and predicted',
     xlab = 'Age (years)', ylab = 'Height (cm)')
for (i in 1:n) {
  points(t[[i]], y[[i]], pch = 19, cex = 0.3, col = rainbow(n)[i])
  lines(t_p, predict_curve(t_p, t[[i]], y[[i]], basis_fct, res$c, warp_fct, res$w[, i],
                           amp_cov, res$amp_cov_par),
        lwd = 0.5, col = rainbow(n)[i])

}
lines(t_p, theta, ylim = range(y), lwd = 2, lty = 2)

# Compute and display growth velocities
plot(t_p, t_p, ylim = c(0, 23), type = 'n', main = 'Predicted growth velocities',
     xlab = 'Age (years)', ylab = 'Growth velocity (cm/year)')
for (i in 1:n) {
  lines(t_p, predict_curve(t_p, t[[i]], y[[i]], basis_fct, res$c, warp_fct, res$w[, i],
                           amp_cov, res$amp_cov_par, deriv = TRUE),
        lwd = 0.5, col = rainbow(n)[i])
}


# Display predicted warping functions
plot(t_p, t_p, type = 'l', lwd = 2, lty = 2, main = 'Warping functions',
     xlab = 'Age (years)', ylab = 'Biological age (years)')
for (i in 1:n) lines(t[[i]], warp_fct(res$w[,i], t[[i]]), lwd = 0.4, col = rainbow(n)[i])
}
\references{
L.L. Raket, S. Sommer, and B. Markussen, "A nonlinear mixed-effects model for simultaneous smoothing and registration of functional data," Pattern Recognition Letters, vol. 38, pp. 1-7, 2014, 10.1016/j.patrec.2013.10.018.
}
\seealso{
See the vignettes for many more examples.
}
\keyword{estimation}
\keyword{likelihood}

