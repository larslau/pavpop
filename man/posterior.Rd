% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood_functions.R
\name{posterior}
\alias{posterior}
\title{Posterior of the data given the random warping parameters}
\usage{
posterior(w, warp_fct, t, y, basis_fct, c, Sinv, Cinv)
}
\arguments{
\item{w}{warping parameters.}

\item{warp_fct}{warping function.}

\item{t}{evaluation points.}

\item{y}{values at evaluation points.}

\item{basis_fct}{basis function to describe the mean function.}

\item{c}{basis coefficients.}

\item{Sinv}{precision matrix for amplitude variation.}

\item{Cinv}{precision matrix for the warping parameters.}
}
\description{
This function calculates the posterior of the data given the random warping parameters
}
\keyword{posterior}
\keyword{warping}

