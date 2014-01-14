\name{WeibullIntegrate}
\alias{WeibullIntegrate}
\title{Function to be integrated in function \code{\link{SurvRegCens}}}
\description{Function to be integrated to compute log-likelihood function for the Weibull survival regression model with a censored covariate.}
\usage{
WeibullIntegrate(x, x_i_noncens = NULL, density, param_y_i,
                 param_delta_i, param_lambda, param_gamma,
                 param_beta, intlimit = 10^-10, ForIntegrate = TRUE)
}
\arguments{
  \item{x}{Value of the censored covariate for observation \eqn{i}.}
  \item{x_i_noncens}{Vector of values of the non-censored covariates for observation \eqn{i}, i.e. one row of the matrix where each column represents a non-censored covariate.}
  \item{density}{Density function of the censored covariate.}
  \item{param_y_i}{Value of the time-to-event for observation \eqn{i}.}
  \item{param_delta_i}{Censoring indicator of time-to-event for observation \eqn{i} (0: censored, 1: not censored).}
  \item{param_lambda}{Scale parameter of the Weibull distribution.}
  \item{param_gamma}{Shape parameter of the Weibull distribution.}
  \item{param_beta}{Regression parameters (i.e. \eqn{\beta}): (betaNonCens1, ..., betaNonCens, betaCens)}
  \item{intlimit}{In computation of integrals, values of the function to be integrated below \code{intlimit} are set to 0. This makes integration results more accurate and speeds up integration. If the data is such that the absolute values of the underlying baseline Weibull density are very small, i.e. in the range of \code{intlimit}, it is advisable to rescale the time variable, e.g. change the scaling from days to years. A very small value of the estimated \eqn{\lambda} is indicative of that situation.}
  \item{ForIntegrate}{\code{logical} indicating whether the function is to be integrated or not.}
}

\note{Function is not intended to be invoked by the user.}

\author{Stanislas Hubeaux, \email{stanislas.hubeaux@epfl.ch}}

\keyword{methods}
