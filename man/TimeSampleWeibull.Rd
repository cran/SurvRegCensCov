\name{TimeSampleWeibull}
\alias{TimeSampleWeibull}
\alias{survtime}
\alias{F_time}
\title{Generate time-to-event data according to a Weibull regression model}
\description{Generates time-to-event data using the transform inverse sampling method, and such that the time-to-event is distributed according to a Weibull distribution induced by censored and/or non-censored covariates. Can be used to set up simulations.}
\usage{
TimeSampleWeibull(covariate_noncens = NULL, covariate_cens, lambda, gamma, beta)
}
\arguments{
  \item{covariate_noncens}{Censored covariate vector.}
  \item{covariate_cens}{Matrix where each column represents a non-censored covariate.}
  \item{lambda}{Scale parameter.}
  \item{gamma}{Shape parameter.}
  \item{beta}{Regression parameters, ordered as \eqn{\beta} linked to the non-censored covariates, \eqn{\beta} linked to the censored covariate.}
}

\note{The use of this function is illustrated in \code{\link{SurvRegCens}}.}

\author{Stanislas Hubeaux, \email{stanislas.hubeaux@epfl.ch}} 

\keyword{methods}
