\name{LoglikWeibullSurvRegCens}
\alias{LoglikWeibullSurvRegCens}
\alias{h_conditional}
\alias{S_conditional}
\title{Log-likelihood function of a Weibull Survival Regression Model allowing for an interval-censored covariate.}
\description{Computes the log-likelihood function of a Weibull Survival Regression Model allowing for an interval-censored covariate.}
\usage{
LoglikWeibullSurvRegCens(x, data_y, data_delta_loglik, data_cov_noncens = NULL, 
                         data_cov_cens, density, data_r_loglik, data_lowerbound,
                         intlimit = 10^-10)
}
\arguments{
  \item{x}{Vector of parameters, ordered as follows: Scale parameter, Shape parameter, regression parameters (i.e. \eqn{\beta}) linked to the non-censored covariates, regression parameter (i.e. \eqn{\beta}) linked to the censored covariate.}
  \item{data_y}{Time-to-event vector.}
  \item{data_delta_loglik}{Censored indicator vector of the time-to-event (0: censored, 1: not censored).}
  \item{data_cov_noncens}{Matrix where each column represents a non-censored covariate.}
  \item{data_cov_cens}{Censored covariate vector.}
  \item{density}{Density function of the censored covariate.}
  \item{data_r_loglik}{Censored indicator vector of the censored covariate (0: censored, 1: not censored).}
  \item{data_lowerbound}{A vector which corresponds to the lower bounds for the interval-censored observations of the censored covariate. If no lower bound is available then put \code{NA}.}
    \item{intlimit}{In computation of integrals, values of the function to be integrated below \code{intlimit} are set to 0. This makes integration results more accurate and speeds up integration. If the data is such that the absolute values of the underlying baseline Weibull density are very small, i.e. in the range of \code{intlimit}, it is advisable to rescale the time variable, e.g. change the scaling from days to years. A very small value of the estimated \eqn{\lambda} is indicative of that situation.}
}

\references{
Hubeaux, S. (2013). 
Parametric Surival Regression Model with left- and/or interval-censored covariate. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.  
	
Sattar, A., Sinha, S. K. and Morris, N. J. (2012). 
A Parametric Survival Model When a Covariate is Subject to Left-Censoring. 
\emph{Biometrics & Biostatistics}, \bold{S3}(2).
}

\author{
Stanislas Hubeaux, \email{stan.hubeaux@bluewin.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}
}

\note{Function not intended to be invoked by the user.}

\keyword{methods}
