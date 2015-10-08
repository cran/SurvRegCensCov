\name{LoglikCens}
\alias{LoglikNormalCens}
\alias{LoglikWeibullCens}
\alias{LoglikLogisticCens}
\alias{LoglikGammaCens}
\title{Log-likelihood functions for estimation of canonical parameters from a censored sample}
\description{Computes the log-likelihood function for a censored sample, according to a specified distributional assumptions. Available distributions are Normal, Weibull, Logistic, and Gamma.}
\usage{
LoglikNormalCens(x, data, lowerbound, vdelta)
LoglikWeibullCens(x, data, lowerbound, vdelta)
LoglikLogisticCens(x, data, lowerbound, vdelta)
LoglikGammaCens(x, data, lowerbound, vdelta)
}
\arguments{
  \item{x}{Two-dimensional vector giving the canonical parameters of the distribution.}
  \item{data}{Observed or censored event times.}
  \item{lowerbound}{A vector that collect lower bounds for the interval-censored observations. If no lower bound is available then put \code{NA}.}
  \item{vdelta}{A vector which indicates censoring (0: censored, 1: not censored).}
}

\references{
Hubeaux, S. (2013). 
Estimation from left- and/or interval-censored samples. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.  

Lynn, H. S. (2001). 
Maximum likelihood inference for left-censored HIV RNA data. 
\emph{Stat. Med.}, \bold{20}, 33--45.
}

\author{
Stanislas Hubeaux, \email{stanislas.hubeaux@alumni.epfl.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}
}

\note{Function not intended to be invoked by the user.}

\keyword{methods}
