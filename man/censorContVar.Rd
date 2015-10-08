\name{censorContVar}
\alias{censorContVar}
\title{Censor a vector of continuous numbers}
\description{Given a vector of realizations of a continuous random variable, interval-, left-, or right-censor these numbers at given boundaries. Useful when setting up simulations involving censored observations.}
\usage{censorContVar(x, LLOD = NA, ULOD = NA)}
\arguments{
  \item{x}{Vector of random numbers.}
  \item{LLOD}{Lower limit where \code{x} should be left-censored at. If no left-censoring, set to \code{NA}.}
  \item{ULOD}{Upper limit where \code{x} should be left-censored at. If no left-censoring, set to \code{NA}.}
}
\value{A \code{data.frame} as specified by \code{code = interval2} in \code{\link{Surv}}.}

\author{
Stanislas Hubeaux, \email{stanislas.hubeaux@alumni.epfl.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}
}

\examples{

## random vector
x <- rnorm(200)

## interval-censor this vector at -1 and 0.5
censorContVar(x, -1, 0.5)
}

\keyword{methods}

