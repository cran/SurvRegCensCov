\name{NormalMeanDiffCens}
\alias{NormalMeanDiffCens}
\title{Maximum Likelihood Estimator for the mean difference between two censored normally distributed samples}
\description{Computes estimates of the parameters of two censored Normal samples, as well as the mean difference between the two samples.}
\usage{NormalMeanDiffCens(censdata1, censdata2, conf.level = 0.95, 
     null.values = c(0, 0, 1, 1))}
\arguments{
   \item{censdata1}{Observations of first sample, format as specified by \code{code = interval2} in \code{\link{Surv}}.}
   \item{censdata2}{Observations of second sample, as specified by \code{code = interval2} in \code{\link{Surv}}.}
   \item{conf.level}{Confidence level for confidence intervals.}  
   \item{null.values}{Fixed values for hypothesis tests. Four-dimensional vector specifying the hypothesis for \eqn{\mu_1}, \eqn{\mu_2}, \eqn{\sigma_1}, \eqn{\sigma_2}.}
}
\value{A table with estimators and inference for the means and standard deviations of both samples, as well as the difference \eqn{\Delta} between the mean of the first and second sample. Hypothesis tests are for the values in \code{null.values} and for the null hypothesis of no mean difference.}

\references{
Hubeaux, S. (2013). 
Estimation from left- and/or interval-censored samples. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.  
	
Lynn, H. S. (2001). 
Maximum likelihood inference for left-censored HIV RNA data. 
\emph{Stat. Med.}, \bold{20}, 33--45.
}

\author{
Stanislas Hubeaux, \email{stan.hubeaux@bluewin.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}
}

\examples{
## example with interval-censored Normal samples
n <- 500
prop.cens <- 0.35
mu <- c(0, 2)
sigma <- c(1, 1)

set.seed(2013)

## Sample 1:
LOD1 <- qnorm(prop.cens, mean = mu[1], sd = sigma[1])
x1 <- rnorm(n, mean = mu[1], sd = sigma[1])
s1 <- censorContVar(x1, LLOD = LOD1)

## Sample 2:
LOD2 <- qnorm(0.35, mean = mu[2], sd = sigma[2])
x2 <- rnorm(n, mean = mu[2], sd = sigma[2])
s2 <- censorContVar(x2, LLOD = LOD2)

## inference on distribution parameters and mean difference:
NormalMeanDiffCens(censdata1 = s1, censdata2 = s2)
}
\keyword{methods}
