\name{ParamSampleCens}
\alias{ParamSampleCens}
\title{Maximum Likelihood Estimator of parameters from a censored sample}
\description{Computes maximum likelihood estimators of the canonical parameters for several distributions, based on a censored sample.}
\usage{ParamSampleCens(censdata, dist = c("normal", "logistic", "gamma", "weibull")[1], 
     null.values = c(0, 1), conf.level = 0.95, initial = NULL)}
\arguments{
  \item{censdata}{Dataframe that contains censored data, format as specified by \code{code = interval2} in \code{\link{Surv}}.}
  \item{dist}{Assumed distribution of the sample.}
  \item{null.values}{Fixed values for hypothesis tests.}
  \item{conf.level}{Confidence level of confidence intervals.}  
  \item{initial}{Initial values for the maximization.}
}

\value{
	\item{coeff}{Estimators, standard errors, confidence intervals, and 2-sided \eqn{p}-values for the null hypothesis as given in \code{null.values}. All this inference is based on maximum likelihood theory, i.e. standard errors are computed using the inverse expected Hessian at the maximum likelihood estimator.}
	\item{percent.cens}{Percentage of censored observations.}
	\item{loglik}{Log likelihood function value at the estimator.}
	\item{info.converg}{Convergence information provided by the function \code{optim}.}
	\item{info.converg.message}{Message provided by the function \code{optim}.}
}

\references{
Hubeaux, S. (2013). 
Estimation from left- and/or interval-censored samples. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.  

Lynn, H. S. (2001). 
Maximum likelihood inference for left-censored HIV RNA data. 
\emph{Stat. Med.}, \bold{20}, 33--45.
}

\note{Functions with similar functionality are provided in the package \pkg{fitdistrplus}.}

\author{
Stanislas Hubeaux, \email{stanislas.hubeaux@epfl.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}
}

\examples{
n <- 500
prop.cens <- 0.35

## example with a left-censored Normally distributed sample
set.seed(2013)
mu <- 3.5
sigma <- 1
LOD <- qnorm(prop.cens, mean = mu, sd = sigma)
x1 <- rnorm(n, mean = mu, sd = sigma)
s1 <- censorContVar(x1, LLOD = LOD)
ParamSampleCens(censdata = s1)


## example with an interval-censored Normal sample
set.seed(2013)
x2 <- rnorm(n, mean = mu, sd = sigma)
LOD <- qnorm(prop.cens / 2, mean = mu, sd = sigma)
UOD <- qnorm(1 - prop.cens / 2, mean = mu, sd = sigma)
s2 <- censorContVar(x2, LLOD = LOD, ULOD = UOD)
ParamSampleCens(censdata = s2)


\dontrun{
## compare to fitdistrplus
library(fitdistrplus)
s2 <- as.data.frame(s2)
colnames(s2) <- c("left", "right")
summary(fitdistcens(censdata = s2, distr = "norm"))
}
}
\keyword{methods}
