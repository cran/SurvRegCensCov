\name{logLik.src}
\alias{logLik.src}
\title{Extract value of log-likelihood at maximum for Weibull regression with an interval-censored covariate}
\description{\code{logLik} method for class \code{"src"}.}
\usage{\method{logLik}{src}(object, ...)}
\arguments{
  \item{object}{An object of class \code{"src"}, usually a result of a call to \code{SurvRegCens}.}
  \item{\dots}{Further arguments.}
}
\value{The function \code{logLik.src} returns the value of the log-likelihood at the maximum likelihood estimate, as well as the corresponding degrees of freedom.}
\references{
Hubeaux, S. (2013). 
Parametric Surival Regression Model with left- and/or interval-censored covariate. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.     

Hubeaux, S. and Rufibach, K. (2014).
\bold{SurvRegCensCov}: Weibull Regression for a Right-Censored Endpoint with a Censored Covariate.
Preprint, \url{http://arxiv.org/abs/1402.0432}.

Sattar, A., Sinha, S. K. and Morris, N. J. (2012). 
A Parametric Survival Model When a Covariate is Subject to Left-Censoring. 
\emph{Biometrics & Biostatistics}, \bold{S3}(2).
}
\author{
Stanislas Hubeaux, \email{stanislas.hubeaux@epfl.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}
}
\examples{
## See help file of function "SurvRegCens".
}
\keyword{htest}
\keyword{nonparametric}
