\name{coef.src}
\alias{coef.src}
\title{Extract coefficients of Weibull regression with an interval-censored covariate}
\description{\code{coef} method for class \code{"src"}.}
\usage{\method{coef}{src}(object, ...)}
\arguments{
  \item{object}{An object of class \code{"src"}, usually a result of a call to \code{SurvRegCens}.}
  \item{\dots}{Further arguments.}
}
\value{The function \code{coef.src} returns the estimated parameters of the Weibull regression when calling \code{\link{SurvRegCens}}.}
\references{
Hubeaux, S. (2013). 
Parametric Surival Regression Model with left- and/or interval-censored covariate. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.     

Hubeaux, S. and Rufibach, K. (2014).
\bold{SurvRegCensCov}: Weibull Regression for a Right-Censored Endpoint with a Censored Covariate.
Preprint, \url{https://arxiv.org/abs/1402.0432}.

Sattar, A., Sinha, S. K. and Morris, N. J. (2012). 
A Parametric Survival Model When a Covariate is Subject to Left-Censoring. 
\emph{Biometrics & Biostatistics}, \bold{S3}(2).
}
\author{
Stanislas Hubeaux, \email{stan.hubeaux@bluewin.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}
}
\examples{
## See help file of function "SurvRegCens".
}
\keyword{htest}
\keyword{nonparametric}
