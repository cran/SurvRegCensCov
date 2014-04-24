\name{summary.src}
\alias{summary.src}
\title{Summarizing Weibull regression with an interval-censored covariate}
\description{\code{summary} method for class \code{"src"}.}
\usage{\method{summary}{src}(object, ...)}
\arguments{
  \item{object}{An object of class \code{"src"}, usually a result of a call to \code{SurvRegCens}.}
  \item{\dots}{Further arguments.}
}
\value{The function \code{summary.src} returns the estimated parameters, incl. statistical inference, of the Weibull regression, incl. AIC, when calling \code{\link{SurvRegCens}}.}
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
