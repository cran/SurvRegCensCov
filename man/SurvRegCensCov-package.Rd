\name{SurvRegCensCov-package}
\alias{SurvRegCensCov.package}
\alias{SurvRegCensCov}
\docType{package}
\title{Weibull Regression for a Right-Censored Endpoint with Interval-Censored Covariates}
\description{The function \code{\link{SurvRegCens}} of this package allows estimation of a Weibull Regression for a right-censored endpoint, one interval-censored covariate, and an arbitrary number of non-censored covariates. Additional functions allow to switch between different parametrizations of Weibull regression used by different \code{R} functions (\code{\link{ConvertWeibull}}, \code{\link{WeibullReg}}, \code{\link{WeibullDiag}}), inference for the mean difference of two arbitrarily censored Normal samples (\code{\link{NormalMeanDiffCens}}), and estimation of canonical parameters from censored samples for several distributional assumptions (\code{\link{ParamSampleCens}}).
}
\details{
\tabular{ll}{
Package: \tab SurvRegCensCov\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab 2022-01-12\cr
License: \tab GPL (>=2) \cr
}
}
\author{
Stanislas Hubeaux (maintainer), \email{stan.hubeaux@bluewin.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}

We thank Sarah Haile for contributing the functions \code{\link{ConvertWeibull}}, \code{\link{WeibullReg}}, \code{\link{WeibullDiag}} to the package.
}

\references{
Hubeaux, S. (2013). 
Estimation from left- and/or interval-censored samples. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.     

Hubeaux, S. (2013). 
Parametric Surival Regression Model with left- and/or interval-censored covariate. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.	

Hubeaux, S. and Rufibach, K. (2014).
\bold{SurvRegCensCov}: Weibull Regression for a Right-Censored Endpoint with a Censored Covariate.
Preprint, \url{https://arxiv.org/abs/1402.0432}.

Lynn, H. S. (2001). 
Maximum likelihood inference for left-censored HIV RNA data. 
\emph{Stat. Med.}, \bold{20}, 33--45.
	
Sattar, A., Sinha, S. K. and Morris, N. J. (2012). 
A Parametric Survival Model When a Covariate is Subject to Left-Censoring. 
\emph{Biometrics & Biostatistics}, \bold{S3}(2).
}

\keyword{package}

\examples{
# The main functions in this package are illustrated in their respective help files.
}
