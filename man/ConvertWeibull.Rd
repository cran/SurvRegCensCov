\name{ConvertWeibull}
\alias{ConvertWeibull}
\title{Transformation of survreg output for the Weibull distribution}
\description{
Transforms output from \code{\link[survival]{survreg}} using the Weibull distribution to a more natural parameterization. See details and the vignette for more information.
}
\usage{ConvertWeibull(model, conf.level = 0.95)}
\arguments{
  \item{model}{A \code{\link[survival]{survreg}} model, with \code{dist = "weibull"}.}
  \item{conf.level}{Confidence level used to produce two-sided \eqn{1-\alpha/2} confidence intervals for the hazard and event time ratios.}
}
\details{    
The \code{\link[survival]{survreg}} function fits a Weibull accelerated failure time model of the form 

\deqn{\log t = \mu + \alpha^T Z + \sigma W,}{log t = \mu + \alpha^T Z + \sigma W,} 

where \eqn{Z} is a matrix of covariates, and \eqn{W} has the extreme value distribution, \eqn{\mu} is the intercept, 
\eqn{\alpha} is a vector of parameters for each of the covariates, and \eqn{\sigma} is the scale. The usual 
parameterization of the model, however, is defined by hazard function 

\deqn{h(t|Z) = \gamma \lambda t^{\gamma - 1} \exp(\beta^T Z).}{h(t|Z) = \gamma \lambda t^{\gamma - 1} exp(\beta^T Z).}

The transformation is as follows: \eqn{\gamma = 1/\sigma}, \eqn{\lambda = \exp(-\mu/\sigma)}{\lambda = exp(-\mu/\sigma)}, and 
\eqn{\beta=-\alpha/\sigma}, and estimates of the standard errors can be found using the delta method.

The Weibull distribution has the advantage of having two separate interpretations. The first, via proportional hazards, leads to a hazard ratio, defined by \eqn{\exp \beta}{exp \beta}. The second, of accelerated failure times, leads to an event time ratio (also known as an acceleration factor), defined by \eqn{\exp (-\beta/\gamma)}{exp (-\beta/\gamma)}.

Further details regarding the transformations of the parameters and their standard errors can be found in Klein and 
Moeschberger (2003, Chapter 12). An explanation of event time ratios for the accelerated failure time interpretation of the model can be found in Carroll (2003). A general overview can be found in the \code{vignette("weibull")} of this package.
}
\value{
    \item{vars}{A matrix containing the values of the transformed parameters and their standard errors}
    \item{HR}{A matrix containing the hazard ratios for the covariates, and \eqn{1 - }\code{conf.level}\eqn{ / 2} confidence intervals.}
    \item{ETR}{A matrix containing the event time ratios for the covariates, and \eqn{1 - }\code{conf.level}\eqn{ / 2} confidence intervals.} 
}
\references{
Carroll, K. (2003). 
On the use and utility of the Weibull model in the analysis of survival data.
\emph{Controlled Clinical Trials}, \bold{24}, 682--701.

Klein, J. and Moeschberger, M. (2003). 
\emph{Survival analysis: techniques for censored and truncated data}. 
2nd edition, Springer.
}

\author{Sarah R. Haile, Epidemiology, Biostatistics and Prevention Institute (EBPI), University of Zurich, \email{sarah.haile@uzh.ch}}
 
\seealso{This function is used by \code{\link{WeibullReg}}.}
\examples{
data(larynx)
ConvertWeibull(survreg(Surv(time, death) ~ stage + age, larynx), conf.level = 0.95)
}
\keyword{survival}
\keyword{regression}
