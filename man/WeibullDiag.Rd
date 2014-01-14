\name{WeibullDiag}
\alias{WeibullDiag}
\title{Diagnostic Plot of Adequacy of Weibull Distribution}
\description{This function constructs a diagnostic plot of the adequacy of the Weibull distribution for survival 
data with respect to one categorical covariate. If the Weibull distribution fits the data well, then the 
lines produced should be linear and parallel.
}
\usage{WeibullDiag(formula, data = parent.frame(), labels = names(m$strata))}
\arguments{
  \item{formula}{A formula containing a \code{\link{Surv}} object, should only contain one categorical 
  predictor, or a set of indicators describing only one predictor.}
  \item{data}{Data set.}
  \item{labels}{A vector containing labels for the plotted lines.}
}
\details{
As discussed in Klein and Moeschberger (2003), one method for checking the adequacy of the Weibull model with a 
categorical covariate is to produce stratified Kaplan-Meier estimates (KM), which can be transformed to estimate 
the log cumulative hazard for each stratum. Then in a plot of \eqn{\log(t)}{log(t)} versus \eqn{\log(-\log(KM))}{log(-log(KM))}, the 
lines should be linear and parallel. This can be seen as the log cumulative hazard for the Weibull distribution 
is 

\deqn{\log H(t) = \log \lambda + \alpha \log t.}{log H(t) = log \lambda + \alpha log t.}
}
\value{Produces a plot of log Time vs. log Estimated Cumulative Hazard for each level of the predictor 
(similarly to what can be obtained using \code{\link{plot.survfit}} and the \code{fun = "cloglog"} option), 
as well as a data set containing that information.}

\references{
Klein, J. and Moeschberger, M. (2003). 
\emph{Survival analysis: techniques for censored and truncated data}. 
2nd edition, Springer.
}

\author{Sarah R. Haile, \email{sarah.haile@gmail.com}} 

\seealso{
Requires \pkg{survival}. A similar plot can be produced using \code{\link{plot.survfit}} and the option \code{fun = "cloglog"}.
}

\examples{
data(larynx)
WeibullDiag(Surv(time, death) ~ stage, data = larynx)
}
\keyword{survival}
\keyword{regression}
