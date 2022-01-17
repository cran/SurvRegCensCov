\name{larynx}
\Rdversion{1.1}
\alias{larynx}
\docType{data}
\title{Survival Times of Larynx Cancer Patients}
\description{
A study of 90 males with laryngeal cancer was performed, comparing survival times. Each patient's age, year of 
diagnosis, and disease stage was noted, see Kardaun (1983) and Klein and Moeschberger (2003).
}
\usage{data(larynx)}
\format{
  A data frame with 90 observations on the following 5 variables.
  \describe{
    \item{\code{stage}}{Disease stage (1-4) from TNM cancer staging classification.}
    \item{\code{time}}{Time from first treatment until death, or end of study.}
    \item{\code{age}}{Age at diagnosis.}
    \item{\code{year}}{Year of diagnosis.}
    \item{\code{death}}{Indicator of death [1, if patient died at time t; 0, otherwise].}
  }
}
\source{
\url{https://www.mcw.edu/-/media/MCW/Departments/Biostatistics/datafromsection18.txt?la=en}
}
\references{
Kardaun, O. (1983).
Statistical survival analysis of male larynx-cancer patients-a case study. 
\emph{Statistica Neerlandica}, \bold{37}, 103--125.

Klein, J. and Moeschberger, M. (2003). 
\emph{Survival analysis: techniques for censored and truncated data}. 
2nd edition, Springer.
}
\examples{
library(survival)
data(larynx)
Surv(larynx$time, larynx$death)
}
\keyword{datasets}
