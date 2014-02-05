\name{SurvRegCens}
\alias{SurvRegCens}
\title{Weibull Survival Regression Model with a censored covariate}
\description{Computes estimators for the shape and scale parameter of the Weibull distribution, as well as for the vector of regression parameters in a parametric survival model with potentially right-censored time-to-event endpoint distributed according to a Weibull distribution. The regression allows for one potentially interval-censored and an arbitrary number of non-censored covariates.}
\usage{
SurvRegCens(time, event, CovariateNonCens = NULL, CovariateCens, 
            Density, initial, conf.level = 0.95, intlimit = 10^-10, 
            namCens = "VarCens", trace = 0, reltol = 10^-8)
}
\arguments{
  \item{time}{Vector of event times.}
  \item{event}{Censoring indicator for the time-to-event endpoint (0: censored, 1: event).}
  \item{CovariateNonCens}{\code{data.frame} containing the non-censored covariates. If no non-censored covariates are to be included in the model, set to \code{NULL}.}
  \item{CovariateCens}{Information on censored covariate. Format as for \code{censdata} in \code{\link{ParamSampleCens}}.}
  \item{Density}{Density function of the censored covariate.}
  \item{initial}{Initial values for the parameters to be optimized over, ordered according to Scale parameter, Shape parameter, regression parameters (i.e. \eqn{\beta}) linked to the non-censored covariates, regression parameter (i.e. \eqn{\beta}) linked to the censored covariate. A straightforward initial vector is based on ignoring the censoring of the censored covariate and using \code{\link{survreg}}, see the example below for an illustration.}
  \item{conf.level}{Confidence level of confidence intervals.}  
  \item{intlimit}{In computation of integrals, values of the function to be integrated below \code{intlimit} are set to 0. This makes integration results more accurate and speeds up integration. If the data is such that the absolute values of the underlying baseline Weibull density are very small, i.e. in the range of \code{intlimit}, it is advisable to rescale the time variable, e.g. change the scaling from days to years. A very small value of the estimated \eqn{\lambda} is indicative of that situation.}
  \item{namCens}{Name of censored variables, to tidy outputs.}
  \item{trace}{\code{trace} argument in \code{optim}, indicates whether to show optimization progress.}
  \item{reltol}{\code{reltol} argument in \code{optim}. By changing this one can modify the relative tolerance in maximization of the likelihood function.}
}

\details{
The time-to-event distributed according to a Weibull distribution, i.e. time-to-event \eqn{\sim}{~} Weibull\eqn{(\lambda,\gamma)}, has conditional density given by,
\deqn{f_{Y_i}(t|\mathbf{x}_i,\boldsymbol{\beta}) =\gamma \lambda t^{\gamma-1} \exp\left(\mathbf{x}_i\boldsymbol{\beta}\right)\exp\left( - \lambda t^{\gamma} \exp\left(\mathbf{x}_i\boldsymbol{\beta}\right) \right),}{f_{Y_i}(t|x_i,\beta) =\gamma * \lambda * t^(\gamma-1) * exp(x_i*\beta) * exp(- \lambda * t^\gamma * exp(x_i*\beta)),}
conditional hazard function given by,
\deqn{h_i(t|\mathbf{x}_i,\boldsymbol{\beta})= \lambda \gamma t^{\gamma-1} \exp\left( \mathbf{x}_i\boldsymbol{\beta} \right),}{h_i(t|x_i,\beta)= \lambda * \gamma * t^(\gamma-1) * exp(x_i*\beta),}
and conditional survival function given by,
\deqn{S_i(t|\mathbf{x}_i,\boldsymbol{\beta}) = \exp\left(-\lambda t^{\gamma} \exp\left(\mathbf{x}_i\boldsymbol{\beta}\right)\right),}{S_i(t|x_i,\beta) = exp(- \lambda * t^\gamma * exp(x_i*\beta),}
where \eqn{\mathbf{x}_i}{x_i} collects the values of each covariate for observation \eqn{i} and \eqn{\boldsymbol{\beta}}{\beta} represents the regression parameters.
}

\value{
	\item{coeff}{Estimators, confidence intervals, \eqn{p}-values for the for the null hypothesis: \{Estimators is equal to 0\}, and this for each of the parameters of the Weibull survival regression model.}
	\item{percent.cens}{Percentage of censored observations in the censored covariate.}
	\item{loglik}{Log-likelihood function value at the estimators.}
	\item{info.converg}{Convergence information provided by the function \code{\link{optim}}.}
	\item{info.converg.message}{Message provided by \code{\link{optim}}.}	
}

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
\dontrun{
library(survival)

## --------------------------------------------------------------
## 1 censored-covariate and 2 non-censored covariates 
## no censoring, to compare result with survival::survreg
## modify prop.cens to introduce censoring of covariate
## --------------------------------------------------------------

set.seed(158)
size <- 100
lambda <- exp(-2)
gamma <- 1.5

## vector of regression parameters: the last entry is the one for the censored covariate
ParamBeta <- c(0.3, -0.2, 0.25) 
true <- c(lambda, gamma, ParamBeta)

## simulate from a Weibull regression model
Censdata <- rnorm(size, mean = 5, sd = 0.5)
min.value <- min(Censdata)
prop.cens <- 0
LOD <- qnorm(prop.cens, mean = 5, sd = 0.5)
Censdata[Censdata <= LOD] <- LOD
index.noncens <- which(Censdata > LOD)
InfoCens <- matrix(nrow = size, ncol = 1)
InfoCens[index.noncens] <- Censdata[index.noncens]
CovariateCens <- cbind(InfoCens, Censdata)

NonCensdata <- data.frame("var1" = rnorm(size, mean = 4, sd = 0.5), 
    "var2" = rnorm(size, mean = 4, sd = 0.5))
index.cens <- which(Censdata <= LOD)
index.cens <- index.cens[index.cens >= 250]
InfoCens[index.cens] <- min.value
time <- TimeSampleWeibull(covariate_noncens = NonCensdata, covariate_cens = Censdata, 
                          lambda = lambda, gamma = gamma, beta=ParamBeta) 

IndicTime <- matrix(1, nrow = size, ncol = 1)
TimeCensVector <- rweibull(size, 1.5, exp(3))
IndicTime[time >= TimeCensVector] <- 0
index.cens.time <- which(IndicTime >= TimeCensVector)
time[index.cens.time] <- TimeCensVector[index.cens.time]

## specify the true density for the censored covariate:
DensityCens <- function(value){return(dnorm(value, mean = 5 , sd = 0.5))}

## use Weibull regression where each censored covariate value is set 
## to LOD ("naive" method)
naive <- survreg(Surv(time, IndicTime) ~ NonCensdata[,1] + 
                    NonCensdata[,2] + Censdata, dist = "weibull")
initial <- as.vector(ConvertWeibull(naive)$vars[, 1])

## use new method that takes into account the left-censoring of one covariate
cens1 <- SurvRegCens(time = time, event = IndicTime, CovariateNonCens = NonCensdata, 
                     CovariateCens = CovariateCens, Density = DensityCens,
                     initial = initial, namCens = "biomarker")

## compare estimates
tab <- data.frame(cbind(true, initial, cens1$coeff[, 1]))
colnames(tab) <- c("true", "naive", "Weibull MLE")
tab

## compare confidence intervals
ConvertWeibull(naive)$HR[, 2:3]
cens1$coeff[, 7:8]



## --------------------------------------------------------------
## model without the non-censored covariates
## --------------------------------------------------------------
naive2 <- survreg(Surv(time, IndicTime) ~ Censdata, dist = "weibull")
initial2 <- as.vector(ConvertWeibull(naive2)$vars[, 1])

## use new method that takes into account the left-censoring of one covariate
cens2 <- SurvRegCens(time = time, event = IndicTime, CovariateNonCens = NULL, 
                     CovariateCens = CovariateCens, Density = DensityCens,
                     initial = initial2, namCens = "biomarker")

## compare estimates
tab <- data.frame(cbind(true[c(1, 2, 5)], initial2, cens2$coeff[, 1]))
colnames(tab) <- c("true", "naive", "Weibull MLE")
tab

## compare confidence intervals
ConvertWeibull(naive2)$HR[, 2:3]
cens2$coeff[, 7:8]
}
}
\keyword{methods}
\keyword{survival}
\keyword{regression}
