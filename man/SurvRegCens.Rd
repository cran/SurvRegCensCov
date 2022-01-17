\name{SurvRegCens}
\alias{SurvRegCens}
\title{Weibull Survival Regression Model with a censored covariate}
\description{Computes estimators for the shape and scale parameter of the Weibull distribution, as well as for the vector of regression parameters in a parametric survival model with potentially right-censored time-to-event endpoint distributed according to a Weibull distribution. The regression allows for one potentially interval-censored and an arbitrary number of non-censored covariates.}
\usage{
SurvRegCens(formula, data = parent.frame(), Density, initial, conf.level = 0.95, 
          intlimit = 10^-10, namCens = "VarCens", trace = 0, reltol = 10^-8)
}
\arguments{
  \item{formula}{A formula expression as for other regression models. The response has to be a survival object for right-censored data, as returned by the \code{\link{Surv}} function. The censored covariate is equally specified using \code{\link{Surv}}. See the examples below and the documentation for \code{\link{Surv}}, \code{\link{lm}} and \code{\link{formula}} for details.}
  \item{data}{A data frame in which to interpret the variables named in the formula argument.}
  \item{Density}{Density function of the censored covariate.}
  \item{initial}{Initial values for the parameters to be optimized over, ordered according to Scale parameter, Shape parameter, regression parameters (i.e. \eqn{\beta}) linked to the non-censored covariates, regression parameter (i.e. \eqn{\beta}) linked to the censored covariate. A straightforward initial vector is based on ignoring the censoring of the censored covariate and using \code{\link{survreg}}, see the example below for an illustration.}
  \item{conf.level}{Confidence level of confidence intervals.}  
  \item{intlimit}{In computation of integrals, values of the function to be integrated below \code{intlimit} are set to 0. This makes integration results more accurate and speeds up integration. If the data is such that the absolute values of the underlying baseline Weibull density are very small, i.e. in the range of \code{intlimit}, it is advisable to rescale the time variable, e.g. change the scaling from days to years. A very small value of the estimated \eqn{\lambda} is indicative of that situation.}
  \item{namCens}{Name of censored covariate, to tidy outputs.}
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
\code{\link{SurvRegCens}} returns an object of class \code{"src"}, a list containing the 
following components:
	\item{coeff}{Estimators, confidence intervals, \eqn{p}-values for the for the null hypothesis: \{Estimators is equal to 0\}, and this for each of the parameters of the Weibull survival regression model.}
	\item{percent.cens}{Percentage of censored observations in the censored covariate.}
	\item{loglik}{Log-likelihood function value at the estimators.}
	\item{info.converg}{Convergence information provided by the function \code{\link{optim}}.}
	\item{info.converg.message}{Message provided by \code{\link{optim}}.}	
The methods \code{\link{print.src}}, \code{\link{summary.src}}, \code{\link{coef.src}}, and \code{\link{logLik.src}} are used to print or obtain a summary, coefficients, or the value of the log-likelihood at the maximum.
}

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
\dontrun{
## --------------------------------------------------------------
## 1 censored-covariate and 2 non-censored covariates 
## no censoring, to compare result with survival::survreg
## modify prop.cens to introduce left-censoring of covariate
## --------------------------------------------------------------

set.seed(158)
n <- 100
lambda <- exp(-2)
gamma <- 1.5

## vector of regression parameters: the last entry is the one for the censored covariate
beta <- c(0.3, -0.2, 0.25) 
true <- c(lambda, gamma, beta)

## non-censored covariates
var1 <- rnorm(n, mean = 4, sd = 0.5)
var2 <- rnorm(n, mean = 4, sd = 0.5)

## Generate censored covariate. 
## For generation of Weibull survival times, do not left-censor it yet.
var3 <- rnorm(n, mean = 5, sd = 0.5)

## simulate from a Weibull regression model
time <- TimeSampleWeibull(covariate_noncens = data.frame(var1, var2), 
          covariate_cens = var3, lambda = lambda, gamma = gamma, beta = beta) 

## left-censor covariate
## prop.cens specifies the proportion of observations that should be left-censored
prop.cens <- 0
LOD <- qnorm(prop.cens, mean = 5, sd = 0.5)
var3.cens <- censorContVar(var3, LLOD = LOD)

## censor survival time
event <- matrix(1, nrow = n, ncol = 1)
time.cens <- rexp(n, rate = 0.5)
ind.time <- (event >= time.cens)
event[ind.time] <- 0
time[ind.time] <- time.cens[ind.time]

## specify the density for the censored covariate:
## For simplicity, we take here the "true" density we simulate from. In an application,
## you might want to use a density with parameters estimated from the censored covariate,
## e.g. using the function ParamSampleCens. See example in Hubeaux & Rufibach (2014).
DensityCens <- function(value){return(dnorm(value, mean = 5, sd = 0.5))}

## use Weibull regression where each censored covariate value is set 
## to LOD ("naive" method)
naive <- survreg(Surv(time, event) ~ var1 + var2 + var3.cens[, 2], dist = "weibull")
initial <- as.vector(ConvertWeibull(naive)$vars[, 1])

## use new method that takes into account the left-censoring of one covariate
data <- data.frame(time, event, var3.cens, var1, var2)
formula <- formula(Surv(time, event) ~  Surv(time = var3.cens[, 1], time2 = var3.cens[, 2], 
                      type = "interval2") + var1 + var2)
cens1 <- SurvRegCens(formula = formula, data = data, Density = DensityCens, initial = initial, 
                      namCens = "biomarker")
summary(cens1)
coef(cens1)
logLik(cens1)

## compare estimates
tab <- data.frame(cbind(true, initial, cens1$coeff[, 1]))
colnames(tab) <- c("true", "naive", "Weibull MLE")
rownames(tab) <- rownames(cens1$coeff)
tab

## compare confidence intervals
ConvertWeibull(naive)$HR[, 2:3]
cens1$coeff[, 7:8]


## --------------------------------------------------------------
## model without the non-censored covariates
## --------------------------------------------------------------
naive2 <- survreg(Surv(time, event) ~ var3.cens[, 2], dist = "weibull")
initial2 <- as.vector(ConvertWeibull(naive2)$vars[, 1])

## use new method that takes into account the left-censoring of one covariate
formula <- formula(Surv(time, event) ~ Surv(time = var3.cens[, 1], time2 = var3.cens[, 2], 
                      type = "interval2"))
cens2 <- SurvRegCens(formula = formula, data = data, Density = DensityCens, initial = initial2, 
                      namCens = "biomarker")
summary(cens2)

## compare estimates
tab <- data.frame(cbind(true[c(1, 2, 5)], initial2, cens2$coeff[, 1]))
colnames(tab) <- c("true", "naive", "Weibull MLE")
rownames(tab) <- rownames(cens2$coeff)
tab

## compare confidence intervals
ConvertWeibull(naive2)$HR[, 2:3]
cens2$coeff[, 7:8]
}
}
\keyword{methods}
\keyword{survival}
\keyword{regression}
