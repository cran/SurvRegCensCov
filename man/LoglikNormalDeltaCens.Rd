\name{LoglikNormalDeltaCens}
\alias{LoglikNormalDeltaCens}
\title{Log likelihood function to compute mean difference between two normally distributed censored samples.}
\description{Reparametrization of the log likelihood function for a normally distributed censored sample such that the mean difference is a parameter of the function, thus allowing to be made inference on. The mean difference is computed as sample 1 - sample 2.}
\usage{
LoglikNormalDeltaCens(x, data1, lowerbound1, vdelta1, data2,
                      lowerbound2, vdelta2)
}
\arguments{
  \item{x}{A vector of four components where the first component corresponds to the mean of the normal distribution of the first sample, the second component corresponds to mean difference between the two samples: sample 1 - sample 2, the third component corresponds to the standard deviation of the normal distribution of the first sample, and the fourth component corresponds to the standard deviation of the normal distribution of the second sample.}
  \item{data1}{A vector of data corresponding to the first sample.}
  \item{lowerbound1}{A vector which corresponds to the lower bounds for the interval-censored observations of the vector of data corresponding to the first sample. If no lower bound is available then put \code{NA}.}
  \item{vdelta1}{A vector which indicates for censoring for the first sample (0: censored, 1: not censored).}
  \item{data2}{A vector of data corresponding to the second sample.}
  \item{lowerbound2}{A vector which corresponds to the lower bounds for the interval-censored observations of the vector of data corresponding to the second sample. If no lower bound is available then put \code{NA}.}
  \item{vdelta2}{A vector which indicates for censoring for the second sample (0: censored, 1: not censored).}
}

\references{
Hubeaux, S. (2013). 
Estimation from left- and/or interval-censored samples. 
Technical report, Biostatistics Oncology, F. Hoffmann-La Roche Ltd.  

Lynn, H. S. (2001). 
Maximum likelihood inference for left-censored HIV RNA data. 
\emph{Stat. Med.}, \bold{20}, 33--45.
}

\author{
Stanislas Hubeaux, \email{stanislas.hubeaux@alumni.epfl.ch} 

Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}
}

\note{Function not intended to be invoked by the user.}

\keyword{methods}
