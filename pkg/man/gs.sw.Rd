\name{gs.sw}
\alias{gs.sw}
\title{Group Sequential SW-CRCT Designs}
\description{Computes the required sample size and boundaries of a one-sided group sequential SW-CRCT with error spending boundaries.}
\usage{gs.sw(X, sigma.e = sqrt(0.5), sigma.c = sqrt(0.02), alpha = 0.05,
      beta = 0.1, delta = 0.2, set.T = 2:5, gamma = 0.5,
      summary = TRUE)}
\arguments{
  \item{X}{Matrix of the X_ij, indicating the number of cluster and time periods and when each switches on to the experimental intervention.}
  \item{sigma.e}{Value of the residual standard deviation. Must be strictly positive. Defaults to sqrt(0.5).}
  \item{sigma.c}{Value of the betweeen cluster standard deviation. Must be strictly positive. Defaults to sqrt(0.02).}
  \item{alpha}{Desired type-I error rate. Must be strictly between 0 and 1. Defaults to 0.05.}
  \item{beta}{Desired type-II error rate. Must be strictly between 0 and 1. Defaults to 0.2.}
  \item{delta}{Value of the clinically relevant difference to power at. Must be strictly positive. Defaults to 0.2.}
  \item{set.T}{Vector of time points after which interim analyses will be conducted. Must be a vector of increasing values, where the maximal value is less than or equal to ncol(X), and the minimal values is greater than or equal to min_j {sum(X[, 1:j]) > 0}. Defaults to 2:5.}
  \item{gamma}{Value for the parameter in the family error spending function. Must be greater than or equal to 0. Defaults to 0.5.}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list of class "gsSW", with the following components:
  \item{alpha}{As per input.}
  \item{beta}{As per input.}
  \item{delta}{As per input.}
  \item{f}{Futility stopping boundaries.}
  \item{gamma}{As per input.}
  \item{I}{Information vector for the identified design.}
  \item{Lambda}{Covariance matrix of the standardised treatment effects for the identified design.}
  \item{n}{Number of patients required in each cluster in each time period by the design.}
  \item{perf.H0}{Performance of the identified design when tau = 0.}
  \item{perf.H1}{Performance of the identified design when tau = delta.}
  \item{set.T}{As per input.}
  \item{sigma.c}{As per input.}
  \item{sigma.e}{As per input.}
  \item{X}{As per input.}
}
\author{Michael J. Grayling}
\seealso{\code{\link{plot.gsSW}}}
\examples{# Compute the design for the Bashour et al. trial
X                 <- matrix(0, nrow = 4, ncol = 5)
for (i in 1:4){
  X[i, (i + 1):5] <- 1
}
bashour           <- gs.sw(X = X)
}
