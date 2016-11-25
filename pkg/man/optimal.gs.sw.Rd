\name{optimal.gs.sw}
\alias{optimal.gs.sw}
\title{Optimal Group Sequential SW-CRCT Designs}
\description{Computes the required sample size and boundaries of optimal one-sided group sequential SW-CRCT designs.}
\usage{optimal.gs.sw(C = 4, Ti = 5, sigma.c = sqrt(0.02),
              sigma.e = sqrt(0.5), alpha = 0.05, beta = 0.1,
              delta = 0.2, set.T = 2:5, eta = -1, penalty = 1,
              seed = Sys.time(), summary = TRUE)}
\arguments{
  \item{C}{Number of clusters in the trial. Must be an integer greater than or equal to 2. Defaults to 4.}
  \item{Ti}{Number of time periods in the trial. Must be an integer greater than or equal to 2. Defaults to 5.}
  \item{sigma.c}{Value of the betweeen cluster standard deviation. Must be strictly positive. Defaults to sqrt(0.02).}
  \item{sigma.e}{Value of the residual standard deviation. Must be strictly positive. Defaults to sqrt(0.5).}
  \item{alpha}{Desired type-I error rate. Must be strictly between 0 and 1. Defaults to 0.05.}
  \item{beta}{Desired type-II error rate. Must be strictly between 0 and 1. Defaults to 0.2.}
  \item{delta}{Value of the clinically relevant difference to power at. Must be strictly positive. Defaults to 0.2.}
  \item{set.T}{Vector of time points after which interim analyses will be conducted. Must be a vector of increasing values, where the maximal value is less than or equal to ncol(X), and the minimal values is greater than or equal to min_j {sum(X[, 1:j]) > 0}. Defaults to 2:5.}
  \item{eta}{Sets the treatment effect we minimise the expected sample size for as tau = eta*delta. Defaults to -1.}
  \item{penalty}{Multiple of the fixed design sample size to use as the penalty constant in the objective function. Must be strictly positive. Defaults to 1.}
  \item{seed}{A seed is set at the beginning of the search to allow replication of results. Must be compatible with set.seed(). Defaults to Sys.time().}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list of class "gsSW", with the following components:
  \item{alpha}{As per input.}
  \item{beta}{As per input.}
  \item{C}{As per input.}
  \item{delta}{As per input.}
  \item{eta}{As per input.}
  \item{f}{Identified optimal futility stopping boundaries.}
  \item{I}{Information vector for the identified design.}
  \item{Lambda}{Covariance matrix of the standardised treatment effects for the identified optimal design.}
  \item{n}{Number of patients required in each cluster in each time period by the optimal design.}
  \item{penalty}{As per input.}
  \item{perf.eta.delta}{Performance of the identified design when tau = eta*delta.}
  \item{perf.H0}{Performance of the identified design when tau = 0.}
  \item{perf.H1}{Performance of the identified design when tau = delta.}
  \item{run.time}{Run time of the programme.}
  \item{seed}{As per input.}
  \item{set.T}{As per input.}
  \item{sigma.c}{As per input.}
  \item{sigma.e}{As per input.}
  \item{Ti}{As per input.}
  \item{X}{As per input.}
}
\author{Michael J. Grayling}
\seealso{\code{\link{plot.gsSW}}}
\examples{# Compute the design for the Bashour et al. trial
  bashour           <- optimal.gs.sw()
}
