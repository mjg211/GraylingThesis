\name{ssre.sw}
\alias{ssre.sw}
\title{Sample Size Re-estimation Procedures for SW-CRCT Designs}
\description{Computes the performance of two types of sample size re-estimation procedure for SW-CRCT designs.}
\usage{ssre.sw(X, sigma.e.tilde = sqrt(0.5), sigma.e = sqrt(0.5),
        sigma.c.tilde = sqrt(0.02), sigma.c = sqrt(0.02),
        alpha = 0.05, beta = 0.1, delta = 0.2, tau = 0,
        pi = rep(0, ncol(X)), lambda.L = 1, n.max = 200, t = 1,
        blinded = TRUE, correction = 0, REML = TRUE,
        parallel = TRUE, cpus = 8, replicates = 10000,
        summary = TRUE)}
\arguments{
 \item{X}{Matrix of the X_ij, indicating the number of cluster and time periods and when each switches on to the experimental intervention.}
 \item{sigma.e.tilde}{Assumed value of the residual standard deviation. Must be strictly positive. Defaults to sqrt(0.5).}
  \item{sigma.e}{True value of the residual standard deviation. Must be strictly positive. Defaults to sqrt(0.5).}
  \item{sigma.c.tilde}{Assumed value of the betweeen cluster standard deviation. Must be strictly positive. Defaults to sqrt(0.02).}
  \item{sigma.c}{True value of the betweeen cluster standard deviation. Must be strictly positive. Defaults to sqrt(0.02).}
  \item{alpha}{Desired type-I error rate. Must be strictly between 0 and 1. Defaults to 0.05.}
  \item{beta}{Desired type-II error rate. Must be strictly between 0 and 1. Defaults to 0.2.}
  \item{delta}{Value of the clinically relevant difference to power at. Must be strictly positive. Defaults to 0.2.}
  \item{tau}{True value of the treatment effect. Defaults to 0.}
  \item{pi}{The fixed effect of each period 1 through D. Should be a vector of length ncol(X). Defaults to rep(0, ncol(X)).}

  \item{lambda.L}{Minimum proportion of the initial assumed required per cluster per period sample size allowed as the re-estimated value. Must be strictly positive. Defaults to 1.}
  \item{n.max}{Maximum per cluster per period sample size allowed after re-estimation.}
  \item{t}{Time period after which to perform the SSRE.}
  \item{blinded}{Logical. Indicates whether to use the blinded or unblinded SSRE procedure. Defaults to TRUE.}
  \item{correction}{Assumed value of tau to use in the formula for estimating the between cluster variance in the blinded re-estimation procedure. Defaults to 0.}
  \item{REML}{Logical. Indicates whether REML estimation should be used at the interim analyses. Defaults to TRUE.}
  \item{parallel}{Specifies whether to use parallel computation to speed up the program. Must be logical. Defaults to TRUE.}
  \item{cpus}{Specifies the number of CPUs to parallelise over if parallel is set to TRUE. Must be an integer greater than or equal to 1. Defaults to 8.}
  \item{replicates}{The number of replicate simulations to use to estimate average performance. As this value is increased the programme will increase in running time. Defaults to 10000.}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list of class "gsSW", with the following components:
  \item{alpha}{As per input.}
  \item{av.results}{Average results across the simulated trials.}
  \item{beta}{As per input.}
  \item{blinded}{As per input.}
  \item{correction}{As per input.}
  \item{cpus}{As per input.}
  \item{delta}{As per input.}
  \item{lambda.L}{As per input.}
  \item{n.max}{As per input.}
  \item{parallel}{As per input.}
  \item{pi}{As per input.}
  \item{REML}{As per input.}
  \item{replicates}{As per input.}
  \item{results}{Results of the individual simulated trials.}
  \item{sigma.c}{As per input.}
  \item{sigma.c.tilde}{As per input.}
  \item{sigma.e}{As per input.}
  \item{sigma.e.tilde}{As per input.}
  \item{t}{As per input.}
  \item{tau}{As per input.}
  \item{X}{As per input.}
}
\author{Michael J. Grayling}
\examples{# Compute the performance of the blinded procedure
  bashour           <- ssre.sw()
}
