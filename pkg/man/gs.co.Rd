\name{gs.co}
\alias{gs.co}
\title{Group Sequential Crossover Trial Designs}
\description{Computes the required sample size and boundaries of a one-sided group sequential crossover trial with power family boundaries.}
\usage{gs.co(D = 4, L = 3, sigma.e = sqrt(6.51), delta = 1.25,
      alpha = 0.05, beta = 0.2, sequence.type = "williams",
      Delta = 0, summary = TRUE)}
\arguments{
  \item{D}{Number of drugs in the trial initially. Must be an integer greater than or equal to 2. Defaults to 4.}
  \item{L}{Maximum number of stages for the trial. Must be an integer greater than or equal to 1. Defaults to 3.}
  \item{sigma.e}{Value of the within person standard deviation. Must be strictly positive. Defaults to sqrt(6.51).}
  \item{delta}{Value of the clinically relevant difference to power at. Must be strictly positive. Defaults to 1.25.}
  \item{alpha}{Desired familywise error rate. Must be strictly between 0 and 1. Defaults to 0.05.}
  \item{beta}{Desired type-II error rate for H_0^(1). Must be strictly between 0 and 1. Defaults to 0.2.}
  \item{sequence.type}{Sequence type to plan for. Must be set to "latin" or "williams". Defaults to "williams".}
  \item{Delta}{Value for the shape parameter of the power family boundaries. Must be less than or equal to 0.5. Defaults to 0.}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list of class "gsCO", with the following components:
  \item{alpha}{As per input.}
  \item{beta}{As per input.}
  \item{D}{As per input.}
  \item{Delta}{As per input.}
  \item{delta}{As per input.}
  \item{design.performance}{Summary of the performance of the design.}
  \item{exact.design.performance}{Summary of the performance of the exact design.}
  \item{e}{Efficacy stopping boundaries.}
  \item{f}{Futility stopping boundaries.}
  \item{I}{Information vector for the rounded identified design.}
  \item{I.1}{Information vector of a single treatment effect for the rounded identified design.}
  \item{I.exact}{Information vector for the exact design.}
  \item{I.1.exact}{Information vector of a single treatment effect for the exact design.}
  \item{L}{As per input.}
  \item{Lambda}{Covariance matrix of the standardised treatment effects for the identified design.}
  \item{n}{Number of patients required in each stage by the rounded design.}
  \item{n.exact}{Number of patients required in each stage by the exact design.}
  \item{scenarios}{The different possible outcomes for the trial and their respective probabilities given the identified design.}
  \item{scenarios.exact}{The different possible outcomes for the trial and their respective probabilities given the identified exact design.}
  \item{sequence.type}{As per input.}
  \item{sigma.e}{As per input.}
}
\author{Michael J. Grayling}
\seealso{\code{\link{plot.gsCO}}}
\examples{# Compute the design for the TOMADO trial parameters with L = 3 stages and Delta = 0:
tomado.3_0 <- gs.co()
# Compute the design for the TOMADO trial parameters with L = 2 stages and Delta = 0.5:
tomado.2_0.5 <- gs.co(L = 2, Delta = 0.5)
}
