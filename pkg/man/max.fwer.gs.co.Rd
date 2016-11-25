\name{max.fwer.gs.co}
\alias{max.fwer.gs.co}
\title{Maximum FWER of a Group Sequential Crossover Trial Design}
\description{Computes the maximum familywise error rate of a specified group sequential crossover trial design.}
\usage{max.fwer.gs.co(D = 4, L = 3, sigma.e = sqrt(6.51),
               sigma.b = sqrt(10.12), alpha = 0.05,
               beta = 0.2, delta = 1.25, sequences,
               e = c(8.53, 3.18, 1.92), f = c(-0.30, 0.59, 1.92),
               n.or.m = "m", n.m = 3, max.calls = 1000,
               lower = NULL, upper = NULL, init = NULL,
               seed = Sys.time(), summary = TRUE)}
\arguments{
  \item{D}{Number of drugs in the trial initially. Must be an integer greater than or equal to 2. Defaults to 4.}
  \item{L}{Maximum number of stages for the trial. Must be an integer greater than or equal to 1. Defaults to 3.}
  \item{sigma.e}{Value of the within person standard deviation. Must be strictly positive. Defaults to sqrt(6.51).}
  \item{sigma.b}{Value of the between person standard deviation. Must be strictly positive. Defaults to sqrt(6.51).}
  \item{alpha}{Desired familywise error rate. Must be strictly between 0 and 1. Defaults to 0.05.}
  \item{beta}{Desired type-II error rate for H_0^(1). Must be strictly between 0 and 1. Defaults to 0.2.}
  \item{delta}{Value of the clinically relevant difference to power at. Must be strictly positive. Defaults to 1.25.}
  \item{sequences}{Sequences to use in the trial for treatment allocation. Must be set to either "williams", "latin", or be a list where elements 2 through D are matrices containing the sequences to use when 2 through D treatments are remaining in the trial. The values in the matrix of element d must be 0:(d - 1).}
  \item{e}{Efficacy stopping boundaries for the interim analyses. Must be a vector of length L. Defaults to c(8.53, 3.18, 1.92).}
  \item{f}{Futility stopping boundaries for the interim analyses. Must be a vector of length L. Defaults to c(-0.30, 0.59, 1.92).}
  \item{n.or.m}{Whether to interpret the value of n.m as n or m, i.e. whether it is an absolute sample size to use in each stage (n), or the multiple of the specified sequences to use as the sample size (m). Must be set to "n" or "m". Defaults to "m".}
  \item{n.m}{Value of n or m. Must be an integer greater than or equal to 1. Defaults to 3.}
  \item{max.calls}{Maximum number of calls to the objective function to use in the search. Must be an integer greater than or equal to 1. Defaults to 1000.}
  \item{lower}{Lower limits for the treatment effects to use in the search. Should be a vector of length D - 1. Defaults to NULL.}
  \item{upper}{Upper limits for the treatment effects to use in the search. Should be a vector of length D - 1. Defaults to NULL.}
  \item{init}{Starting point for the treatment effects to use in the search. Should be a vector of length D - 1. Defaults to NULL.}
  \item{seed}{A seed is set at the beginning of the search to allow replication of results. Must be compatible with set.seed(). Defaults to Sys.time().}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list with the following components:
  \item{alpha}{As per input.}
  \item{beta}{As per input.}
  \item{D}{As per input.}
  \item{delta}{As per input.}
  \item{e}{As per input.}
  \item{f}{As per input.}
  \item{fwer}{The identified maximal FWER.}
  \item{init}{Initial starting values of the treatment effects used in the search.}
  \item{L}{As per input.}
  \item{lower}{Lower limits for the treatment effects used in the search.}
  \item{max.fwer.search}{Object containing details of the search.}
  \item{n.m}{As per input.}
  \item{n.or.m}{As per input.}
  \item{seed}{As per input.}
  \item{sequences}{Sequences used for patient treatment allocation.}
  \item{sigma.b}{As per input.}
  \item{sigma.e}{As per input.}
  \item{tau}{Vector of treatment effects conveying the maximal FWER.}
  \item{upper}{Upper limits for the treatment effects used in the search.}
}
\author{Michael J. Grayling}
\seealso{\code{\link{gs.co}}}
\examples{# Example
example <- max.fwer.gs.co()
}
