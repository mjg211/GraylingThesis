\name{ssre.be}
\alias{ssre.be}
\title{Sample Size Re-estimation Procedures for Multi-formulation Bioequivalence Trials}
\description{Determines the performance of several sample size re-estimation procedures for multi-formulation bioequivalence trials.}
\usage{ssre.be(method = "A", p = 0.5, D = 4,
        sigma.e.tilde = sqrt(log(0.3^2 + 1)),
        sigma.e = sqrt(log(0.3^2 + 1)),
        sigma.b = sqrt(2*sigma.e^2), alpha = 0.05, beta = 0.2,
        B = log(1.25), delta = 0, theta = 0, pi = rep(0, D - 1),
        n.max = 48, current = FALSE, REML = TRUE,
        t.boundaries = TRUE, parallel = TRUE, cpus = 8,
        replicates = 10000, summary = TRUE)}
\arguments{
  \item{method}{SSRE procedure to use. Must be set to "A", "B" or "C". Defaults to "A".}
  \item{p}{Proportion of initial estimated required sample size to perform the interim analysis at. Must be strictly positive. Defaults to 0.5.}
  \item{D}{Number of drugs in the trial initially. Must be an integer greater than or equal to 2. Defaults to 4.}
  \item{sigma.e.tilde}{Assumed value of the within person standard deviation. Must be strictly positive. Defaults to sqrt(log(0.3^2 + 1)).}
  \item{sigma.e}{True value of the within person standard deviation. Must be strictly positive. Defaults to sqrt(log(0.3^2 + 1)).}
  \item{sigma.b}{True value of the between person standard deviation. Must be strictly positive. Defaults to sqrt(2*sigma.e^2).}
  \item{alpha}{Desired familywise error rate. Must be strictly between 0 and 1. Defaults to 0.05.}
  \item{beta}{Desired type-II error rate for H_0^(1). Must be strictly between 0 and 1. Defaults to 0.2.}
  \item{B}{Bioequivalence margin. Must be strictly positive. Defaults to log(1.25).}
  \item{delta}{Bioequivalence margin to power at. Must be greater than or equal to 0 and strictly less than B. Defaults to 0.}
  \item{theta}{True fixed effect of each experimental treatment. Defaults to 0.}
  \item{pi}{The fixed effect of each period 2 through D, relative to period 1. Should be a vector of length D-1. Defaults to rep(0, D - 1).}
  \item{n.max}{Maximum allowed sample size for the trial. Must be an integer greater than or equal to 2. Defaults to 48.}
  \item{current}{Logical. Indicates whether rejection boundaries should be computed based on the current number of treatments remaining in the trial at that analysis, or the initial number in the trial. Defaults to FALSE.}
  \item{REML}{Logical. Indicates whether REML estimation should be used at the interim analyses. Defaults to TRUE.}
  \item{t.boundaries}{Logical. Indicates whether t or normal distribution boundaries should be used at the interim analyses. Defaults to TRUE.}
  \item{parallel}{Specifies whether to use parallel computation to speed up the program. Must be logical. Defaults to TRUE.}
  \item{cpus}{Specifies the number of CPUs to parallelise over if parallel is set to TRUE. Must be an integer greater than or equal to 1. Defaults to 8.}
  \item{replicates}{The number of replicate simulations to use to estimate average performance. As this value is increased the programme will increase in running time. Defaults to 10000.}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list with the following components:
  \item{all.results}{The result of each simulated trial.}
  \item{alpha}{As per input.}
  \item{alpha.BC}{Determined stage-wise significance level to use in method "B" or "C".}
  \item{av.results}{The average result across the simulated trials.}
  \item{B}{As per input.}
  \item{beta}{As per input.}
  \item{cpus}{As per input.}
  \item{current}{As per input.}
  \item{D}{As per input.}
  \item{delta}{As per input.}
  \item{method}{As per input.}
  \item{n.max}{As per input.}
  \item{p}{As per input.}
  \item{parallel}{As per input.}
  \item{pi}{As per input.}
  \item{REML}{As per input.}
  \item{replicates}{As per input.}
  \item{sigma.b}{As per input.}
  \item{sigma.e}{As per input.}
  \item{sigma.e.tilde}{As per input.}
  \item{t.boundaries}{As per input.}
  \item{theta}{As per input.}
}
\author{Michael J. Grayling}
\seealso{\code{\link{optimal.gs.be}}}
\examples{# Determine performance of method "B" for defaults parameters:
method.B <- ssre.be(method = "B")
}
