\name{simulate.gs.be}
\alias{simulate.gs.be}
\title{Simulation of Group Sequential Bioequivalence Trial Designs}
\description{Simulates the outcome of a group sequential bioequivalence trial design in order to compute its realised operating characteristics.}
\usage{simulate.gs.be(D = 3, L = 3, sigma.e = sqrt(log(0.3^2 + 1)),
               sigma.b = sqrt(2*sigma.e^2), n = 12,
               a = c(-0.413, 0.623, 1.899), B = log(1.25),
               theta = 0, sequence.type = "williams",
               REML = TRUE, adjust = TRUE, mu.0 = 0,
               pi = rep(0, D - 1), parallel = TRUE, cpus = 8,
               replicates = 10000, summary = TRUE)}
\arguments{
  \item{D}{Number of drugs in the trial initially. Must be an integer greater than or equal to 2. Defaults to 3.}
  \item{L}{Maximum number of stages for the trial. Must be an integer greater than or equal to 1. Defaults to 3.}
  \item{sigma.e}{Value of the within person standard deviation. Must be strictly positive. Defaults to sqrt(log(0.3^2 + 1)).}
  \item{sigma.b}{Value of the between person standard deviation. Must be strictly positive. Defaults to sqrt(2*sigma.e^2).}
  \item{n}{Number of patients to recruit in each stage of the trial. Must be an integer greater than or equal to 1. Defaults to 12.}
  \item{a}{Acceptance boundaries for each interim analysis. Should be a vector of length L.}
  \item{B}{Bioequivalence margin. Must be strictly positive. Defaults to log(1.25).}
  \item{theta}{Range of values for theta, the true fixed effect of each experimental treatment. Defaults to 0.}
  \item{sequence.type}{Sequences to use. Must be set to "latin" or "williams". Defaults to "williams".}
  \item{REML}{Logical. Indicates whether REML estimation should be used at the interim analyses. Defaults to TRUE.}
  \item{adjust}{Logical. Indicates whether the supplied efficacy and futility boundaries should be adjusted using the approach of Whitehead (2009). Defaults to TRUE.}
  \item{mu.0}{The baseline response. Should be a single number. Defaults to 0.}
  \item{pi}{The fixed effect of each period 2 through D, relative to period 1. Should be a vector of length D-1. Defaults to rep(0, D - 1).}
  \item{parallel}{Specifies whether to use parallel computation to speed up the program. Must be logical. Defaults to TRUE.}
  \item{cpus}{Specifies the number of CPUs to parallelise over if parallel is set to TRUE. Must be an integer greater than or equal to 1. Defaults to 8.}
  \item{replicates}{The number of replicate simulations to use for each values of theta to estimate average performance. As this value is increased the programme will increase in running time. Defaults to 10000.}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list with the following components:
  \item{a}{As per input.}
  \item{adjust}{As per input.}
  \item{B}{As per input.}
  \item{cpus}{As per input.}
  \item{D}{As per input.}
  \item{L}{As per input.}
  \item{mu.0}{As per input.}
  \item{n}{As per input.}
  \item{parallel}{As per input.}
  \item{performance}{A matrix with first column theta, second column P(Reject H_01 | theta), third column P(Reject H_0d for some d | theta), fourth column E(N | theta) and fifth column E(O | theta).}
  \item{pi}{As per input.}
  \item{REML}{As per input.}
  \item{replicates}{As per input.}
  \item{sequence.type}{As per input.}
  \item{sigma.b}{As per input.}
  \item{sigma.e}{As per input.}
  \item{theta}{As per input.}
}
\author{Michael J. Grayling}
\seealso{\code{\link{optimal.gs.be}}}
\examples{# Compute an example optimal design:
example <- optimal.gs.be()
# Compute performance of example optimal design under H_G and H_A:
simulated.performance <- simulate.gs.be(example, theta = c(-log(1.25), 0),
                                        adjust = FALSE, REML = FALSE)
}
