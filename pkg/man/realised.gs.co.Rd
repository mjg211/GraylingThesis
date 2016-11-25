\name{realised.gs.co}
\alias{realised.gs.co}
\title{Realised performance of Group Sequential Crossover Trial Designs}
\description{Simulates group sequential crossover trial designs under different analysis procedures, recruitment rates, and period lengths, in order to compute its realised operating characteristics.}
\usage{realised.gs.co(trial.type = "1A", all.data = TRUE, D = 4,
               sigma.e = sqrt(6.51), sigma.b = sqrt(10.12),
               sequence.type = "williams", lambda = 1, w = 10,
               p = 10, n.vec = c(36, 72, 108),
               f = c(0.35, 1.37, 2.08), e = c(2.90, 2.35, 2.08),
               mu.0 = 10.65, tau = rep(0, D - 1),
               pi = c(-0.77, -0.96, -0.55), REML = TRUE,
               parallel = TRUE, cpus = 8, replicates = 10000,
               seed = Sys.time(), summary = TRUE)}
\arguments{
  \item{trial.type}{Trial procedure to use. Must be set to "1A", "1B" or "2". Defaults to "1A".}
  \item{all.data}{Logical. Indicates whether all available or only data marked for a particular stage should be used in analyses. Defaults to TRUE.}
  \item{D}{Number of drugs in the trial initially. Must be an integer greater than or equal to 2. Defaults to 4.}
  \item{sigma.e}{Value of the within person standard deviation. Must be strictly positive. Defaults to sqrt(6.51).}
  \item{sigma.b}{Value of the between person standard deviation. Must be strictly positive. Defaults to sqrt(10.12).}
  \item{sequence.type}{Sequences to use. Must be set to "latin" or "williams". Defaults to "williams".}
  \item{lambda}{Recruitment rate per unit of time. Must be strictly positive. Defaults to 1.}
  \item{w}{Washout period length. Must be an integer greater than or equal to 1. Defaults to 10.}
  \item{p}{Length of each period on a treatment. Must be an integer greater than or equal to 1. Defaults to 10.}
  \item{n.vec}{Vector containing sample sizes after which interim analyses should be conducted. Should contain only increasing integer values. Defaults to c(36, 72, 108).}
  \item{f}{Futility boundaries for each interim analysis. Should be a vector of length length(n.vec).}
  \item{e}{Efficacy boundaries for each interim analysis. Should be a vector of length length(n.vec).}
  \item{mu.0}{The baseline response. Should be a single number. Defaults to 0.}
  \item{tau}{The vector of true treatment effects in each experimental arm. Should be a vector of length D-1. Defaults to rep(0, D - 1).}
  \item{pi}{The fixed effect of each period 2 through D, relative to period 1. Should be a vector of length D-1. Defaults to rep(0, D - 1).}
  \item{REML}{Logical. Indicates whether REML estimation should be used at the interim analyses. Defaults to TRUE.}
  \item{parallel}{Specifies whether to use parallel computation to speed up the program. Must be logical. Defaults to TRUE.}
  \item{cpus}{Specifies the number of CPUs to parallelise over if parallel is set to TRUE. Must be an integer greater than or equal to 1. Defaults to 8.}
  \item{replicates}{The number of replicate simulations to use for each values of theta to estimate average performance. As this value is increased the programme will increase in running time. Defaults to 10000.}
  \item{seed}{A seed is set at the beginning of the search to allow replication of results. Must be compatible with set.seed(). Defaults to Sys.time().}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list with the following components:
  \item{all.data}{As per input.}
  \item{av.results}{Average performance across the replicates.}
  \item{cpus}{As per input.}
  \item{D}{As per input.}
  \item{e}{As per input.}
  \item{f}{As per input.}
  \item{lambda}{As per input.}
  \item{mu.0}{As per input.}
  \item{n.vec}{As per input.}
  \item{p}{As per input.}
  \item{parallel}{As per input.}
  \item{pi}{As per input.}
  \item{REML}{As per input.}
  \item{replicates}{As per input.}
  \item{results}{Matric of the results from each individual trial simulation.}
  \item{seed}{As per input.}
  \item{sequence.type}{As per input.}
  \item{sigma.b}{As per input.}
  \item{sigma.e}{As per input.}
  \item{tau}{As per input.}
  \item{trial.type}{As per input.}
  \item{w}{As per input.}
}
\author{Michael J. Grayling}
\seealso{\code{\link{gs.co}}}
\examples{# Example:
example <- realised.gs.co()
}
