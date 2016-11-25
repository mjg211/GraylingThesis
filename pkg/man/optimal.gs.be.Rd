\name{optimal.gs.be}
\alias{optimal.gs.be}
\title{Determination of Optimal Group Sequential Bioequivalence Trial Designs}
\description{Determines the optimal group size and boundaries for a group sequential bioequivalence trial design.}
\usage{optimal.gs.be(algorithm = "DEoptim", D = 3, L = 2,
              sigma.e = sqrt(log(0.3^2 + 1)), alpha = 0.05,
              beta = 0.2, B = log(1.25), delta = 0,
              sequence.type = "williams", w = rep(1/3, 3),
              opt.crit = "O", penalty = 1, num.runs = 1,
              find.constrained = TRUE, max.calls = 10000,
              lower = NULL, upper = NULL, init = NULL,
              init.file.ext = 1, seed = Sys.time(),
              summary = TRUE)}
\arguments{
  \item{D}{Algorithm to use for determining the optimal design. Must be set to "DEopt", "DEoptim", "ga", "GenSA", "hydroPSO", "NLOPT_GN_CRS2_LM", "NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L", "NLOPT_GN_ISRES", "PSopt", "psoptim", "SANN" or "soma". Defaults to "DEoptim".}
  \item{D}{Number of drugs in the trial initially. Must be an integer greater than or equal to 2. Defaults to 3.}
  \item{L}{Maximum number of stages for the trial. Must be an integer greater than or equal to 1. Defaults to 2.}
  \item{sigma.e}{Value of the within person standard deviation. Must be strictly positive. Defaults to sqrt(log(0.3^2 + 1)).}
  \item{alpha}{Desired familywise error rate. Must be strictly between 0 and 1. Defaults to 0.05.}
  \item{beta}{Desired type-II error rate for H_0^(1). Must be strictly between 0 and 1. Defaults to 0.2.}
  \item{B}{Bioequivalence margin. Must be strictly positive. Defaults to log(1.25).}
  \item{delta}{Bioequivalence margin to power at. Must be greater than or equal to 0 and strictly less than B. Defaults to 0.}
  \item{sequence.type}{Sequences to use. Must be set to "latin" or "williams". Defaults to "williams".}
  \item{w}{Vector or weights to use in objective function. Defaults to rep(1/3, 3).}
  \item{opt.crit}{Optimality criteria to use. Must be set to "N" or "O". Defaults to "O".}
  \item{penalty}{Multiple of the single stage sample size to use as the penalty constant in the objective function. Must be strictly positive. Defaults to 1.}
  \item{num.runs}{Number of times the optimal design should be searched for. Must be an integer greater than or equal to 1. Defaults to 1.}
  \item{find.constrained}{Logical. Determines if the constrained optimal design should be searched for, or just the exact optimal design. Defaults to TRUE.}
  \item{max.calls}{Maximum number of calls to the objective function to use in the search. Must be an integer greater than or equal to 1. Defaults to 10000.}
  \item{lower}{Lower limits for the parameters to use in the search. Should be a vector of length L + 1. Defaults to NULL.}
  \item{upper}{Upper limits for the parameters to use in the search. Should be a vector of length L + 1. Defaults to NULL.}
  \item{init}{Starting point to use in the search. Should be a vector of length L + 1. Defaults to NULL.}
  \item{init.file.ext}{.csv files containing the optimal designs will be written to the current working directory with the name algorithm_BLANK.csv. BLANK takes the values init.file.ext, ... , init.file.ext + num.runs. Must be an integer. Defaults to 1.}
  \item{seed}{A seed is set at the beginning of the search to allow replication of results. Must be compatible with set.seed(). Defaults to Sys.time().}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
}
\value{Returns a list of class "gsBE", with the following components:
  \item{a}{Acceptance stopping boundaries of the constained optimal design.}
  \item{a.exact}{Acceptance stopping boundaries of the exact optimal design.}
  \item{all.run.information}{Matrix containing the results of the performed searches.}
  \item{alpha}{As per input.}
  \item{B}{As per input.}
  \item{beta}{As per input.}
  \item{D}{As per input.}
  \item{Delta}{As per input.}
  \item{delta}{As per input.}
  \item{I}{Information vector for the constrained optimal design.}
  \item{I.1}{Information vector of a single treatment effect for the constrained optimal design.}
  \item{I.exact}{Information vector for the exact optimal design.}
  \item{I.1.exact}{Information vector of a single treatment effect for the exact design.}
  \item{init}{Starting values of the parameters used in the searches.}
  \item{init.file.ext}{As per input.}
  \item{L}{As per input.}
  \item{Lambda}{Covariance matrix of the standardised treatment effects for the identified design.}
  \item{lower}{Lower limits for the parameters used in the searches.}
  \item{max.calls}{As per input.}
  \item{n}{Number of patients required in each stage by the constrained optimal design.}
  \item{n.exact}{Number of patients required in each stage by the exact optimal design.}
  \item{num.runs}{As per input.}
  \item{opt.crit}{As per input.}
  \item{penalty}{As per input.}
  \item{scenarios.optimal}{The different possible outcomes for the trial and their respective probabilities for the constrained optimal design.}
  \item{scenarios.exact}{The different possible outcomes for the trial and their respective probabilities for the exact optimal design.}
  \item{seed}{As per input.}
  \item{sequence.type}{As per input.}
  \item{sigma.e}{As per input.}
  \item{upper}{Upper limits for the parameters used in the searches.}
  \item{w}{As per input.}
}
\author{Michael J. Grayling}
\seealso{\code{\link{plot.gsBE}}, \code{\link{simulate.gs.be}}}
\examples{# Compute an example optimal design:
example <- optimal.gs.be()
}
