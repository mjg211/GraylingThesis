\name{plot.gsBE}
\alias{plot.gsBE}
\title{S3 Plotting Method for Group Sequential Bioequivalence Trial Designs}
\description{Plots the performance of a one-sided group sequential crossover trial.}
\usage{plot(design, theta = NULL, to.plot = "EN", col = "black",
     add = FALSE, type = "l", lty = 1, lwd = 2, xlim = NULL,
     ylim = NULL, return.perf = FALSE, exact.design = FALSE,
     parallel = TRUE, cpus = 8, summary = TRUE, ...)}
\arguments{
  \item{design}{A design as returned by gs.be or optimal.gs.be}
  \item{theta}{Range of values for theta; plotted across the x-axis. Defaults (internally) to seq(from = -3*design$B, to = 3*design$B, length.out = 200).}
  \item{to.plot}{Determines what is plotted on the y-axis. Must be set to "Power.1", "Power.Any", "EN" or "EO". Defaults to "EN".}
  \item{col}{Colour to use for the plotted curve. Defaults to "black".}
  \item{add}{Detemines whether to add the curve to a new plot or on to an existing one. Must be logical. Defaults to FALSE.}
  \item{type}{The line type. Passed to plot. Defaults to "l".}
  \item{lty}{The line type. Passed to plot. Defaults to 1.}
  \item{lwd}{The line width. Passed to plot. Defaults to 2.}
  \item{xlim}{The limits of the x-axis of the plot. Either a vector of length 2, or NULL in which case appropriate values are chosen by the program. Defaults to NULL.}
  \item{ylim}{The limits of the y-axis of the plot. Either a vector of length 2, or NULL in which case appropriate values are chosen by the program. Defaults to NULL.}
  \item{return.perf}{Specifies whether to return a matrix containing the performance of the design across
the specified range for theta. Must be logical. Defaults to FALSE.}
  \item{exact.design}{Specifies whether to plot the performance of the exact or rounded design. Must be logical. Defaults to FALSE.}
  \item{parallel}{Specifies whether to use parallel computation to speed up the program. Must be logical. Defaults to TRUE.}
  \item{cpus}{Specifies the number of CPUs to parallelise over if parallel is set to TRUE. Must be an integer greater than or equal to 1. Defaults to 8.}
  \item{summary}{Logical. Indicates whether to print a summary of the programs progress. Defaults to TRUE.}
  \item{\dots}{Additional arguments to pass to plot.}
}
\value{If return.perf is set to TRUE, a matrix is returned containing the performance of the design.}
\author{Michael J. Grayling}
\seealso{\code{\link{optimal.gs.be}}, \code{\link{plot}}}
\examples{# Compute an example optimal design:
example <- optimal.gs.be()
# Plot its expected sample size curve:
plot(example)
}
