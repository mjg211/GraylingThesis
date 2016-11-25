plot.gsSW <- function(design, theta = NULL, to.plot = "EN", col = "black",
                      add = FALSE, type = "l", lty = 1, lwd = 2, xlim = NULL,
                      ylim = NULL, return.perf = FALSE, summary = TRUE, ...){

  ##### ERROR CHECKING ########################################################

  if (!is.null(theta)){
    if (!is.vector(theta)){
      stop("theta must be a vector of strictly increasing values.")
    }
    for (i in 1:(length(theta) - 1)){
      if (theta[i] >= theta[i + 1]){
        stop("theta must be a vector of strictly increasing values.")
      }
    }
  } else {
    theta <- seq(from = -3*design$delta, to = 3*design$delta, length.out = 200)
  }
  if (!(to.plot %in% c("Power", "EN"))){
    stop("to.plot must be set to one of \"Power\" or \"EN\".")
  }
  if (!is.logical(add)){
    stop("add must be set to TRUE or FALSE.")
  }
  if (!is.logical(return.perf)){
    stop("return.perf must be set to TRUE or FALSE.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }

  ##### FUNCTION INITIALISATION ###############################################

  perfEF <- function(theta, L, e, f, I, nvec, Sigma){
    PE <- numeric(L)
    PF <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        PE[l] <- pmvnorm(lower = e[l], upper = Inf,
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
        PF[l] <- pmvnorm(lower = -Inf, upper = f[l],
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
      } else {
        PE[l] <- pmvnorm(lower = c(f[1:(l - 1)], e[l]),
                         upper = c(e[1:(l - 1)], Inf),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
        PF[l] <- pmvnorm(lower = c(f[1:(l - 1)], -Inf),
                         upper = c(e[1:(l - 1)], f[l]),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
      }
    }
    PR <- sum(PE)
    EN <- sum(nvec*(PE + PF))
    return(c(PR, EN))
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  design.performance      <- matrix(0, nrow = length(theta), ncol = 3)
  design.performance[, 1] <- theta
  colnames(design.performance) <- c("theta", "P(RH0|theta)", "E(N|theta)")

  L           <- length(design$set.T)
  Lambda      <- design$Lambda
  n           <- design$n
  f        <- design$f
  e        <- design$e
  I           <- design$I
  nvec      <- design$n*design$set.T*nrow(design$X)

  if (summary == TRUE){
    print("Determining performance across theta...")
  }

  for (i in 1:length(theta)){
    design.performance[i, 2:3] <- perfEF(theta[i], L, e, f, I, nvec, Lambda)
  }

  if (summary == TRUE){
    print("Plotting design performance...")
  }

  if (add == TRUE){
    if (to.plot == "Power"){
      lines(theta, design.performance[, 2], col = col, type = type, lty = lty, lwd = lwd,
            xlab = expression(theta),
            ylab = expression(paste("P(Reject ", H[0], " | ", theta, ")", sep = "")), ...)
    } else if (to.plot == "EN"){
      lines(theta, design.performance[, 3], col = col, type = type, lty = lty, lwd = lwd,
            xlab = expression(tau),
            ylab = expression(paste("E(N | ", theta, ")", sep = "")), ...)
    }
  } else {
    par(mar = c(5.5, 4.5, 4.5, 2.5))
    if (is.null(xlim)){
      xlim = c(theta[1], theta[length(theta)])
    }
    if (to.plot == "Power"){
      if (is.null(ylim)){
        ylim = c(0, 1)
      }
      plot(theta, design.performance[, 2], col = col, type = type, lty = lty, lwd = lwd,
           xlim = xlim, ylim = ylim, xlab = expression(theta),
           ylab = expression(paste("P(Reject ", H[0], " | ", theta, ")", sep = "")), ...)
    } else if (to.plot == "EN"){
      if (is.null(ylim)){
        ylim = c(0.9*design$n*nrow(design$X)*design$set.T[1],
                 1.1*design$n*nrow(design$X)*ncol(design$X))
      }
      plot(theta, design.performance[, 3], col = col, type = type, lty = lty, lwd = lwd,
           xlim = xlim, ylim = ylim, xlab = expression(theta),
           ylab = expression(paste("E(N | ", theta, ")", sep = "")))
    }
  }
  if (return.perf == TRUE){
    return(design.performance)
  }
}
