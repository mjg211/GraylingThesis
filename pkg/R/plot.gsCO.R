plot.gsCO <- function(design, theta = NULL, to.plot = "EN", col = "black",
                      add = FALSE, type = "l", lty = 1, lwd = 2, xlim = NULL,
                      ylim = NULL, return.perf = FALSE, exact.design = FALSE,
                      parallel = TRUE, cpus = 8, summary = TRUE, ...){

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
    theta <- seq(from = -3, to = 5, by = 0.01)
  }
  if (!(to.plot %in% c("Power.1", "Power.Any", "EN", "EO"))){
    stop("to.plot must be set to one of \"Power.1\", \"Power.Any\", \"EN\" or \"EO\".")
  }
  if (!is.logical(add)){
    stop("add must be set to TRUE or FALSE.")
  }
  if (!is.logical(return.perf)){
    stop("return.perf must be set to TRUE or FALSE.")
  }
  if (!is.logical(exact.design)){
    stop("exact.design must be set to TRUE or FALSE.")
  }
  if (!is.logical(parallel)){
    stop("parallel must be set to TRUE or FALSE.")
  }
  if ((cpus%%1 != 0) | (cpus < 1)){
    stop("cpus must be a whole number greater than or equal to 1.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }

  ##### FUNCTION INITIALISATION ###############################################

  covarianceCO <- function(D, L) {
    Lambda <- matrix(0, (D - 1)*L, (D - 1)*L)
    for (i in 1:L){
      Lambda[(1 + (i - 1)*(D - 1)):(i*(D - 1)),
             (1 + (i - 1)*(D - 1)):(i*(D - 1))] <- matrix(0.5, D - 1, D - 1) +
        diag(0.5, D - 1, D - 1)
    }
    if (L > 1){
      for (i in 2:L){
        for (j in 1:(i - 1)){
          Cov.ij <- diag(0.5*sqrt(j/i), D - 1, D - 1) +
            matrix(0.5*sqrt(j/i), D - 1, D - 1)
          Lambda[(1 + (j - 1)*(D - 1)):(j*(D - 1)),
                 (1 + (i - 1)*(D - 1)):(i*(D - 1))] <- Cov.ij
          Lambda[(1 + (i - 1)*(D - 1)):(i*(D - 1)),
                 (1 + (j - 1)*(D - 1)):(j*(D - 1))] <- Cov.ij
        }
      }
    }
    return(Lambda)
  }

  operatingCharacteristics <- function(tau, D, L, e, f, I, Lambda, scenarios){
    ind.lower <- numeric(L)
    ind.upper <- numeric(L)
    lower     <- numeric((D - 1)*L)
    upper     <- numeric((D - 1)*L)
    for (i in 1:nrow(scenarios)){
      for (j in 1:(D - 1)){
        if ((scenarios[i, j] == 1) & (scenarios[i, D - 1 + j] == 0)){
          ind.lower <- rep(-Inf, L)
          ind.upper <- c(f[1], rep(Inf, L - 1))
        } else if ((scenarios[i, j] == 1) & (scenarios[i, D - 1 + j] == 1)){
          ind.lower <- c(e[1], rep(-Inf, L - 1))
          ind.upper <- rep(Inf, L)
        } else if ((scenarios[i, j] %in% 2:(L - 1)) &
                   (scenarios[i, D - 1 + j] == 0)){
          ind.lower <- c(f[1:(scenarios[i, j] - 1)], -Inf,
                         rep(-Inf, L - scenarios[i, j]))
          ind.upper <- c(e[1:(scenarios[i, j] - 1)], f[scenarios[i, j]],
                         rep(Inf, L - scenarios[i, j]))
        } else if ((scenarios[i, j] %in% 2:(L - 1)) &
                   (scenarios[i, D - 1 + j] == 1)){
          ind.lower <- c(f[1:(scenarios[i, j] - 1)], e[scenarios[i, j]],
                         rep(-Inf, L - scenarios[i, j]))
          ind.upper <- c(e[1:(scenarios[i, j] - 1)], Inf,
                         rep(Inf, L - scenarios[i, j]))
        } else if ((scenarios[i, j] == L) & (scenarios[i, D - 1 + j] == 0)){
          ind.lower <- c(f[1:(scenarios[i, j] - 1)], -Inf)
          ind.upper <- e
        } else if ((scenarios[i, j] == L) & (scenarios[i, D - 1 + j] == 1)){
          ind.lower <- f
          ind.upper <- c(e[1:(scenarios[i, j] - 1)], Inf)
        }
        lower[seq(from = j, to = j + (D - 1)*(L - 1),
                  length.out = L)] <- ind.lower
        upper[seq(from = j, to = j + (D - 1)*(L - 1),
                  length.out = L)] <- ind.upper
      }
      scenarios[i, 2*(D - 1) + 6] <- pmvnorm(lower = lower,
                                             upper = upper,
                                             mean = rep(tau, L)*sqrt(I),
                                             sigma = Lambda)[1]
    }
    power.any.tau <- sum(scenarios[, 2*(D - 1) + 2]*scenarios[, 2*(D - 1) + 6])
    if (power.any.tau > 1){
      power.any.tau <- 1
    }
    power.one.tau <- sum(scenarios[, 2*(D - 1) + 3]*scenarios[, 2*(D - 1) + 6])
    if (power.one.tau > 1){
      power.one.tau <- 1
    }
    EN.tau        <- sum(scenarios[, 2*(D - 1) + 1]*
                           scenarios[, 2*(D - 1) + 6]*
                           scenarios[, 2*(D - 1) + 4])
    EO.tau        <- sum(scenarios[, 2*(D - 1) + 1]*
                           scenarios[, 2*(D - 1) + 6]*
                           scenarios[, 2*(D - 1) + 5])
    return(c(power.one.tau, power.any.tau, EN.tau, EO.tau))
  }

  wrapper <- function(i){
    performance <- operatingCharacteristics(rep(theta[i], D - 1), D, L, e, f,
                                            I, Lambda, scenarios)
    return(performance)
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  design.performance           <- matrix(0, nrow = length(theta), ncol = 5)
  design.performance[, 1]      <- theta
  colnames(design.performance) <- c("theta", "P(RH01|theta)", "P(RH0d|theta)",
                                    "E(N|theta)", "E(O|theta)")

  D           <- design$D
  L           <- design$L
  e           <- design$e
  f           <- design$f
  Lambda      <- design$Lambda
  if (exact.design == FALSE){
    I         <- design$I
    n         <- design$n
    Lambda    <- design$Lambda
    scenarios <- design$scenarios
  } else {
    I         <- design$I.exact
    n         <- design$n.exact
    scenarios <- design$scenarios.exact
  }

  if (summary == TRUE){
    print("Determining performance across theta...")
  }

  sink("NULL")
  suppressMessages(sfInit(parallel = parallel, cpus = cpus))
  suppressMessages(sfLibrary(mvtnorm))
  sink()
  sfExport("D", "L", "e", "f", "I", "Lambda", "scenarios", "theta",
           "operatingCharacteristics")
  results <- sfLapply(1:length(theta), wrapper)
  for (i in 1:length(theta)){
    design.performance[i, 2:5] <- results[[i]]
  }

  if (summary == TRUE){
    print("Plotting design performance...")
  }

  if (add == TRUE){
    if (to.plot == "Power.1"){
      lines(theta, design.performance[, 2], col = col, type = type, lty = lty, lwd = lwd,
            xlab = expression(theta),
            ylab = expression(paste("P(Reject ", H[0]^1, " | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")), ...)
    } else if (to.plot == "Power.Any"){
      lines(theta, design.performance[, 3], col = col, type = type, lty = lty, lwd = lwd,
            xlab = expression(theta),
            ylab = expression(paste("P(Reject ", H[0]^i, " for some i | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")), ...)
    } else if (to.plot == "EN"){
      lines(theta, design.performance[, 4], col = col, type = type, lty = lty, lwd = lwd,
            xlab = expression(tau),
            ylab = expression(paste("E(N | ", tau[1], " = ... = ", tau[D - 1], " = ", tau, ")", sep = "")), ...)
    } else if (to.plot == "EO"){
      lines(theta, design.performance[, 5], col = col, type = type, lty = lty, lwd = lwd,
            xlab = expression(tau),
            ylab = expression(paste("E(O | ", tau[1], " = ... = ", tau[D - 1], " = ", tau, ")", sep = "")), ...)
    }
  } else {
    par(mar = c(5.5, 4.5, 4.5, 2.5))
    if (is.null(xlim)){
      xlim = c(theta[1], theta[length(theta)])
    }
    if (to.plot == "Power.1"){
      if (is.null(ylim)){
        ylim = c(0, 1)
      }
      plot(theta, design.performance[, 3], col = col, type = type, lty = lty, lwd = lwd,
           xlim = xlim, ylim = ylim, xlab = expression(theta),
           ylab = expression(paste("P(Reject ", H[0]^1, " | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")), ...)
    } else if (to.plot == "Power.Any"){
      if (is.null(ylim)){
        ylim = c(0, 1)
      }
      plot(theta, design.performance[, 2], col = col, type = type, lty = lty, lwd = lwd,
           xlim = xlim, ylim = ylim, xlab = expression(theta),
           ylab = expression(paste("P(Reject ", H[0]^i, " for some i | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")), ...)
    } else if (to.plot == "EN"){
      if (is.null(ylim)){
        ylim = c(0.9*n, 1.1*n*L)
      }
      plot(theta, design.performance[, 4], col = col, type = type, lty = lty, lwd = lwd,
           xlim = xlim, ylim = ylim, xlab = expression(theta),
           ylab = expression(paste("E(N | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")))
    } else if (to.plot == "EO"){
      if (is.null(ylim)){
        ylim = c(n*D, n*D*L)
      }
      plot(theta, design.performance[, 5], col = col, type = type, lty = lty, lwd = lwd,
           xlim = xlim, ylim = ylim, xlab = expression(theta),
           ylab = expression(paste("E(O | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")), ...)
    }
  }
  if (return.perf == TRUE){
    return(design.performance)
  }
}
