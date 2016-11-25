plot.gsBE <- function(design, theta = NULL, to.plot = "EN", col = "black",
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
    theta <- seq(from = -3*design$B, to = 3*design$B, length.out = 200)
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

  covariance <- function(D, L){
    Lambda <- matrix(0, 2*(D - 1)*L , 2*(D - 1)*L)
    for (i in 1:L){
      Lambda[(1 + (i - 1)*2*(D - 1)):(i*2*(D - 1)),
             (1 + (i - 1)*2*(D - 1)):(i*2*(D - 1))] <- matrix(0.5, 2*(D - 1), 2*(D - 1))
    }
    for (i in 1:((D - 1)*L)){
      Lambda[(1 + (i - 1)*2):(i*2),
             (1 + (i - 1)*2):(i*2)] <- Lambda[(1 + (i - 1)*2):(i*2),
                                              (1 + (i - 1)*2):(i*2)] + matrix(0.5, 2, 2)
    }
    if (L > 1){
      for (i in 2:L){
        for (j in 1:(i - 1)){
          Cov.ij <- matrix(0.5*sqrt(j/i), 2*(D - 1), 2*(D - 1))
          for (k in 1:(D - 1)){
            Cov.ij[(1 + 2*(k - 1)):(2*k), (1 + 2*(k - 1)):(2*k)] <- Cov.ij[(1 + 2*(k - 1)):(2*k), (1 + 2*(k - 1)):(2*k)] + matrix(0.5*sqrt(j/i), 2, 2)
          }
          Lambda[(1 + (j - 1)*2*(D - 1)):(j*2*(D - 1)),
                 (1 + (i - 1)*2*(D - 1)):(i*2*(D - 1))] <- Cov.ij
          Lambda[(1 + (i - 1)*2*(D - 1)):(i*2*(D - 1)),
                 (1 + (j - 1)*2*(D - 1)):(j*2*(D - 1))] <- Cov.ij
        }
      }
    }
    return(Lambda)
  }

  operatingCharacteristics <- function(tau, D, L, B, a, I, Lambda, scenarios){
    ind.lower.minus <- numeric(L)
    ind.upper.minus <- numeric(L)
    ind.lower.plus  <- numeric(L)
    ind.upper.plus  <- numeric(L)
    lower           <- numeric(2*(D - 1)*L)
    upper           <- numeric(2*(D - 1)*L)
    for (i in 1:nrow(scenarios)){
      for (j in 1:(D - 1)){
        if (scenarios[i, j] == 1){
          if (scenarios[i, D - 1 + j] == 0){
            ind.lower.plus <- rep(-Inf, L)
            ind.upper.plus <- c(a[1], rep(Inf, L - 1))
          } else {
            ind.lower.plus <- c(a[1], rep(-Inf, L - 1))
            ind.upper.plus <- rep(Inf, L)
          }
          if (scenarios[i, 2*(D - 1) + j] == 0){
            ind.lower.minus <- c(-a[1], rep(-Inf, L - 1))
            ind.upper.minus <- rep(Inf, L)
          } else {
            ind.lower.minus <- rep(-Inf, L)
            ind.upper.minus <- c(-a[1], rep(Inf, L - 1))
          }
        } else if (scenarios[i, j] > 1 & scenarios[i, j] < L){
          if (scenarios[i, D - 1 + j] == 0){
            ind.lower.plus <- c(a[1:(scenarios[i, j] - 1)], rep(-Inf, L - (scenarios[i, j] - 1)))
            ind.upper.plus <- c(rep(Inf, scenarios[i, j] - 1), a[scenarios[i, j]],
                                rep(Inf, L - scenarios[i, j]))
          } else {
            ind.lower.plus <- c(a[1:scenarios[i, j]], rep(-Inf, L - scenarios[i, j]))
            ind.upper.plus <- rep(Inf, L)
          }
          if (scenarios[i, 2*(D - 1) + j] == 0){
            ind.lower.minus <- c(rep(-Inf, scenarios[i, j] - 1), -a[scenarios[i, j]],
                                 rep(-Inf, L - scenarios[i, j]))
            ind.upper.minus <- c(-a[1:(scenarios[i, j] - 1)], rep(Inf, L - (scenarios[i, j] - 1)))
          } else {
            ind.lower.minus <- rep(-Inf, L)
            ind.upper.minus <- c(-a[1:scenarios[i, j]], rep(Inf, L - scenarios[i, j]))
          }
        } else {
          if (scenarios[i, D - 1 + j] == 0){
            ind.lower.plus <- c(a[1:(L - 1)], -Inf)
            ind.upper.plus <- c(rep(Inf, L - 1), a[L])
          } else {
            ind.lower.plus <- a
            ind.upper.plus <- rep(Inf, L)
          }
          if (scenarios[i, 2*(D - 1) + j] == 0){
            ind.lower.minus <- c(rep(-Inf, L - 1), -a[L])
            ind.upper.minus <- c(-a[1:(L - 1)], Inf)
          } else {
            ind.lower.minus <- rep(-Inf, L)
            ind.upper.minus <- -a
          }
        }
        lower[seq(from = 1 + 2*(j - 1), by = 2*(D - 1),
                  length.out = L)] <- ind.lower.minus
        lower[seq(from = 2 + 2*(j - 1), by = 2*(D - 1),
                  length.out = L)] <- ind.lower.plus
        upper[seq(from = 1 + 2*(j - 1), by = 2*(D - 1),
                  length.out = L)] <- ind.upper.minus
        upper[seq(from = 2 + 2*(j - 1), by = 2*(D - 1),
                  length.out = L)] <- ind.upper.plus
      }
      scenarios[i, 3*(D - 1) + 6] <- pmvnorm(lower = lower,
                                             upper = upper,
                                             mean = (rep(tau, L) + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                             sigma = Lambda)[1]
    }
    power.any.tau <- sum(scenarios[, 3*(D - 1) + 2]*scenarios[, 3*(D - 1) + 6])
    power.one.tau <- sum(scenarios[, 3*(D - 1) + 3]*scenarios[, 3*(D - 1) + 6])
    if (power.any.tau > 1){
      power.any.tau <- 1
    }
    if (power.one.tau > 1){
      power.one.tau <- 1
    }
    EN.tau  <- sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                      scenarios[, 3*(D - 1) + 6])
    EO.tau  <- sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                      scenarios[, 3*(D - 1) + 6])

    return(c(power.one.tau, power.any.tau, EN.tau, EO.tau))
  }

  wrapper <- function(i){
    performance <- operatingCharacteristics(rep(theta[i], D - 1), D, L, B, a,
                                            I, Lambda, scenarios)
    return(performance)
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  design.performance      <- matrix(0, nrow = length(theta), ncol = 5)
  design.performance[, 1] <- theta
  colnames(design.performance) <- c("theta", "P(RH01|theta)", "P(RH0d|theta)",
                                    "E(N|theta)", "E(O|theta)")

  D           <- design$D
  L           <- design$L
  Lambda      <- design$Lambda
  B           <- design$B
  if (exact.design == TRUE){
    I         <- design$I.exact
    a         <- design$a.exact
    n         <- design$n.exact
    scenarios <- design$scenarios.exact
  } else {
    I         <- design$I
    a         <- design$a
    n         <- design$n
    scenarios <- design$scenarios
  }

  if (summary == TRUE){
    print("Determining performance across theta...")
  }

  sink("NULL")
  suppressMessages(sfInit(parallel = parallel, cpus = cpus))
  suppressMessages(sfLibrary(mvtnorm))
  sink()
  sfExport("D", "L", "B", "a", "I", "Lambda", "scenarios", "theta",
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
            ylab = expression(paste("P(Reject ", H[0]^{(d)}, " for some d | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")), ...)
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
      plot(theta, design.performance[, 2], col = col, type = type, lty = lty, lwd = lwd,
           xlim = xlim, ylim = ylim, xlab = expression(theta),
           ylab = expression(paste("P(Reject ", H[0]^1, " | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")), ...)
    } else if (to.plot == "Power.Any"){
      if (is.null(ylim)){
        ylim = c(0, 1)
      }
      plot(theta, design.performance[, 3], col = col, type = type, lty = lty, lwd = lwd,
           xlim = xlim, ylim = ylim, xlab = expression(theta),
           ylab = expression(paste("P(Reject ", H[0]^{(d)}, " for some d | ", tau[1], " = ... = ", tau[D - 1], " = ", theta, ")", sep = "")), ...)
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
