optimal.gs.sw <- function(C = 4, Ti = 5, sigma.c = sqrt(0.02),
                          sigma.e = sqrt(0.5), alpha = 0.05, beta = 0.1,
                          delta = 0.2, set.T = 2:5, eta = -1, penalty = 1,
                          seed = Sys.time(), summary = TRUE){

  start.time <- Sys.time()

  ##### ERROR CHECKING ########################################################

  if ((C%%1 != 0) | (C < 2)){
    stop("C must be a whole number greater than or equal to 2.")
  }
  if ((Ti%%1 != 0) | (Ti < 2)){
    stop("Ti must be a whole number greater than or equal to 2.")
  }
  if (sigma.e <= 0){
    stop("Within person standard deviation sigma.e must be strictly positive.")
  }
  if (sigma.c <= 0){
    stop("Between cluster standard deviation sigma.e must be strictly positive.")
  }
  if ((alpha <= 0) | (alpha >= 1)){
    stop("Type-I error rate alpha must be strictly between 0 and 1.")
  }
  if ((beta <= 0) | (beta >= 1)){
    stop("Type-II error rate beta must be strictly between 0 and 1.")
  }
  if (!is.numeric(eta)){
    stop("eta must be numeric.")
  }
  if (penalty <= 0){
    stop("penalty must be strictly positive.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }

  ##### FUNCTION INITIALISATION ###############################################

  boundariesScore <- function(n1.f, C, Ti, switches, delta, alpha, beta, sigma.c,
                              sigma.e, set.T, eta, penalty){
    if (sum(sort(switches) == switches) != C){
      Score                 <- Inf
    } else if (min(switches) > set.T[1]){
      Score                 <- Inf
    } else if (length(unique(switches)) == 1){
      Score                 <- Inf
    } else {
      X                     <- matrix(0, nrow = C, ncol = Ti)
      for (c in 1:C){
        X[c, switches[c]:Ti] <- 1
      }
      n1                    <- n1.f[1]
      f                     <- n1.f[2:length(n1.f)]
      n                     <- n1*C*set.T
      sigma2                <- sigma.e^2/n1
      L                     <- length(set.T)
      I                     <- numeric(L)
      for (i in 1:length(set.T)){
        U    <- sum(X[, 1:set.T[i]])
        if (set.T[i] > 1){
          V  <- sum(rowSums(X[, 1:set.T[i]])^2)
          W  <- sum(colSums(X[, 1:set.T[i]])^2)
        } else {
          V  <- sum(X[, set.T[i]]^2)
          W  <- sum(X[, set.T[i]])^2
        }
        I[i] <- ((C*U - W)*sigma2 + (U^2 + C*set.T[i]*U - set.T[i]*W - C*V)*sigma.c^2)/
          (C*sigma2*(sigma2 + set.T[i]*sigma.c^2))
      }
      Sigma <- diag(1, L, L)
      for (i in 2:L){
        for (j in 1:(i - 1)){
          Sigma[i, j] <- sqrt(I[j]/I[i])
          Sigma[j, i] <- Sigma[i, j]
        }
      }
      if (any(as.vector(is.nan(Sigma)))){
        Score <- Inf
      } else {
        P.rej.H0        <- pmvnorm(lower = f, upper = rep(Inf, L),
                                   mean = numeric(L), sigma = Sigma)[1]
        P.n.rej.H1      <- 1 - pmvnorm(lower = f, upper = rep(Inf, L),
                                       mean = delta*sqrt(I), sigma = Sigma)[1]
        P.F.delta.eta <- numeric(L)
        for (l in 1:L){
          if (l == 1){
            P.F.delta.eta[l] <- pmvnorm(lower = -Inf, upper = f[l],
                                        mean = delta*eta*sqrt(I[l]),
                                        sigma = 1)[1]
          } else {
            P.F.delta.eta[l] <- pmvnorm(lower = c(f[1:(l - 1)], -Inf),
                                        upper = c(rep(Inf, l - 1), f[l]),
                                        mean = delta*eta*sqrt(I[1:l]),
                                        sigma = Sigma[1:l, 1:l])[1]
          }
        }
        EN.delta.eta <- sum(n[1:(L - 1)]*P.F.delta.eta[1:(L - 1)]) +
          n[L]*(1 - sum(P.F.delta.eta[1:(L - 1)]))
        Score        <- EN.delta.eta + penalty*(as.numeric(P.rej.H0 > alpha)*
                                                  (P.rej.H0 - alpha)/alpha +
                                                  as.numeric(P.n.rej.H1 > beta)*
                                                  (P.n.rej.H1 - beta)/beta)
      }
    }
    return(Score)
  }

  nCRCT <- function(n, X, delta, beta, sigma.e, sigma.c, e){
    I     <- informationCRCT(n, X, ncol(X), sigma.e, sigma.c)
    PnotR <- pnorm(q = e, mean = delta*sqrt(I))
    Score <- (beta - PnotR)^2
    return(Score)
  }

  informationCRCT <- function(n, X, set.T, sigma.e, sigma.c){
    if (length(set.T) == 1 && set.T == 1){
      C      <- length(X)
      sigma2 <- sigma.e^2/n
      U      <- sum(X)
      V      <- sum(X^2)
      W      <- sum(X)^2
      I      <- ((C*U - W)*sigma2 + (U^2 + C*U - W - C*V)*sigma.c^2)/
        (C*sigma2*(sigma2 + sigma.c^2))
    } else {
      C      <- nrow(X)
      sigma2 <- sigma.e^2/n
      I      <- numeric(length(set.T))
      for (i in 1:length(set.T)){
        U    <- sum(X[, 1:set.T[i]])
        if (set.T[i] > 1){
          V  <- sum(rowSums(X[, 1:set.T[i]])^2)
          W  <- sum(colSums(X[, 1:set.T[i]])^2)
        } else {
          V  <- sum(X[, set.T[i]]^2)
          W  <- sum(X[, set.T[i]])^2
        }
        I[i] <- ((C*U - W)*sigma2 + (U^2 + C*set.T[i]*U - set.T[i]*W - C*V)*sigma.c^2)/
          (C*sigma2*(sigma2 + set.T[i]*sigma.c^2))
      }
    }
    return(I)
  }

  perfF <- function(theta, L, f, I, nvec, Sigma){
    PF <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        PF[l] <- pmvnorm(lower = -Inf, upper = f[l],
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
      } else {
        PF[l] <- pmvnorm(lower = c(f[1:(l - 1)], -Inf),
                         upper = c(rep(Inf, l - 1), f[l]),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
      }
    }
    PR <- 1 - sum(PF)
    EN <- sum(nvec[1:(L - 1)]*PF[1:(L - 1)]) +
      nvec[L]*(1 - sum(PF[1:(L - 1)]))
    return(c(PR, EN))
  }

  objFn <- function(n1.f, switches, C, Ti, delta, alpha, beta, sigma.c, sigma.e, set.T, eta, penalty){
    switches.int          <- switches + 2L
    Score                 <- boundariesScore(n1.f, C, Ti, switches.int, delta, alpha, beta,
                                             sigma.c, sigma.e, set.T, eta, penalty)
    return(Score)
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  set.seed(seed)

  switches <- numeric(C)
  for (c in 1:C){
    switches[c] <- (c + 1)%%Ti
    if (switches[c] == 0){
      switches[c] <- Ti
    }
    if (c >= Ti){
      switches[c] <- switches[c] + 1
    }
  }
  e                             <- qnorm(1 - alpha)
  switches.initial              <- switches
  X.initial                     <- matrix(0, nrow = C, ncol = Ti)
  for (c in 1:C){
    X.initial[c, switches.initial[c]:Ti] <- 1
  }
  n.initial <- suppressWarnings(optim(par = 10^-2, fn = nCRCT,
                                      X = X.initial, delta = delta,
                                      sigma.c = sigma.c,
                                      sigma.e = sigma.e, beta = beta,
                                      e = e)$par)
  penalty    <- penalty*n.initial*C*Ti
  contMean     <- c(n.initial, rep(0, length(set.T)))
  contSD       <- c(n.initial, rep(10, length(set.T)))
  contConstMat <- matrix(0, length(set.T) + 1, length(set.T) + 1)
  contConstMat[1, 1] <- -1
  contConstVec <- numeric(length(set.T) + 1)
  discCat      <- as.integer(rep(Ti - 1, C))
  discSmooth   <- 0.5

  if (summary == TRUE){
    print("Searching for optimal GS design...")
  }

  optimal.design <- CEoptim(f = objFn, f.arg = list(C = C, Ti = Ti, delta = delta, alpha = alpha,
                                                      beta = beta, sigma.c = sigma.c, sigma.e = sigma.e,
                                                      set.T = set.T, eta = eta, penalty = penalty),
                            continuous = list(mean = contMean, sd = contSD,
                                              conMat = contConstMat, conVec = contConstVec,
                                              smoothMean = 0.5, smoothSd = 0.5),
                            discrete = list(categories = discCat,
                                            smoothProb = discSmooth), N = 10000L, rho = 0.001,
                            verbose = TRUE)
  score    <- optimal.design$optimum
  switches <- optimal.design$optimizer$discrete + 2
  n        <- ceiling(optimal.design$optimizer$continuous[1])
  f        <- optimal.design$optimizer$continuous[2:(length(set.T) + 1)]

  X                  <- matrix(0, nrow = C, ncol = Ti)
  for (c in 1:C){
    X[c, switches[c]:Ti] <- 1
  }
  nvec <- n*C*set.T
  I    <- informationCRCT(n, X, set.T, sigma.e, sigma.c)
  L <- length(set.T)
  Lambda <- diag(1, L, L)
  for (i in 2:L){
    for (j in 1:(i - 1)){
      Lambda[i, j] <- sqrt(I[j]/I[i])
      Lambda[j, i] <- Lambda[i, j]
    }
  }
  perf.H0        <- perfF(0, length(set.T), f, I, nvec, Lambda)
  perf.H1        <- perfF(delta, length(set.T), f, I, nvec, Lambda)
  perf.delta.eta <- perfF(delta*eta, length(set.T), f, I, nvec, Lambda)

  end.time <- Sys.time()
  run.time <- as.numeric(difftime(end.time, start.time, units = "min"))

  if (summary == TRUE){
    print("Outputting...")
  }

  output        <- list(alpha = alpha, beta = beta, C = C, delta = delta, eta = eta,
                        f = f, I = I, Lambda = Lambda, n = n, penalty = penalty,
                        perf.eta.delta = perf.eta.delta, perf.H0 = perf.H0,
                        perf.H1 = perf.H1, run.time = run.time, seed = seed,
                        set.T = set.T, sigma.c = sigma.c, sigma.e = sigma.e,
                        Ti = Ti, X = X)
  class(output) <- "gsSW"
  return(output)
}
