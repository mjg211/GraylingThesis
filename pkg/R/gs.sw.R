gs.sw <- function(X, sigma.e = sqrt(0.5), sigma.c = sqrt(0.02), alpha = 0.05,
                  beta = 0.1, delta = 0.2, set.T = 2:5, gamma = 0.5,
                  summary = TRUE){

  ##### ERROR CHECKING ########################################################

  if (!is.matrix(X) | ncol(X) == 1 | nrow(X) == 1 | !all(X %in% c(0, 1))){
    stop("X must be an indicator matrix with at least 2 rows and 2 columns.")
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
  if (delta <= 0){
    stop("Clinically relevant difference delta to power for must be strictly positive.")
  }
  if (!is.vector(set.T) | !all(set.T %in% 1:ncol(X)) | sum(X[, 1:set.T[1]]) == 0){
    stop("set.T must be a vector with elements in (min_j sum(X[, 1:j]) > 0):ncol(X).")
  }
  if (gamma < 0){
    stop("gamma must be greater than or equal to 0.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }

  ##### FUNCTION INITIALISATION ###############################################

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

  boundaryF <- function(crit, pi2l, prefbounds, preebounds, currSigma, currI,
                        delta){
    integral <- pmvnorm(lower = c(prefbounds, -Inf),
                        upper = c(preebounds, crit),
                        mean = rep(delta, length(currI))*sqrt(currI),
                        sigma = currSigma)[1]
    return((pi2l - integral)^2)
  }

  boundaryE <- function(crit, pi1l, prefbounds, preebounds, currSigma){
    integral <- pmvnorm(lower = c(prefbounds, crit),
                        upper = c(preebounds, Inf),
                        mean = numeric(length(prefbounds) + 1),
                        sigma = currSigma)[1]
    return((pi1l - integral)^2)
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

  covariance <- function(L, I){
    Sigma <- diag(1, L, L)
    for (l1 in 2:L){
      for (l2 in 1:(l1 - 1)){
        Sigma[l1, l2] <- sqrt(I[l2]/I[l1])
        Sigma[l2, l1] <- Sigma[l1, l2]
      }
    }
    return(Sigma)
  }

  family <- function(L, I, spend, gamma){
    pi      <- numeric(L)
    pi[1]   <- min(spend*(I[1]/I[L])^gamma, spend)
    for (l in 2:L){
      pi[l] <- min(spend*(I[l]/I[L])^gamma, spend) -
        min(spend*(I[l - 1]/I[L])^gamma, spend)
    }
    return(pi)
  }

  nSeqCRCTF <- function(n, X, delta, alpha, beta, sigma.e, sigma.c, set.T, gamma){
    C     <- nrow(X)
    L     <- length(set.T)
    I     <- informationCRCT(n, X, set.T, sigma.e, sigma.c)
    Sigma <- covariance(L, I)
    pi2   <- family(L, I, beta, gamma)
    f     <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        f[l] <- qnorm(pi2[l], mean = delta*sqrt(I[l]))
      } else if ((l > 1) & (l < L)){
        f[l] <- suppressWarnings(optim(par = qnorm(1 - pi2[l],
                                                   mean = delta*sqrt(I[l])),
                                       fn = boundaryF, pi2l = pi2[l],
                                       prefbounds = f[1:(l - 1)],
                                       preebounds = rep(Inf, l - 1),
                                       currSigma = Sigma[1:l, 1:l],
                                       currI = I[1:l], delta = delta)$par)
      } else {
        f[l] <- suppressWarnings(optim(par = qnorm(1 - alpha/L),
                                       fn = boundaryE, pi1l = alpha,
                                       prefbounds = f[1:(l - 1)],
                                       preebounds = rep(Inf, L - 1),
                                       currSigma = Sigma)$par)
      }
    }
    nvec <- n*C*set.T
    PR   <- perfF(delta, L, f, I, nvec, Sigma)[1]
    return((beta - (1 - PR))^2)
  }

  nCRCT <- function(n, X, delta, beta, sigma.e, sigma.c, e){
    I     <- informationCRCT(n, X, ncol(X), sigma.e, sigma.c)
    PnotR <- pnorm(q = e, mean = delta*sqrt(I))
    Score <- (beta - PnotR)^2
    return(Score)
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  n.sw     <- suppressWarnings(optim(par = 10^-2, fn = nCRCT, X = X,
                                     delta = delta, beta = beta,
                                     sigma.e = sigma.e, sigma.c = sigma.c,
                                     e = qnorm(1 - alpha))$par)

  if (summary == TRUE){
    print("Determining exact GS design...")
  }

  n.seq.sw <- suppressWarnings(optim(par = n.sw, fn = nSeqCRCTF, X = X,
                                     delta = delta, alpha = alpha,
                                     beta = beta, sigma.e = sigma.e,
                                     sigma.c = sigma.c, set.T = set.T,
                                     gamma = gamma)$par)
  C        <- nrow(X)
  L        <- length(set.T)
  n        <- ceiling(n.seq.sw)
  nvec     <- n*C*set.T
  I        <- informationCRCT(n, X, set.T, sigma.e, sigma.c)
  Lambda   <- covariance(L, I)
  pi2      <- family(L, I, beta, gamma)
  f        <- numeric(L)
  for (l in 1:L){
    if (l == 1){
      f[l] <- qnorm(pi2[l], mean = delta*sqrt(I[l]))
    } else if ((l > 1) & (l < L)){
      f[l] <- suppressWarnings(optim(par = qnorm(1 - pi2[l],
                                                 mean = delta*sqrt(I[l])),
                                     fn = boundaryF, pi2l = pi2[l],
                                     prefbounds = f[1:(l - 1)],
                                     preebounds = rep(Inf, l - 1),
                                     currSigma = Lambda[1:l, 1:l],
                                     currI = I[1:l], delta = delta)$par)
    } else {
      f[l] <- suppressWarnings(optim(par = qnorm(1 - alpha/L), fn = boundaryE,
                                     pi1l = alpha, prefbounds = f[1:(l - 1)],
                                     preebounds = rep(Inf, L - 1),
                                     currSigma = Lambda)$par)
    }
  }
  perf.H0        <- perfF(0, L, f, I, nvec, Lambda)
  perf.H1        <- perfF(delta, L, f, I, nvec, Lambda)

  if (summary == TRUE){
    print("Outputting...")
  }

  output <- list(alpha = alpha, beta = beta, delta = delta, f = f,
                 gamma = gamma, I = I, Lambda = Lambda, n = n,
                 perf.H0 = perf.H0, perf.H1 = perf.H1, set.T = set.T,
                 sigma.c = sigma.c, sigma.e = sigma.e, X = X)
  class(output) <- "gsSW"
  return(output)
}
