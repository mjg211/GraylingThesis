optimal.gs.be <- function(algorithm = "DEoptim", D = 3, L = 2,
                          sigma.e = sqrt(log(0.3^2 + 1)), alpha = 0.05,
                          beta = 0.2, B = log(1.25), delta = 0,
                          sequence.type = "williams", w = rep(1/3, 3),
                          opt.crit = "O", penalty = 1, num.runs = 1,
                          find.constrained = TRUE, max.calls = 10000,
                          lower = NULL, upper = NULL, init = NULL,
                          init.file.ext = 1, seed = Sys.time(),
                          summary = TRUE){

  ##### ERROR CHECKING ########################################################

  if (!(algorithm %in% c("DEopt", "DEoptim", "ga", "GenSA", "hydroPSO",
                         "NLOPT_GN_CRS2_LM", "NLOPT_GN_DIRECT",
                         "NLOPT_GN_DIRECT_L", "NLOPT_GN_ISRES", "PSopt",
                         "psoptim", "SANN", "soma"))){
    stop("algorithm must be set to \"DEopt\", \"DEoptim\", \"ga\", \"GenSA\", \"hydroPSO\", \"NLOPT_GN_CRS2_LM\", \"NLOPT_GN_DIRECT\", \"NLOPT_GN_DIRECT_L\", \"NLOPT_GN_ISRES\", \"PSopt\", \"psoptim\", \"SANN\", or \"soma\".")
  }
  if ((D%%1 != 0) | (D < 2)){
    stop("D must be a whole number greater than or equal to 2.")
  }
  if ((L%%1 != 0) | (L < 1)){
    stop("L must be a whole number greater than or equal to 1.")
  }
  if (sigma.e <= 0){
    stop("Within person standard deviation sigma.e must be strictly positive.")
  }
  if ((alpha <= 0) | (alpha >= 1)){
    stop("FWER alpha must be strictly between 0 and 1.")
  }
  if ((beta <= 0) | (beta >= 1)){
    stop("Type-II error rate beta must be strictly between 0 and 1.")
  }
  if (B <= 0){
    stop("BE margin B must be strictly positive.")
  }
  if (delta < 0 | delta >= B){
    stop("BE margin delta to power for must be in [0,B).")
  }
  if (!(sequence.type %in% c("latin", "williams"))){
    stop("sequence.type must be set to \"latin\" or \"williams\".")
  }
  if (!is.vector(w) | length(w) != 3 | any(w < 0) | sum(w[1:2]) == 0){
    stop("w must be a vector of length 3, with all elements greater than or equal to 0, and at least one of the first 2 elements strictly positive.")
  }
  if (!(opt.crit %in% c("N", "O"))){
    stop("opt.crit must be set to \"opt.crit\" or \"opt.crit\".")
  }
  if (penalty <= 0){
    stop("penalty must be strictly positive.")
  }
  if ((num.runs%%1 != 0) | (num.runs < 1)){
    stop("num.runs must be a whole number greater than or equal to 1.")
  }
  if (!is.logical(find.constrained)){
    stop("find.constrained must be set to TRUE or FALSE.")
  }
  if ((max.calls%%1 != 0) | (max.calls < 1)){
    stop("max.calls must be a whole number greater than or equal to 1.")
  }
  if (!is.null(lower)){
    if (!(is.vector(lower)) | length(lower) != L + 1){
      stop("lower must be a vector of length L + 1.")
    }
  }
  if (!is.null(upper)){
    if (!(is.vector(upper)) | length(upper) != L + 1){
      stop("upper must be a vector of length L + 1.")
    }
  }
  if (!is.null(init)){
    if (!is.matrix(init) | nrow(init) != num.runs | ncol(init) != L + 1){
      stop("init must be a matrix of dimension num.runs x (L + 1).")
    }
  }
  if (init.file.ext%%1 != 0){
    stop("init.file.ext must be a whole number.")
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

  objFn <- function(parameters, D, L, B, delta, alpha, beta, sigma.e, w, Lambda,
                    penalty, scenarios, opt.crit){
    n              <- parameters[1]
    a              <- parameters[2:(L + 1)]
    I.1            <- numeric(L)
    I.1[L]         <- n*L/(2*sigma.e^2)
    I.1[1:(L - 1)] <- seq_len(L - 1)*I.1[L]/L
    I              <- numeric(2*(D - 1)*L)
    for (l in 1:L){
      I[(1 + (l - 1)*2*(D - 1)):(l*2*(D - 1))] <- I.1[l]
    }
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
                                             mean = (B + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                             sigma = Lambda)[1]
      scenarios[i, 3*(D - 1) + 7] <- pmvnorm(lower = lower,
                                             upper = upper,
                                             mean = (delta + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                             sigma = Lambda)[1]
    }
    fwer   <- sum(scenarios[, 3*(D - 1) + 2]*scenarios[, 3*(D - 1) + 6])
    type.II  <- 1 - sum(scenarios[, 3*(D - 1) + 3]*scenarios[, 3*(D - 1) + 7])
    if (opt.crit == "N"){
      EN.H0  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                        scenarios[, 3*(D - 1) + 6])
      EN.H1  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                        scenarios[, 3*(D - 1) + 7])
      maxN   <- n*L
      scores <- c(EN.H0, EN.H1, maxN)
    } else {
      EO.H0  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                        scenarios[, 3*(D - 1) + 6])
      EO.H1  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                        scenarios[, 3*(D - 1) + 7])
      maxO   <- n*L*D
      scores <- c(EO.H0, EO.H1, maxO)
    }
    Score <- sum(w*scores) + penalty*(as.numeric(fwer > alpha)*(fwer - alpha)/alpha +
                                        as.numeric(type.II > beta)*
                                        (type.II - beta)/beta)
    return(Score)
  }

  constrainedObjFn <- function(parameters, D, L, B, delta, alpha, beta, w, n,
                               I, Lambda, penalty, scenarios, opt.crit){
    a               <- parameters[1:L]
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
                                             mean = (B + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                             sigma = Lambda)[1]
      scenarios[i, 3*(D - 1) + 7] <- pmvnorm(lower = lower,
                                             upper = upper,
                                             mean = (delta + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                             sigma = Lambda)[1]
    }
    fwer     <- sum(scenarios[, 3*(D - 1) + 2]*scenarios[, 3*(D - 1) + 6])
    type.II  <- 1 - sum(scenarios[, 3*(D - 1) + 3]*scenarios[, 3*(D - 1) + 7])
    if (opt.crit == "N"){
      EN.H0  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                        scenarios[, 3*(D - 1) + 6])
      EN.H1  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                        scenarios[, 3*(D - 1) + 7])
      maxN   <- n*L
      scores <- c(EN.H0, EN.H1, maxN)
    } else {
      EO.H0  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                        scenarios[, 3*(D - 1) + 6])
      EO.H1  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                        scenarios[, 3*(D - 1) + 7])
      maxO   <- n*L*D
      scores <- c(EO.H0, EO.H1, maxO)
    }
    Score <- sum(w*scores) + penalty*(as.numeric(fwer > alpha)*(fwer - alpha)/alpha +
                                        as.numeric(type.II > beta)*
                                        (type.II - beta)/beta)
    return(Score)
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  if (delta == 0){
    n.be <- ceiling(2*(qmvnorm(1 - alpha,
                               sigma = matrix(0.5, D - 1, D - 1) +
                                 diag(0.5, D - 1, D - 1))$quantile +
                         qnorm(1 - beta/2))^2*sigma.e^2/B^2)
  } else {
    n.be <- ceiling(2*(qmvnorm(1 - alpha,
                               sigma = matrix(0.5, D - 1, D - 1) +
                                 diag(0.5, D - 1, D - 1))$quantile +
                         qnorm(1 - beta))^2*sigma.e^2/(delta - B)^2)
  }

  if (opt.crit == "N"){
    penalty <- penalty*n.be
  } else {
    penalty <- penalty*n.be*D
  }

  all.run.information      <- matrix(0, ncol = L + 14, nrow = (1 + 2*(find.constrained == TRUE))*num.runs)
  colnames(all.run.information) <- c("Run", "Type", "Algorithm", "Run Time (mins)",
                                     "Score", "E(N|H_N)", "E(O|H_N)", "FWER",
                                     "E(N|H_A)", "E(O|H_A)", "Power", "max N", "max O", "n",
                                     paste("a_", seq_len(L), sep = ""))

  all.omega     <- permutations(n = L, r = D - 1, repeats.allowed = TRUE)
  all.psi.plus  <- permutations(n = 2, r = D - 1, v = 0:1,
                                repeats.allowed = TRUE)
  all.psi.minus  <- permutations(n = 2, r = D - 1, v = 0:1,
                                 repeats.allowed = TRUE)

  all.scenarios <- matrix(0, nrow = L^(D - 1)*2^(D - 1)*2^(D - 1),
                          ncol = 3*(D - 1))
  for (i in 1:(L^(D - 1))){
    for (j in 1:(2^(D - 1))){
      all.scenarios[(1 + (j - 1)*2^(D - 1) + (i - 1)*2^(D - 1)*2^(D - 1)):(j*2^(D - 1) + (i - 1)*2^(D - 1)*2^(D - 1)), ] <- cbind(matrix(all.omega[i, ], nrow = 2^(D - 1),
                                                                                                                                         ncol = D - 1, byrow = TRUE),
                                                                                                                                  matrix(all.psi.plus[j, ], nrow = 2^(D - 1),
                                                                                                                                         ncol = D - 1, byrow = TRUE),
                                                                                                                                  all.psi.minus)
    }
  }
  scenarios.to.drop <- NULL
  for (i in 1:nrow(all.scenarios)){
    for (j in 1:(D - 1)){
      if (all.scenarios[i, j] < L & all.scenarios[i, D - 1 + j] == 1 & all.scenarios[i, 2*(D - 1) + j] == 1){
        scenarios.to.drop <- c(scenarios.to.drop, i)
        break
      }
    }
  }
  all.scenarios <- all.scenarios[-scenarios.to.drop, ]

  degeneracy <- numeric(nrow(all.scenarios))
  degenerate <- NULL
  for (i in 1:(nrow(all.scenarios) - 1)){
    zero.zero.i <- NULL
    zero.one.i  <- NULL
    one.zero.i  <- NULL
    one.one.i   <- NULL
    for (j in 1:(D - 1)){
      if (all.scenarios[i, D - 1 + j] == 0 & all.scenarios[i, 2*(D - 1) + j] == 0){
        zero.zero.i <- c(zero.zero.i, j)
      } else if (all.scenarios[i, D - 1 + j] == 0 & all.scenarios[i, 2*(D - 1) + j] == 1){
        zero.one.i <- c(zero.one.i, j)
      } else if (all.scenarios[i, D - 1 + j] == 1 & all.scenarios[i, 2*(D - 1) + j] == 0){
        one.zero.i <- c(one.zero.i, j)
      } else if (all.scenarios[i, D - 1 + j] == 1 & all.scenarios[i, 2*(D - 1) + j] == 1){
        one.one.i <- c(one.one.i, j)
      }
    }
    for (j in (i + 1):nrow(all.scenarios)){
      if (sum(c(sort(all.scenarios[i, 1:(D - 1)]),
                sort(all.scenarios[i, D:(2*(D - 1))]),
                sort(all.scenarios[i, (2*(D - 1) + 1):(3*(D - 1))])) ==
              c(sort(all.scenarios[j, 1:(D - 1)]),
                sort(all.scenarios[j, D:(2*(D - 1))]),
                sort(all.scenarios[j, (2*(D - 1) + 1):(3*(D - 1))]))) == 3*(D - 1)){
        zero.zero.j <- NULL
        zero.one.j  <- NULL
        one.zero.j  <- NULL
        one.one.j   <- NULL
        for (k in 1:(D - 1)){
          if (all.scenarios[j, D - 1 + k] == 0 & all.scenarios[j, 2*(D - 1) + k] == 0){
            zero.zero.j <- c(zero.zero.j, k)
          } else if (all.scenarios[j, D - 1 + k] == 0 & all.scenarios[j, 2*(D - 1) + k] == 1){
            zero.one.j <- c(zero.one.j, k)
          } else if (all.scenarios[j, D - 1 + k] == 1 & all.scenarios[j, 2*(D - 1) + k] == 0){
            one.zero.j <- c(one.zero.j, k)
          } else if (all.scenarios[j, D - 1 + k] == 1 & all.scenarios[j, 2*(D - 1) + k] == 1){
            one.one.j <- c(one.one.j, k)
          }
        }
        check <- 1
        if (length(zero.zero.i) != length(zero.zero.j) |
            length(zero.one.i) != length(zero.one.j) |
            length(one.zero.i) != length(one.zero.j) |
            length(one.one.i) != length(one.one.j)){
          check <- 0
        }
        if (sum(sort(all.scenarios[i, zero.zero.i]) == sort(all.scenarios[j, zero.zero.j])) != length(zero.zero.i) |
            sum(sort(all.scenarios[i, zero.one.i]) == sort(all.scenarios[j, zero.one.j])) != length(zero.one.i) |
            sum(sort(all.scenarios[i, one.zero.i]) == sort(all.scenarios[j, one.zero.j])) != length(one.zero.i) |
            sum(sort(all.scenarios[i, one.one.i]) == sort(all.scenarios[j, one.one.j])) != length(one.one.i)){
          check <- 0
        }
        if (check == 1){
          degeneracy[i] <- degeneracy[i] + 1
          degenerate    <- c(degenerate, j)
        }
      }
    }
  }
  degeneracy     <- degeneracy + 1
  if (!is.null(degenerate)){
    degeneracy   <- degeneracy[-unique(degenerate)]
    scenarios    <- cbind(all.scenarios[-unique(degenerate), ], degeneracy)
  } else {
    scenarios    <- cbind(all.scenarios, rep(1, nrow(all.scenarios)))
  }
  if (D == 2){
    degeneracy.1 <- numeric(nrow(scenarios))
    for (i in 1:nrow(scenarios)){
      if (scenarios[i, D] == 1 & scenarios[i, 2*(D - 1) + 1] == 1){
        degeneracy.1[i] <- 1
      }
    }
    degeneracy.any <- degeneracy.1
  } else {
    degeneracy.any <- numeric(nrow(scenarios))
    degeneracy.1   <- numeric(nrow(scenarios))
    for (i in 1:nrow(scenarios)){
      scenarios[i, 1:(3*(D - 1))] <- c(rev(scenarios[i, 1:(D - 1)]),
                                       rev(scenarios[i, D:(2*(D - 1))]),
                                       rev(scenarios[i, (2*(D - 1) + 1):(3*(D - 1))]))
    }
    for (i in 1:nrow(scenarios)){
      check <- 0
      for (j in 1:(D - 1)){
        if (scenarios[i, D - 1 + j] == 1 & scenarios[i, 2*(D - 1) + j] == 1){
          check <- 1
          break
        }
      }
      if (check == 1){
        degeneracy.any[i]    <- scenarios[i, 3*(D - 1) + 1]

        poss.ind.omega.perms <- NULL
        poss.ind.psi.plus.perms   <- NULL
        poss.ind.psi.minus.perms   <- NULL
        drop.eff             <- NULL
        for (j in 1:(D - 1)){
          if (scenarios[i, D - 1 + j] == 1 & scenarios[i, 2*(D - 1) + j] == 1){
            drop.eff <- c(drop.eff, j)
          }
        }
        all.poss.scenarios <- NULL
        for (j in drop.eff){
          poss.ind.omega.perms <- as.matrix(unique(permutations(n = D - 2, r = D - 2,
                                                                v = scenarios[i, -c(j, D:(3*(D - 1)))],
                                                                set = FALSE,
                                                                repeats.allowed = FALSE)))
          poss.ind.psi.plus.perms   <- as.matrix(unique(permutations(n = D - 2, r = D - 2,
                                                                     v = scenarios[i, -c(1:(D - 1), D - 1 + j, (2*(D - 1) + 1):(3*(D - 1)))],
                                                                     set = FALSE,
                                                                     repeats.allowed = FALSE)))
          poss.ind.psi.minus.perms   <- as.matrix(unique(permutations(n = D - 2, r = D - 2,
                                                                      v = scenarios[i, -c(1:(2*(D - 1)), 2*(D - 1) + j)],
                                                                      set = FALSE,
                                                                      repeats.allowed = FALSE)))
          poss.scenarios <- matrix(0, nrow = nrow(poss.ind.omega.perms)*nrow(poss.ind.psi.plus.perms)*
                                     nrow(poss.ind.psi.minus.perms),
                                   ncol = 3*(D - 2))
          for (k in 1:nrow(poss.ind.omega.perms)){
            for (l in 1:nrow(poss.ind.psi.plus.perms)){
              poss.scenarios[(1 + (l - 1)*nrow(poss.ind.psi.minus.perms) + (k - 1)*nrow(poss.ind.psi.plus.perms)*nrow(poss.ind.psi.minus.perms)):
                               (l*nrow(poss.ind.psi.minus.perms) + (k - 1)*nrow(poss.ind.psi.plus.perms)*nrow(poss.ind.psi.minus.perms)), ] <- cbind(matrix(poss.ind.omega.perms[k, ],
                                                                                                                                                            nrow = nrow(poss.ind.psi.minus.perms),
                                                                                                                                                            ncol = D - 2, byrow = T), matrix(poss.ind.psi.plus.perms[l, ],
                                                                                                                                                                                             nrow = nrow(poss.ind.psi.minus.perms),
                                                                                                                                                                                             ncol = D - 2, byrow = T),
                                                                                                                                                     poss.ind.psi.minus.perms)
            }
          }
          if (nrow(poss.scenarios) > 1){
            poss.scenarios <- cbind(rep(scenarios[i, j], nrow(poss.scenarios)), poss.scenarios[, 1:(D - 2)],
                                    rep(scenarios[i, D - 1 + j], nrow(poss.scenarios)), poss.scenarios[, (D - 1):(2*(D - 2))],
                                    rep(scenarios[i, 2*(D - 1) + j], nrow(poss.scenarios)), poss.scenarios[, (2*(D - 2) + 1):(3*(D - 2))])
          } else {
            poss.scenarios <- c(rep(scenarios[i, j], nrow(poss.scenarios)), poss.scenarios[, 1:(D - 2)],
                                rep(scenarios[i, D - 1 + j], nrow(poss.scenarios)), poss.scenarios[, (D - 1):(2*(D - 2))],
                                rep(scenarios[i, 2*(D - 1) + j], nrow(poss.scenarios)), poss.scenarios[, (2*(D - 2) + 1):(3*(D - 2))])
          }
          all.poss.scenarios <- rbind(all.poss.scenarios, poss.scenarios)
        }
        if (nrow(all.poss.scenarios) > 1){
          degenerate <- NULL
          for (k in 1:(nrow(all.poss.scenarios) - 1)){
            for (j in (k + 1):nrow(all.poss.scenarios)){
              if (sum(all.poss.scenarios[k, ] ==
                      all.poss.scenarios[j, ]) == 3*(D - 1)){
                degenerate    <- c(degenerate, j)
              }
            }
          }
          if (!is.null(degenerate)){
            all.poss.scenarios  <- all.poss.scenarios[-unique(degenerate), ]
          }
        }
        if (is.vector(all.poss.scenarios)){
          if (all.poss.scenarios[D] == 1 & all.poss.scenarios[2*(D - 1) + 1] == 1){
            degeneracy.1[i] <- degeneracy.1[i] + 1
          }
        } else {
          zero.zero.i <- NULL
          zero.one.i  <- NULL
          one.zero.i  <- NULL
          one.one.i   <- NULL
          for (j in 1:(D - 1)){
            if (scenarios[i, D - 1 + j] == 0 & scenarios[i, 2*(D - 1) + j] == 0){
              zero.zero.i <- c(zero.zero.i, j)
            } else if (scenarios[i, D - 1 + j] == 0 & scenarios[i, 2*(D - 1) + j] == 1){
              zero.one.i <- c(zero.one.i, j)
            } else if (scenarios[i, D - 1 + j] == 1 & scenarios[i, 2*(D - 1) + j] == 0){
              one.zero.i <- c(one.zero.i, j)
            } else if (scenarios[i, D - 1 + j] == 1 & scenarios[i, 2*(D - 1) + j] == 1){
              one.one.i <- c(one.one.i, j)
            }
          }
          for (j in 1:nrow(all.poss.scenarios)){
            if (all.poss.scenarios[j, D] == 1 & all.poss.scenarios[j, 2*(D - 1) + 1] == 1){
              zero.zero.j <- NULL
              zero.one.j  <- NULL
              one.zero.j  <- NULL
              one.one.j   <- NULL
              for (k in 1:(D - 1)){
                if (all.poss.scenarios[j, D - 1 + k] == 0 & all.poss.scenarios[j, 2*(D - 1) + k] == 0){
                  zero.zero.j <- c(zero.zero.j, k)
                } else if (all.poss.scenarios[j, D - 1 + k] == 0 & all.poss.scenarios[j, 2*(D - 1) + k] == 1){
                  zero.one.j <- c(zero.one.j, k)
                } else if (all.poss.scenarios[j, D - 1 + k] == 1 & all.poss.scenarios[j, 2*(D - 1) + k] == 0){
                  one.zero.j <- c(one.zero.j, k)
                } else if (all.poss.scenarios[j, D - 1 + k] == 1 & all.poss.scenarios[j, 2*(D - 1) + k] == 1){
                  one.one.j <- c(one.one.j, k)
                }
              }
              check <- 1
              if (length(zero.zero.i) != length(zero.zero.j) |
                  length(zero.one.i) != length(zero.one.j) |
                  length(one.zero.i) != length(one.zero.j) |
                  length(one.one.i) != length(one.one.j)){
                check <- 0
              }
              if (check == 1){
                if (sum(sort(scenarios[i, zero.zero.i]) == sort(all.poss.scenarios[j, zero.zero.j])) != length(zero.zero.i) |
                    sum(sort(scenarios[i, zero.one.i]) == sort(all.poss.scenarios[j, zero.one.j])) != length(zero.one.i) |
                    sum(sort(scenarios[i, one.zero.i]) == sort(all.poss.scenarios[j, one.zero.j])) != length(one.zero.i) |
                    sum(sort(scenarios[i, one.one.i]) == sort(all.poss.scenarios[j, one.one.j])) != length(one.one.i)){
                  check <- 0
                }
              }
              if (check == 1){
                degeneracy.1[i] <- degeneracy.1[i] + 1
              }
            }
          }
        }
      }
    }
  }

  scenarios           <- cbind(scenarios, degeneracy.any, degeneracy.1,
                               matrix(0, nrow = nrow(scenarios), ncol = 4))
  colnames(scenarios) <- c(paste("omega_", 1:(D - 1), sep = ""),
                           paste("psi(+)_", 1:(D - 1), sep = ""),
                           paste("psi(-)_", 1:(D - 1), sep = ""),
                           "deg", "deg.any", "deg.1", "N", "O",
                           "P(H0)", "P(H1)")
  for (i in 1:nrow(scenarios)){
    scenarios[i, 3*(D - 1) + 4] <- max(scenarios[i, 1:(D - 1)])
    for (j in 1:L){
      if (length(which(scenarios[i, 1:(D - 1)] >= j)) > 0){
        scenarios[i, 3*(D - 1) + 5] <- scenarios[i, 3*(D - 1) + 5] +
          length(which(scenarios[i, 1:(D - 1)] >= j)) + 1
      }
    }
  }

  if (is.null(lower)){
    lower.b <- c(0, rep(-20, L))
  } else {
    lower.b <- lower
  }
  if (is.null(upper)){
    upper.b <- c(n.be, rep(20, L))
  } else {
    upper.b <- upper
  }
  if (is.null(init)){
    initial.values            <- matrix(0, ncol = L + 1, nrow = num.runs)
    for (k in 1:(L + 1)){
      initial.values[, k]     <- runif(num.runs, min = lower.b[k], max = upper.b[k])
    }
  } else {
    initial.values <- init
  }
  word <- c("lower", "upper")

  Lambda <- covariance(D, L)
  lcm       <- 2
  sink("NULL")
  if (D > 2) {
    for (i in 3:D) {
      S.i   <- i
      if ((sequence.type == "williams") & (is.odd(i))) {
        S.i <- 2*i
      }
      lcm   <- scm(lcm, S.i)
    }
  }
  sink()

  counter <- 0
  for (i in 1:num.runs){

    if (summary == TRUE){
      print(paste("Beginning Run ", i, ". Output from chosen algorithm may follow...", sep = ""))
    }

    main.start.time.i    <- Sys.time()
    if (algorithm == "DEoptim"){
      ff         <- 10*length(lower.b)
      mi         <- ceiling(max.calls/ff)
      DEoptim.i  <- try(DEoptim(fn = objFn, lower = lower.b, upper = upper.b,
                                control = list(itermax = mi), D = D, L = L, B = B,
                                delta = delta,
                                sigma.e = sigma.e, alpha = alpha, beta = beta,
                                Lambda = Lambda, w = w, penalty = penalty,
                                scenarios = scenarios, opt.crit = opt.crit))
      if (class(DEoptim.i) == "try-error"){
        Score <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- DEoptim.i$optim$bestval
        opt.parameters <- DEoptim.i$optim$bestmem
      }
    } else if (algorithm == "DEopt"){
      DEopt.i <- try(DEopt(OF = objFn,
                       algo = list(min=lower.b, max=upper.b,
                                   nG = ceiling(max.calls/50), minmaxConstr = TRUE,
                                   printDetail = summary,
                                   printBar = summary), D = D, L = L,
                       delta = delta, B = B, sigma.e = sigma.e, alpha = alpha,
                       beta = beta, Lambda = Lambda, w = w, penalty = penalty,
                       scenarios = scenarios,
                       opt.crit = opt.crit))
      if (class(DEopt.i) == "try-error"){
        Score <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- DEopt.i$OFvalue
        opt.parameters <- DEopt.i$xbest
      }
    } else if (algorithm == "ga"){
      fitness <- function(...){
        -objFn(...)
      }
      ga.i       <- try(ga(type = "real-valued", fitness = fitness, min = lower.b,
                           max = upper.b, maxiter = ceiling(max.calls*1230/50000), D = D, L = L,  B = B, delta = delta,
                           sigma.e = sigma.e, alpha = alpha, beta = beta,
                           Lambda = Lambda, w = w, penalty = penalty, scenarios = scenarios,
                           opt.crit = opt.crit))
      if (class(ga.i) == "try-error"){
        Score <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score          <- -ga.i@fitnessValue
        opt.parameters <- ga.i@solution
      }
    } else if (algorithm == "GenSA"){
      start.time <- Sys.time()
      GenSA.i    <- try(GenSA(par = initial.values[i, ], fn = objFn, lower = lower.b,
                              upper = upper.b, control = list(max.call = max.calls,
                                                              verbose = summary),
                              D = D, L = L,  B = B, delta = delta, sigma.e = sigma.e, alpha = alpha,
                              beta = beta, Lambda = Lambda, w = w, penalty = penalty,
                              scenarios = scenarios,
                              opt.crit = opt.crit))
      end.time <- Sys.time()
      if (class(GenSA.i) == "try-error"){
        Score <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- GenSA.i$value
        opt.parameters <- GenSA.i$par
      }
    } else if (algorithm == "hydroPSO"){
      hydroPSO.i <- try(hydroPSO(par = initial.values[i, ], fn = objFn, lower = lower.b,
                                 upper = upper.b, control = list(maxfn = max.calls,
                                                                 verbose = summary),
                                 D=D, L = L,  B = B, delta=delta, sigma.e = sigma.e, alpha = alpha,
                                 beta = beta, Lambda = Lambda, w = w,
                                 penalty = penalty, scenarios = scenarios,
                                 opt.crit = opt.crit))
      if (class(hydroPSO.i) == "try-error"){
        Score <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- hydroPSO.i$value
        opt.parameters <- hydroPSO.i$par
      }
    } else if (algorithm == "PSopt"){
      PSopt.i <- try(PSopt(OF = objFn,
                           algo = list(min=lower.b, max=upper.b,
                                       nG = ceiling(max.calls/100),
                                       printDetail = summary,
                                       printBar = summary), D = D, L = L,
                           delta = delta, B = B, sigma.e = sigma.e, alpha = alpha,
                           beta = beta, Lambda = Lambda, w = w, penalty = penalty,
                           scenarios = scenarios,
                           opt.crit = opt.crit))
      if (class(PSopt.i) == "try-error"){
        Score <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- PSopt.i$OFvalue
        opt.parameters <- PSopt.i$xbest
      }
    } else if (algorithm == "NLOPT_GD_STOGO_RAND"){
      NLOPT_GD_STOGO_RAND.i <- try(stogo(x0 = initial.values[i, ], fn = objFn,
                                         lower = lower.b, upper = upper.b,
                                         maxeval = max.calls, D = D, L = L,  B = B, delta = delta,
                                         sigma.e = sigma.e, alpha = alpha,
                                         beta = beta, Lambda = Lambda, w = w,
                                         penalty = penalty, scenarios = scenarios,
                                         opt.crit = opt.crit, nl.info = summary))
      if (class(NLOPT_GD_STOGO_RAND.i) == "try-error"){
        Score <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- NLOPT_GD_STOGO_RAND.i$value
        opt.parameters <- NLOPT_GD_STOGO_RAND.i$par
      }
    } else if (algorithm == "NLOPT_GN_CRS2_LM"){
      NLOPT_GN_CRS2_LM.i <- try(crs2lm(x0 = initial.values[i, ], fn = objFn,
                                       lower = lower.b, upper = upper.b,
                                       maxeval = max.calls, D = D, L = L,  B = B,
                                       delta = delta,
                                       sigma.e = sigma.e, alpha = alpha,
                                       beta = beta, Lambda = Lambda, w = w,
                                       penalty = penalty, scenarios = scenarios,
                                       opt.crit = opt.crit, nl.info = summary))
      if (class(NLOPT_GN_CRS2_LM.i) == "try-error"){
        Score          <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- NLOPT_GN_CRS2_LM.i$value
        opt.parameters  <- NLOPT_GN_CRS2_LM.i$par
      }
    } else if (algorithm == "NLOPT_GN_DIRECT"){
      NLOPT_GN_DIRECT.i <- try(nloptr(x0 = initial.values[i, ], eval_f = objFn,
                                      lb = lower.b, ub = upper.b,
                                      opts = list(algorithm = "NLOPT_GN_DIRECT",
                                                  maxeval = max.calls,
                                                  print_level = as.numeric(summary)),
                                      D = D, L = L,  B = B, delta = delta, sigma.e = sigma.e,
                                      alpha = alpha, beta = beta, Lambda = Lambda,
                                      w = w, penalty = penalty, scenarios = scenarios,
                                      opt.crit = opt.crit))
      if (class(NLOPT_GN_DIRECT.i) == "try-error"){
        Score          <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- NLOPT_GN_DIRECT.i$objective
        opt.parameters <- NLOPT_GN_DIRECT.i$solution
      }
    } else if (algorithm == "NLOPT_GN_DIRECT_L"){
      NLOPT_GN_DIRECT_L.i <- try(nloptr(x0 = initial.values[i, ], eval_f = objFn,
                                        lb = lower.b, ub = upper.b,
                                        opts = list(algorithm = "NLOPT_GN_DIRECT_L",
                                                    maxeval = max.calls,
                                                    print_level = as.numeric(summary)),
                                        D = D, L = L,  B = B,
                                        delta = delta, sigma.e = sigma.e, alpha = alpha,
                                        beta = beta, Lambda = Lambda, w = w,
                                        penalty = penalty, scenarios = scenarios,
                                        opt.crit = opt.crit))
      if (class(NLOPT_GN_DIRECT_L.i) == "try-error"){
        Score          <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- NLOPT_GN_DIRECT_L.i$objective
        opt.parameters <- NLOPT_GN_DIRECT_L.i$solution
      }
    } else if (algorithm == "NLOPT_GN_ISRES"){
      NLOPT_GN_ISRES.i <- try(nloptr(x0 = initial.values[i, ], eval_f = objFn,
                                     lb = lower.b, ub = upper.b,
                                     opts = list(algorithm = "NLOPT_GN_ISRES",
                                                 maxeval = max.calls,
                                                 print_level = as.numeric(summary)), D = D, L = L,  B = B,
                                     delta = delta, sigma.e = sigma.e, alpha = alpha,
                                     beta = beta, Lambda = Lambda, w = w,
                                     penalty = penalty, scenarios = scenarios,
                                     opt.crit = opt.crit))
      if (class(NLOPT_GN_ISRES.i) == "try-error"){
        Score          <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- NLOPT_GN_ISRES.i$objective
        opt.parameters <- NLOPT_GN_ISRES.i$solution
      }
    } else if (algorithm == "psoptim"){
      psoptim.i  <- try(psoptim(par = initial.values[i, ], fn = objFn, lower = lower.b,
                                upper = upper.b, control = list(maxit = ceiling(max.calls*4170/50000),
                                                                trace.stats = summary),
                                D = D, L = L, B = B,
                                delta = delta, sigma.e = sigma.e, alpha = alpha,
                                beta = beta, Lambda = Lambda, w = w, penalty = penalty,
                                scenarios = scenarios,
                                opt.crit = opt.crit))
      if (class(psoptim.i) == "try-error"){
        Score          <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- psoptim.i$value
        opt.parameters  <- psoptim.i$par
      }
    } else if (algorithm == "SANN"){
      SANN.i     <- try(optim(par = initial.values[i, ], fn = objFn, method = "SANN",
                              control = list(maxit = max.calls), D = D, L = L, B = B,
                              delta = delta,
                              sigma.e = sigma.e, alpha = alpha, beta = beta,
                              Lambda = Lambda, w = w, penalty = penalty, scenarios = scenarios,
                              opt.crit = opt.crit))
      if (class(SANN.i) == "try-error"){
        Score          <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- SANN.i$value
        opt.parameters <- SANN.i$par
      }
    } else if (algorithm == "soma"){
      soma.i     <- try(soma(costFunction = objFn,
                             bounds = list(min = lower.b, max = upper.b),
                             options = list(nMigrations = ceiling(max.calls*36/50000)), D = D, L = L,  B = B,
                             delta = delta,
                             sigma.e = sigma.e, alpha = alpha, beta = beta,
                             Lambda = Lambda, w = w, penalty = penalty, scenarios = scenarios,
                             opt.crit = opt.crit))
      if (class(soma.i) == "try-error"){
        Score          <- NA
        opt.parameters <- rep(NA, length(lower.b))
      } else {
        Score           <- soma.i$cost[soma.i$leader]
        opt.parameters <- soma.i$population[, soma.i$leader]
      }
    }
    if (!is.na(Score)){
      n              <- opt.parameters[1]
      I.1            <- numeric(L)
      I.1[L]         <- n*L/(2*sigma.e^2)
      I.1[1:(L - 1)] <- seq_len(L - 1)*I.1[L]/L
      I              <- numeric(2*(D - 1)*L)
      for (l in 1:L){
        I[(1 + (l - 1)*2*(D - 1)):(l*2*(D - 1))] <- I.1[l]
      }
      a              <- opt.parameters[2:(L + 1)]
      n.exact   <- n
      I.exact   <- I
      I.1.exact <- I.1
      a.exact   <- a
      ind.lower.minus <- numeric(L)
      ind.upper.minus <- numeric(L)
      ind.lower.plus  <- numeric(L)
      ind.upper.plus  <- numeric(L)
      lower           <- numeric(2*(D - 1)*L)
      upper           <- numeric(2*(D - 1)*L)
      for (k in 1:nrow(scenarios)){
        for (j in 1:(D - 1)){
          if (scenarios[k, j] == 1){
            if (scenarios[k, D - 1 + j] == 0){
              ind.lower.plus <- rep(-Inf, L)
              ind.upper.plus <- c(a[1], rep(Inf, L - 1))
            } else {
              ind.lower.plus <- c(a[1], rep(-Inf, L - 1))
              ind.upper.plus <- rep(Inf, L)
            }
            if (scenarios[k, 2*(D - 1) + j] == 0){
              ind.lower.minus <- c(-a[1], rep(-Inf, L - 1))
              ind.upper.minus <- rep(Inf, L)
            } else {
              ind.lower.minus <- rep(-Inf, L)
              ind.upper.minus <- c(-a[1], rep(Inf, L - 1))
            }
          } else if (scenarios[k, j] > 1 & scenarios[k, j] < L){
            if (scenarios[k, D - 1 + j] == 0){
              ind.lower.plus <- c(a[1:(scenarios[k, j] - 1)], rep(-Inf, L - (scenarios[k, j] - 1)))
              ind.upper.plus <- c(rep(Inf, scenarios[k, j] - 1), a[scenarios[k, j]],
                                  rep(Inf, L - scenarios[k, j]))
            } else {
              ind.lower.plus <- c(a[1:scenarios[k, j]], rep(-Inf, L - scenarios[k, j]))
              ind.upper.plus <- rep(Inf, L)
            }
            if (scenarios[k, 2*(D - 1) + j] == 0){
              ind.lower.minus <- c(rep(-Inf, scenarios[k, j] - 1), -a[scenarios[k, j]],
                                   rep(-Inf, L - scenarios[k, j]))
              ind.upper.minus <- c(-a[1:(scenarios[k, j] - 1)], rep(Inf, L - (scenarios[k, j] - 1)))
            } else {
              ind.lower.minus <- rep(-Inf, L)
              ind.upper.minus <- c(-a[1:scenarios[k, j]], rep(Inf, L - scenarios[k, j]))
            }
          } else {
            if (scenarios[k, D - 1 + j] == 0){
              ind.lower.plus <- c(a[1:(L - 1)], -Inf)
              ind.upper.plus <- c(rep(Inf, L - 1), a[L])
            } else {
              ind.lower.plus <- a
              ind.upper.plus <- rep(Inf, L)
            }
            if (scenarios[k, 2*(D - 1) + j] == 0){
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
        scenarios[k, 3*(D - 1) + 6] <- pmvnorm(lower = lower,
                                               upper = upper,
                                               mean = (B + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                               sigma = Lambda)[1]
        scenarios[k, 3*(D - 1) + 7] <- pmvnorm(lower = lower,
                                               upper = upper,
                                               mean = (delta + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                               sigma = Lambda)[1]
      }
      fwer   <- sum(scenarios[, 3*(D - 1) + 2]*scenarios[, 3*(D - 1) + 6])
      power  <- sum(scenarios[, 3*(D - 1) + 3]*scenarios[, 3*(D - 1) + 7])
      EN.H0  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                          scenarios[, 3*(D - 1) + 6])
      EN.H1  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                          scenarios[, 3*(D - 1) + 7])
      maxN   <- n*L
      EO.H0  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                        scenarios[, 3*(D - 1) + 6])
      EO.H1  <- n*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                       scenarios[, 3*(D - 1) + 7])
      maxO   <- n*L*D
      perf.H0 <- c(EN.H0, EO.H0, fwer)
      perf.H1 <- c(EN.H1, EO.H1, power)
      scenarios.exact <- scenarios
      scenarios.exact[, c(3*(D - 1) + 4, 3*(D - 1) + 5)] <- n*scenarios.exact[, c(3*(D - 1) + 4, 3*(D - 1) + 5)]
    } else {
      perf.H0 <- rep(NA, 3)
      perf.H1 <- rep(NA, 3)
      maxN    <- NA
      maxO    <- NA
      n       <- NA
      a       <- rep(NA, L)
      scenarios.exact <- NA
      n.exact   <- NA
      I.exact   <- NA
      I.1.exact <- NA
      a.exact <- NA
    }
    main.end.time.i    <- Sys.time()
    all.run.information[1 + (i - 1)*(1 + 2*(find.constrained == TRUE)), ] <- c(i, "Unconstrained", algorithm, as.numeric(difftime(main.end.time.i,
                                                                               main.start.time.i,
                                                                               units = "mins")),
                                              Score, perf.H0, perf.H1, maxN, maxO, n, a)
    write.csv(all.run.information, paste(algorithm, init.file.ext - 1 + i, ".csv", sep = ""))
    if (find.constrained == TRUE){
      scores <- NULL
      if (!is.na(Score)){
        n      <- opt.parameters[1]
        poss.n <- c(floor(n), ceiling(n))
        if (poss.n[1] != poss.n[2]){
          while (poss.n[1]%%lcm != 0){
            poss.n[1] <- poss.n[1] - 1
          }
          while (poss.n[2]%%lcm != 0){
            poss.n[2] <- poss.n[2] + 1
          }
        }
        print(poss.n)
        if (poss.n[1] == 0){
          poss.n <- poss.n[2]
        }
        initial.values.constrained <- opt.parameters[-1]
        lower.constrained <- lower.b[-1]
        upper.constrained <- upper.b[-1]
        possible.opt.param <- matrix(0, nrow = 2, ncol = length(initial.values.constrained))
        scores             <- numeric(2)
        times              <- numeric(2)
        for (j in 1:length(poss.n)){
          start.time.j   <- Sys.time()
          I.1            <- numeric(L)
          I.1[L]         <- poss.n[j]*L/(2*sigma.e^2)
          I.1[1:(L - 1)] <- seq_len(L - 1)*I.1[L]/L
          I              <- numeric(2*(D - 1)*L)
          for (l in 1:L){
            I[(1 + (l - 1)*2*(D - 1)):(l*2*(D - 1))] <- I.1[l]
          }
          if (algorithm == "DEopt"){
            DEopt.i <- try(DEopt(OF = objFn,
                                 algo = list(min=lower.b, max=upper.b,
                                             nG = ceiling(max.calls/50), minmaxConstr = TRUE,
                                             printDetail = summary,
                                             printBar = summary), D = D, L = L, B = B,
                                 delta = delta, n = poss.n[j], I = I,
                                 alpha = alpha, beta = beta,
                                 Lambda = Lambda, w = w, penalty = penalty,
                                 scenarios = scenarios, opt.crit = opt.crit))
            if (class(DEopt.i) == "try-error"){
              Score <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- DEopt.i$OFvalue
              opt.parameters <- DEopt.i$xbest
            }
          } else if (algorithm == "DEoptim"){
            ff         <- 10*length(lower)
            mi         <- round(max.calls/ff)
            DEoptim.i  <- try(DEoptim(fn = constrainedObjFn, lower = lower.constrained, upper = upper.constrained,
                                      control = list(itermax = mi), D = D, L = L, B = B,
                                      delta = delta, n = poss.n[j], I = I,
                                      alpha = alpha, beta = beta,
                                      Lambda = Lambda, w = w, penalty = penalty,
                                      scenarios = scenarios, opt.crit = opt.crit))
            if (class(DEoptim.i) == "try-error"){
              Score <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- DEoptim.i$optim$bestval
              opt.parameters <- DEoptim.i$optim$bestmem
            }
          } else if (algorithm == "ga"){
            fitness <- function(...){
              -constrainedObjFn(...)
            }
            ga.i       <- try(ga(type = "real-valued", fitness = fitness, min = lower.constrained,
                                 max = upper.constrained, maxiter = ceiling(max.calls*1230/50000), D = D, L = L,  B = B, delta = delta,
                                 alpha = alpha, beta = beta, n = poss.n[j], I = I,
                                 Lambda = Lambda, w = w, penalty = penalty, scenarios = scenarios,
                                 opt.crit = opt.crit))
            if (class(ga.i) == "try-error"){
              Score <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score          <- -ga.i@fitnessValue
              opt.parameters <- ga.i@solution
            }
          } else if (algorithm == "GenSA"){
            GenSA.i    <- try(GenSA(par = initial.values.constrained, fn = constrainedObjFn,
                                    lower = lower.constrained,
                                    upper = upper.constrained, control = list(max.call = max.calls,
                                                                              verbose = summary),
                                    D = D, L = L,  B = B, delta = delta, alpha = alpha, n = poss.n[j], I = I,
                                    beta = beta, Lambda = Lambda, w = w, penalty = penalty,
                                    scenarios = scenarios,
                                    opt.crit = opt.crit))
            if (class(GenSA.i) == "try-error"){
              Score <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- GenSA.i$value
              opt.parameters <- GenSA.i$par
            }
          } else if (algorithm == "hydroPSO"){
            hydroPSO.i <- try(hydroPSO(par = initial.values.constrained, fn = constrainedObjFn,
                                       lower = lower.constrained,
                                       upper = upper.constrained, control = list(maxfn = max.calls,
                                                                                 verbose = summary),
                                       D=D, L = L,  B = B, delta=delta, alpha = alpha, n = poss.n[j], I = I,
                                       beta = beta, Lambda = Lambda, w = w,
                                       penalty = penalty, scenarios = scenarios,
                                       opt.crit = opt.crit))
            if (class(hydroPSO.i) == "try-error"){
              Score <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- hydroPSO.i$value
              opt.parameters <- hydroPSO.i$par
            }
          } else if (algorithm == "NLOPT_GD_STOGO_RAND"){
            NLOPT_GD_STOGO_RAND.i <- try(stogo(x0 = initial.values.constrained,
                                               fn = constrainedObjFn,
                                               lower = lower.constrained, upper = upper.constrained,
                                               maxeval = max.calls, D = D, L = L,  B = B, delta = delta,
                                               alpha = alpha, n = poss.n[j], I = I,
                                               beta = beta, Lambda = Lambda, w = w,
                                               penalty = penalty, scenarios = scenarios,
                                               opt.crit = opt.crit, nl.info = summary))
            if (class(NLOPT_GD_STOGO_RAND.i) == "try-error"){
              Score <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- NLOPT_GD_STOGO_RAND.i$value
              opt.parameters <- NLOPT_GD_STOGO_RAND.i$par
            }
          } else if (algorithm == "NLOPT_GN_CRS2_LM"){
            NLOPT_GN_CRS2_LM.i <- try(crs2lm(x0 = initial.values.constrained,
                                             fn = constrainedObjFn,
                                             lower = lower.constrained, upper = upper.constrained,
                                             maxeval = max.calls, D = D, L = L,  B = B,
                                             delta = delta, n = poss.n[j], I = I,
                                             alpha = alpha,
                                             beta = beta, Lambda = Lambda, w = w,
                                             penalty = penalty, scenarios = scenarios,
                                             opt.crit = opt.crit, nl.info = summary))
            if (class(NLOPT_GN_CRS2_LM.i) == "try-error"){
              Score          <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- NLOPT_GN_CRS2_LM.i$value
              opt.parameters  <- NLOPT_GN_CRS2_LM.i$par
            }
          } else if (algorithm == "NLOPT_GN_DIRECT"){
            NLOPT_GN_DIRECT.i <- try(nloptr(x0 = initial.values.constrained,
                                            eval_f = constrainedObjFn,
                                            lb = lower.constrained, ub = upper.constrained,
                                            opts = list(algorithm = "NLOPT_GN_DIRECT",
                                                        maxeval = max.calls,
                                                        print_level = as.numeric(summary)),
                                            D = D, L = L,  B = B, delta = delta, n = poss.n[j], I = I,
                                            alpha = alpha, beta = beta, Lambda = Lambda,
                                            w = w, penalty = penalty, scenarios = scenarios,
                                            opt.crit = opt.crit))
            if (class(NLOPT_GN_DIRECT.i) == "try-error"){
              Score          <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- NLOPT_GN_DIRECT.i$objective
              opt.parameters <- NLOPT_GN_DIRECT.i$solution
            }
          } else if (algorithm == "NLOPT_GN_DIRECT_L"){
            NLOPT_GN_DIRECT_L.i <- try(nloptr(x0 = initial.values.constrained,
                                              eval_f = constrainedObjFn,
                                              lb = lower.constrained, ub = upper.constrained,
                                              opts = list(algorithm = "NLOPT_GN_DIRECT_L",
                                                          maxeval = max.calls,
                                                          print_level = as.numeric(summary)), D = D, L = L,  B = B,
                                              delta = delta, alpha = alpha, n = poss.n[j], I = I,
                                              beta = beta, Lambda = Lambda, w = w,
                                              penalty = penalty, scenarios = scenarios,
                                              opt.crit = opt.crit))
            if (class(NLOPT_GN_DIRECT_L.i) == "try-error"){
              Score          <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- NLOPT_GN_DIRECT_L.i$objective
              opt.parameters <- NLOPT_GN_DIRECT_L.i$solution
            }
          } else if (algorithm == "NLOPT_GN_ISRES"){
            NLOPT_GN_ISRES.i <- try(nloptr(x0 = initial.values.constrained,
                                           eval_f = constrainedObjFn,
                                           lb = lower.constrained, ub = upper.constrained,
                                           opts = list(algorithm = "NLOPT_GN_ISRES",
                                                       maxeval = max.calls,
                                                       print_level = as.numeric(summary)), D = D, L = L,  B = B,
                                           delta = delta, alpha = alpha, n = poss.n[j], I = I,
                                           beta = beta, Lambda = Lambda, w = w,
                                           penalty = penalty, scenarios = scenarios,
                                           opt.crit = opt.crit))
            if (class(NLOPT_GN_ISRES.i) == "try-error"){
              Score          <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- NLOPT_GN_ISRES.i$objective
              opt.parameters <- NLOPT_GN_ISRES.i$solution
            }
          } else if (algorithm == "PSopt"){
            PSopt.i <- try(PSopt(OF = objFn,
                                 algo = list(min=lower.b, max=upper.b,
                                             nG = ceiling(max.calls/100),
                                             printDetail = summary,
                                             printBar = summary), D = D, L = L, B = B,
                                 delta = delta, alpha = alpha, n = poss.n[j], I = I,
                                 beta = beta, Lambda = Lambda, w = w, penalty = penalty,
                                 scenarios = scenarios,
                                 opt.crit = opt.crit))
            if (class(PSopt.i) == "try-error"){
              Score <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score          <- PSopt.i$OFvalue
              opt.parameters <- PSopt.i$xbest
            }
          } else if (algorithm == "psoptim"){
            psoptim.i  <- try(psoptim(par = initial.values.constrained,
                                      fn = constrainedObjFn,
                                      lower = lower.constrained,
                                      upper = upper.constrained, control = list(maxit = ceiling(max.calls*4170/50000),
                                                                                trace.stats = summary),
                                      D = D, L = L, B = B,
                                      delta = delta, alpha = alpha, n = poss.n[j], I = I,
                                      beta = beta, Lambda = Lambda, w = w, penalty = penalty,
                                      scenarios = scenarios,
                                      opt.crit = opt.crit))
            if (class(psoptim.i) == "try-error"){
              Score          <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- psoptim.i$value
              opt.parameters  <- psoptim.i$par
            }
          } else if (algorithm == "SANN"){
            SANN.i     <- try(optim(par = initial.values.constrained,
                                    fn = constrainedObjFn, method = "SANN",
                                    control = list(maxit = max.calls), D = D, L = L, B = B,
                                    delta = delta,
                                    alpha = alpha, beta = beta, n = poss.n[j], I = I,
                                    Lambda = Lambda, w = w, penalty = penalty, scenarios = scenarios,
                                    opt.crit = opt.crit))
            if (class(SANN.i) == "try-error"){
              Score          <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- SANN.i$value
              opt.parameters <- SANN.i$par
            }
          } else if (algorithm == "soma"){
            soma.i     <- try(soma(costFunction = constrainedObjFn,
                                   bounds = list(min = lower.constrained, max = upper.constrained),
                                   options = list(nMigrations = 36), D = D, L = L,  B = B,
                                   delta = delta, n = poss.n[j], I = I,
                                   alpha = alpha, beta = beta,
                                   Lambda = Lambda, w = w, penalty = penalty, scenarios = scenarios,
                                   opt.crit = opt.crit))
            if (class(soma.i) == "try-error"){
              Score          <- NA
              opt.parameters <- rep(NA, length(lower.b) - 1)
            } else {
              Score           <- soma.i$cost[soma.i$leader]
              opt.parameters <- soma.i$population[, soma.i$leader]
            }
          }

          possible.opt.param[j, ] <- opt.parameters
          scores[j]               <- Score
          if (!is.na(Score)){
            a <- opt.parameters
            ind.lower.minus <- numeric(L)
            ind.upper.minus <- numeric(L)
            ind.lower.plus  <- numeric(L)
            ind.upper.plus  <- numeric(L)
            lower           <- numeric(2*(D - 1)*L)
            upper           <- numeric(2*(D - 1)*L)
            for (k in 1:nrow(scenarios)){
              for (l in 1:(D - 1)){
                if (scenarios[k, l] == 1){
                  if (scenarios[k, D - 1 + l] == 0){
                    ind.lower.plus <- rep(-Inf, L)
                    ind.upper.plus <- c(a[1], rep(Inf, L - 1))
                  } else {
                    ind.lower.plus <- c(a[1], rep(-Inf, L - 1))
                    ind.upper.plus <- rep(Inf, L)
                  }
                  if (scenarios[k, 2*(D - 1) + l] == 0){
                    ind.lower.minus <- c(-a[1], rep(-Inf, L - 1))
                    ind.upper.minus <- rep(Inf, L)
                  } else {
                    ind.lower.minus <- rep(-Inf, L)
                    ind.upper.minus <- c(-a[1], rep(Inf, L - 1))
                  }
                } else if (scenarios[k, l] > 1 & scenarios[k, l] < L){
                  if (scenarios[k, D - 1 + l] == 0){
                    ind.lower.plus <- c(a[1:(scenarios[k, l] - 1)], rep(-Inf, L - (scenarios[k, l] - 1)))
                    ind.upper.plus <- c(rep(Inf, scenarios[k, l] - 1), a[scenarios[k, l]],
                                        rep(Inf, L - scenarios[k, l]))
                  } else {
                    ind.lower.plus <- c(a[1:scenarios[k, l]], rep(-Inf, L - scenarios[k, l]))
                    ind.upper.plus <- rep(Inf, L)
                  }
                  if (scenarios[k, 2*(D - 1) + l] == 0){
                    ind.lower.minus <- c(rep(-Inf, scenarios[k, l] - 1), -a[scenarios[k, l]],
                                         rep(-Inf, L - scenarios[k, l]))
                    ind.upper.minus <- c(-a[1:(scenarios[k, l] - 1)], rep(Inf, L - (scenarios[k, l] - 1)))
                  } else {
                    ind.lower.minus <- rep(-Inf, L)
                    ind.upper.minus <- c(-a[1:scenarios[k, l]], rep(Inf, L - scenarios[k, l]))
                  }
                } else {
                  if (scenarios[k, D - 1 + l] == 0){
                    ind.lower.plus <- c(a[1:(L - 1)], -Inf)
                    ind.upper.plus <- c(rep(Inf, L - 1), a[L])
                  } else {
                    ind.lower.plus <- a
                    ind.upper.plus <- rep(Inf, L)
                  }
                  if (scenarios[k, 2*(D - 1) + l] == 0){
                    ind.lower.minus <- c(rep(-Inf, L - 1), -a[L])
                    ind.upper.minus <- c(-a[1:(L - 1)], Inf)
                  } else {
                    ind.lower.minus <- rep(-Inf, L)
                    ind.upper.minus <- -a
                  }
                }
                lower[seq(from = 1 + 2*(l - 1), by = 2*(D - 1),
                          length.out = L)] <- ind.lower.minus
                lower[seq(from = 2 + 2*(l - 1), by = 2*(D - 1),
                          length.out = L)] <- ind.lower.plus
                upper[seq(from = 1 + 2*(l - 1), by = 2*(D - 1),
                          length.out = L)] <- ind.upper.minus
                upper[seq(from = 2 + 2*(l - 1), by = 2*(D - 1),
                          length.out = L)] <- ind.upper.plus
              }
              scenarios[k, 3*(D - 1) + 6] <- pmvnorm(lower = lower,
                                                     upper = upper,
                                                     mean = (B + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                                     sigma = Lambda)[1]
              scenarios[k, 3*(D - 1) + 7] <- pmvnorm(lower = lower,
                                                     upper = upper,
                                                     mean = (delta + rep(c(-B, B), (D - 1)*L))*sqrt(I),
                                                     sigma = Lambda)[1]
            }
            fwer   <- sum(scenarios[, 3*(D - 1) + 2]*scenarios[, 3*(D - 1) + 6])
            power    <- sum(scenarios[, 3*(D - 1) + 3]*scenarios[, 3*(D - 1) + 7])
            EN.H0  <- poss.n[j]*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                              scenarios[, 3*(D - 1) + 6])
            EN.H1  <- poss.n[j]*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 4]*
                              scenarios[, 3*(D - 1) + 7])
            maxN   <- poss.n[j]*L
            EO.H0  <- poss.n[j]*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                              scenarios[, 3*(D - 1) + 6])
            EO.H1  <- poss.n[j]*sum(scenarios[, 3*(D - 1) + 1]*scenarios[, 3*(D - 1) + 5]*
                              scenarios[, 3*(D - 1) + 7])
            maxO   <- poss.n[j]*L*D
            perf.H0 <- c(EN.H0, EO.H0, fwer)
            perf.H1 <- c(EN.H1, EO.H1, power)
            end.time.j              <- Sys.time()
            times[j]                <- as.numeric(difftime(end.time.j,
                                                           start.time.j,
                                                           units = "mins"))
            all.run.information[1 + j + (i - 1)*3, ] <- c(i, paste("Constrained ", word[j], sep = ""), algorithm, times[j],
                                                          Score, perf.H0, perf.H1, maxN, maxO, poss.n[j], a)
          } else {
            end.time.j              <- Sys.time()
            times[j]                <- as.numeric(difftime(end.time.j,
                                                           start.time.j,
                                                           units = "mins"))
            all.run.information[1 + j + (i - 1)*3, ] <- NA
            all.run.information[1 + j + (i - 1)*3, 1:3] <- c(paste(i, ", constrained ", word[j], sep = ""), algorithm, times[j])
          }
        }
      } else {
        all.run.information[2 + (i - 1)*3, ] <- NA
        all.run.information[2 + (i - 1)*3, 1:3] <- c(paste(i, ", constrained lower", sep = ""), algorithm, 0)
        all.run.information[3 + (i - 1)*3, ] <- NA
        all.run.information[3 + (i - 1)*3, 1:3] <- c(paste(i, ", constrained upper", sep = ""), algorithm, 0)
      }
      write.csv(all.run.information, paste(algorithm, init.file.ext - 1 + i, " - with constrained", ".csv", sep = ""))
      if (is.null(scores)){
        scenarios.optimal <- NA
        I.optimal   <- NA
        I.1.optimal <- NA
        n.optimal  <- NA
        a.optimal <- NA
      } else {
        if (all(is.numeric(scores))){
          optimal <- which.min(scores)
          n.optimal       <- poss.n[optimal]
          a.optimal       <- possible.opt.param[optimal, ]
          scenarios.optimal <- scenarios
          scenarios.optimal[, c(3*(D - 1) + 4, 3*(D - 1) + 5)] <- n.optimal*scenarios.optimal[, c(3*(D - 1) + 4, 3*(D - 1) + 5)]
          I.1.optimal            <- numeric(L)
          I.1.optimal[L]         <- n.optimal*L/(2*sigma.e^2)
          I.1.optimal[1:(L - 1)] <- seq_len(L - 1)*I.1.optimal[L]/L
          I.optimal              <- numeric(2*(D - 1)*L)
          for (l in 1:L){
            I.optimal[(1 + (l - 1)*2*(D - 1)):(l*2*(D - 1))] <- I.1.optimal[l]
          }
        } else if (any(is.numeric(scores)) & any(is.na(scores))){
          optimal <- which(is.na(scores))
          n.optimal       <- poss.n[optimal]
          a.optimal       <- possible.opt.param[optimal, ]
          scenarios.optimal <- scenarios
          scenarios.optimal[, c(3*(D - 1) + 4, 3*(D - 1) + 5)] <- n.optimal*scenarios.optimal[, c(3*(D - 1) + 4, 3*(D - 1) + 5)]
          I.1.optimal            <- numeric(L)
          I.1.optimal[L]         <- n.optimal*L/(2*sigma.e^2)
          I.1.optimal[1:(L - 1)] <- seq_len(L - 1)*I.1.optimal[L]/L
          I.optimal              <- numeric(2*(D - 1)*L)
          for (l in 1:L){
            I.optimal[(1 + (l - 1)*2*(D - 1)):(l*2*(D - 1))] <- I.1.optimal[l]
          }
        } else {
          scenarios.optimal <- NA
          n.optimal <- NA
          a.optimal <- NA
        }
      }
    }
  }

  if (summary == TRUE){
    print("Outputting...")
  }

  output <- list(a = as.vector(a.optimal), a.exact = as.vector(a.exact),
                 all.run.information = as.data.frame(all.run.information),
                 alpha = alpha, B = B, beta = beta, D = D, delta = delta,
                 I = I.optimal, I.1 = I.1.optimal,  I.exact = I.exact,
                 I.1.exact = I.1.exact, init = init,
                 init.file.ext = init.file.ext, L = L, Lambda = Lambda,
                 lower = lower, max.calls = max.calls, n = n.optimal,
                 n.exact = as.numeric(n.exact), num.runs = num.runs,
                 opt.crit = opt.crit, penalty = penalty,
                 scenarios = scenarios.optimal, scenarios.exact = scenarios.exact,
                 seed = seed, sequence.type = sequence.type,  sigma.e = sigma.e,
                 upper = upper, w = w)
  class(output) <- "gsBE"
  return(output)
}
