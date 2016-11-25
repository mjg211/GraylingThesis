gs.co <- function(D = 4, L = 3, sigma.e = sqrt(6.51), delta = 1.25,
                  alpha = 0.05, beta = 0.2, sequence.type = "williams",
                  Delta = 0, summary = TRUE){

  ##### ERROR CHECKING ########################################################

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
  if (delta <= 0){
    stop("Clinically relevant difference delta to power for must be strictly positive.")
  }
  if (!(sequence.type %in% c("latin", "williams"))){
    stop("sequence.type must be set to \"latin\" or \"williams\".")
  }
  if (Delta > 0.5){
    stop("Only Delta less than or equal to 0.5 is supported.")
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

  scenarioProbabilities <- function(D, L, delta, e, f, I, Lambda, scenarios){
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
        } else if ((scenarios[i, j] > 1 & scenarios[i, j] < L) &
                     (scenarios[i, D - 1 + j] == 0)){
          ind.lower <- c(f[1:(scenarios[i, j] - 1)], -Inf,
                         rep(-Inf, L - scenarios[i, j]))
          ind.upper <- c(e[1:(scenarios[i, j] - 1)], f[scenarios[i, j]],
                         rep(Inf, L - scenarios[i, j]))
        } else if ((scenarios[i, j] > 1 & scenarios[i, j] < L) &
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
                                             mean = rep(0,
                                                        L*(D - 1))*sqrt(I),
                                             sigma = Lambda)[1]
      scenarios[i, 2*(D - 1) + 7] <- pmvnorm(lower = lower,
                                             upper = upper,
                                             mean = rep(delta,
                                                        L*(D - 1))*sqrt(I),
                                             sigma = Lambda)[1]
    }
    return(scenarios)
  }

  powerFamilyCO <- function(parameters, D, L, alpha, beta, delta, Delta,
                            Lambda, scenarios){
    Cf             <- parameters[1]
    Ce             <- parameters[2]
    I.1            <- numeric(L)
    I.1[L]         <- ((Ce + Cf)/delta)^2
    I.1[1:(L - 1)] <- seq_len(L - 1)*I.1[L]/L
    I              <- numeric((D - 1)*L)
    for (i in 1:L){
      I[(1 + (i - 1)*(D - 1)):(i*(D - 1))] <- I.1[i]
    }
    e         <- Ce*(seq_len(L)/L)^(Delta - 0.5)
    f         <- delta*sqrt(I.1) - Cf*(seq_len(L)/L)^(Delta - 0.5)
    scenarios <- scenarioProbabilities(D, L, delta, e, f, I, Lambda, scenarios)
    fwer      <- sum(scenarios[, 2*(D - 1) + 2]*scenarios[, 2*(D - 1) + 6])
    type.II   <- 1 - sum(scenarios[, 2*(D - 1) + 3]*scenarios[, 2*(D - 1) + 7])
    return((alpha - fwer)^2 + (beta - type.II)^2)
  }

  informationCO <- function(n, D, L, sigma.e){
    I.1            <- numeric(L)
    I.1[L]         <- n*L/(2*sigma.e^2)
    if (L > 1){
      I.1[1:(L - 1)] <- seq_len(L - 1)*I.1[L]/L
    }
    I              <- numeric((D - 1)*L)
    for (j in 1:L){
      I[(1 + (j - 1)*(D - 1)):(j*(D - 1))] <- I.1[j]
    }
    return(list(I = I, I.1 = I.1))
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  all.omega     <- permutations(n = L, r = D - 1, repeats.allowed = TRUE)
  all.psi       <- permutations(n = 2, r = D - 1, v = 0:1,
                                repeats.allowed = TRUE)
  all.scenarios <- matrix(0, nrow = L^(D - 1)*2^(D - 1), ncol = 2*(D - 1))
  for (i in 1:(L^(D - 1))){
    all.scenarios[(1 + (i - 1)*2^(D - 1)):(i*2^(D - 1)),
                  ] <- cbind(matrix(all.omega[i, ], nrow = 2^(D - 1),
                                    ncol = D - 1, byrow = TRUE),
                             all.psi)
  }
  degeneracy <- numeric(nrow(all.scenarios))
  degenerate <- NULL
  for (i in 1:(nrow(all.scenarios) - 1)){
    for (j in (i + 1):nrow(all.scenarios)){
      if (sum(c(sort(all.scenarios[i, 1:(D - 1)]),
                sort(all.scenarios[i, D:(2*(D - 1))])) ==
                c(sort(all.scenarios[j, 1:(D - 1)]),
                  sort(all.scenarios[j, D:(2*(D - 1))]))) == 2*(D - 1)){
        if (sum(sort(all.scenarios[i, which(all.scenarios[i, D:(2*(D - 1))] == 1)]) ==
                  sort(all.scenarios[j, which(all.scenarios[j, D:(2*(D - 1))] == 1)])) ==
              length(which(all.scenarios[i, D:(2*(D - 1))] == 1))){
          degeneracy[i] <- degeneracy[i] + 1
          degenerate    <- c(degenerate, j)
        }
      }
    }
  }
  degeneracy   <- degeneracy + 1
  if (!is.null(degenerate)){
    degeneracy <- degeneracy[-unique(degenerate)]
    scenarios  <- cbind(all.scenarios[-unique(degenerate), ], degeneracy)
  } else {
    scenarios  <- cbind(all.scenarios, rep(1, nrow(all.scenarios)))
  }
  if (D == 2){
    degeneracy.1 <- numeric(nrow(scenarios))
    for (i in 1:nrow(scenarios)){
      if (scenarios[i, D] == 1){
        degeneracy.1[i] <- 1
      }
    }
    degeneracy.any <- degeneracy.1
  } else {
    degeneracy.any <- numeric(nrow(scenarios))
    degeneracy.1   <- numeric(nrow(scenarios))
    for (i in 1:nrow(scenarios)){
      scenarios[i, 1:(2*(D - 1))] <- c(rev(scenarios[i, 1:(D - 1)]),
                                       rev(scenarios[i, D:(2*(D - 1))]))
    }
    for (i in 1:nrow(scenarios)){
      if (any(scenarios[i, D:(2*(D - 1))] == 1)){
        degeneracy.any[i]    <- scenarios[i, 2*(D - 1) + 1]
        poss.ind.omega.perms <- NULL
        poss.ind.psi.perms   <- NULL
        drop.eff             <- which(scenarios[i, D:(2*(D - 1))] == 1)
        all.poss.scenarios   <- NULL
        for (j in drop.eff){
          poss.ind.omega.perms <- as.matrix(unique(permutations(n = D - 2, r = D - 2,
                                                                v = scenarios[i, -c(j, D:(2*(D - 1)))],
                                                                set = FALSE,
                                                                repeats.allowed = FALSE)))
          poss.ind.psi.perms   <- as.matrix(unique(permutations(n = D - 2, r = D - 2,
                                                                v = scenarios[i, -c(1:(D - 1), D - 1 + j)],
                                                                set = FALSE,
                                                                repeats.allowed = FALSE)))
          poss.scenarios       <- matrix(0, nrow = nrow(poss.ind.omega.perms)*
                                                     nrow(poss.ind.psi.perms),
                                         ncol = 2*(D - 2))
          for (k in 1:nrow(poss.ind.omega.perms)){
            poss.scenarios[(1 + (k - 1)*nrow(poss.ind.psi.perms)):
                             (k*nrow(poss.ind.psi.perms)), ] <- cbind(matrix(poss.ind.omega.perms[k, ],
                                                                             nrow = nrow(poss.ind.psi.perms),
                                                                             ncol = D - 2, byrow = T),
                                                                      poss.ind.psi.perms)
          }
          if (nrow(poss.scenarios) > 1){
            poss.scenarios   <- cbind(rep(scenarios[i, j], nrow(poss.scenarios)),
                                      poss.scenarios[, 1:(D - 2)],
                                      rep(scenarios[i, D - 1 + j],
                                          nrow(poss.scenarios)),
                                      poss.scenarios[, (D - 1):(2*(D - 2))])
          } else {
            poss.scenarios   <- c(rep(scenarios[i, j], nrow(poss.scenarios)),
                                  poss.scenarios[, 1:(D - 2)],
                                  rep(scenarios[i, D - 1 + j],
                                      nrow(poss.scenarios)),
                                  poss.scenarios[, (D - 1):(2*(D - 2))])
          }
          all.poss.scenarios <- rbind(all.poss.scenarios, poss.scenarios)
        }
        if (nrow(all.poss.scenarios) > 1){
          degenerate         <- NULL
          for (k in 1:(nrow(all.poss.scenarios) - 1)){
            for (j in (k + 1):nrow(all.poss.scenarios)){
              if (sum(all.poss.scenarios[k, ] ==
                      all.poss.scenarios[j, ]) == 2*(D - 1)){
                degenerate   <- c(degenerate, j)
              }
            }
          }
          if (!is.null(degenerate)){
            all.poss.scenarios <- all.poss.scenarios[-unique(degenerate), ]
          }
        }
        if (is.vector(all.poss.scenarios)){
          if (all.poss.scenarios[D] == 1){
            degeneracy.1[i] <- degeneracy.1[i] + 1
          }
        } else {
          if (nrow(all.poss.scenarios) == 1){
            if (all.poss.scenarios[1, D] == 1){
              degeneracy.1[i] <- degeneracy.1[i] + 1
            }
          } else {
            zero.i <- NULL
            one.i  <- NULL
            for (j in 1:(D - 1)){
              if (scenarios[i, D - 1 + j] == 0){
                zero.i <- c(zero.i, j)
              } else if (scenarios[i, D - 1 + j] == 1){
                one.i <- c(one.i, j)
              }
            }
            for (j in 1:nrow(all.poss.scenarios)){
              if (all.poss.scenarios[j, D] == 1){
                zero.j <- NULL
                one.j  <- NULL
                for (k in 1:(D - 1)){
                  if (all.poss.scenarios[j, D - 1 + k] == 0){
                    zero.j <- c(zero.j, k)
                  } else if (all.poss.scenarios[j, D - 1 + k] == 1){
                    one.j <- c(one.j, k)
                  }
                }
                check <- 1
                if (length(zero.i) != length(zero.j) |
                    length(one.i) != length(one.j)){
                  check <- 0
                }
                if (check == 1){
                  if (sum(sort(scenarios[i, zero.i]) == sort(all.poss.scenarios[j, zero.j])) !=
                      length(zero.i) |
                      sum(sort(scenarios[i, one.i]) == sort(all.poss.scenarios[j, one.j])) !=
                      length(one.i)){
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
  }
  scenarios           <- cbind(scenarios, degeneracy.any, degeneracy.1,
                               matrix(0, nrow = nrow(scenarios), ncol = 4))
  colnames(scenarios) <- c(paste("omega_", 1:(D - 1), sep = ""),
                           paste("psi_", 1:(D - 1), sep = ""),
                           "deg", "deg.any", "deg.1", "N", "O",
                           "P(H0)", "P(H1)")
  for (i in 1:nrow(scenarios)){
    scenarios[i, 2*(D - 1) + 4] <- max(scenarios[i, 1:(D - 1)])
    for (j in 1:L){
      if (length(which(scenarios[i, 1:(D - 1)] >= j)) > 0){
        scenarios[i, 2*(D - 1) + 5] <- scenarios[i, 2*(D - 1) + 5] +
          length(which(scenarios[i, 1:(D - 1)] >= j)) + 1
      }
    }
  }
  Lambda <- covarianceCO(D, L)

  if (summary == TRUE){
    print("Determining exact GS design...")
  }

  if (L == 1){
    # Compute number n must be divisible by in a single stage design
    if ((sequence.type == "williams") & (is.odd(D))){
      lcm <- 2*D
    } else {
      lcm <- D
    }
    # Determine exact required sample size of a single stage design using
    # the Dunnett correction
    e          <- qmvnorm(1 - alpha, sigma = Lambda)$quantile
    f          <- e
    alpha.star <- pnorm(e, lower.tail = FALSE)
    n.exact    <- (qnorm(1 - alpha.star) +
                     qnorm(1 - beta))^2*2*(delta/sigma.e)^(-2)
  } else {
    # Determine number n must be divisible by in an L stage design
    sink("NULL")
    lcm       <- 2
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
    # Determine exact required sample size of an L stage design
    n              <- (qnorm(1 - alpha/(D - 1)) +
                         qnorm(1 - beta))^2*2*(delta/sigma.e)^(-2)
    optimal.values <- optim(par = c(sqrt(n*delta^2/(8*sigma.e^2)),
                                    sqrt(n*delta^2/(8*sigma.e^2))),
                            fn = powerFamilyCO, D = D, L = L, alpha = alpha,
                            beta = beta, delta = delta, Delta = Delta,
                            Lambda = Lambda, scenarios = scenarios)
    Cf             <- optimal.values$par[1]
    Ce             <- optimal.values$par[2]
    I.1            <- numeric(L)
    I.1[L]         <- ((Ce + Cf)/delta)^2
    I.1[1:(L - 1)] <- seq_len(L - 1)*I.1[L]/L
    e              <- Ce*(seq_len(L)/L)^(Delta - 0.5)
    f              <- delta*sqrt(I.1) - Cf*(seq_len(L)/L)^(Delta - 0.5)
    n.exact        <- 2*sigma.e^2*(I.1[1])
  }

  scenarios.exact                  <- scenarios
  scenarios.exact[, 2*(D - 1) + 4] <- n.exact*scenarios.exact[, 2*(D - 1) + 4]
  scenarios.exact[, 2*(D - 1) + 5] <- n.exact*scenarios.exact[, 2*(D - 1) + 5]
  information.exact <- informationCO(n.exact, D, L, sigma.e)
  I.1.exact         <- information.exact$I.1
  I.exact           <- information.exact$I
  scenarios.exact   <- scenarioProbabilities(D, L, delta, e, f, I.exact, Lambda, scenarios.exact)
  type.I.exact      <- sum(scenarios.exact[, 2*(D - 1) + 2]*scenarios.exact[, 2*(D - 1) + 6])
  type.II.exact     <- 1 - sum(scenarios.exact[, 2*(D - 1) + 3]*scenarios.exact[, 2*(D - 1) + 7])
  if (L > 1){
    EN.H0.exact     <- sum(scenarios.exact[, 2*(D - 1) + 1]*
                       scenarios.exact[, 2*(D - 1) + 6]*scenarios.exact[, 2*(D - 1) + 4])
    EN.H1.exact     <- sum(scenarios.exact[, 2*(D - 1) + 1]*
                       scenarios.exact[, 2*(D - 1) + 7]*scenarios.exact[, 2*(D - 1) + 4])
    maxN.exact      <- max(scenarios.exact[, 2*(D - 1) + 4])
    EO.H0.exact     <- sum(scenarios.exact[, 2*(D - 1) + 1]*
                       scenarios.exact[, 2*(D - 1) + 6]*scenarios.exact[, 2*(D - 1) + 5])
    EO.H1.exact     <- sum(scenarios.exact[, 2*(D - 1) + 1]*
                       scenarios.exact[, 2*(D - 1) + 7]*scenarios.exact[, 2*(D - 1) + 5])
    maxO.exact      <- max(scenarios.exact[, 2*(D - 1) + 5])
  } else {
    EN.H0.exact     <- n.exact
    EN.H1.exact     <- n.exact
    maxN.exact      <- n.exact
    EO.H0.exact     <- n.exact*D
    EO.H1.exact     <- n.exact*D
    maxO.exact      <- n.exact*D
  }
  exact.design.performance        <- c(type.I.exact, 1 - type.II.exact,
                                       EN.H0.exact, EN.H1.exact, maxN.exact,
                                       EO.H0.exact, EO.H1.exact, maxO.exact)
  names(exact.design.performance) <- c("FWER", "Power", "E(N|H_N)", "E(N|H_A)",
                                       "maxN", "E(O|H_N)", "E(O|H_A)", "maxO")

  if (summary == TRUE){
    print("Examining rounded GS design...")
  }

  n   <- ceiling(n.exact)
  while (n%%lcm != 0) {
    n <- n + 1
  }
  scenarios[, 2*(D - 1) + 4] <- n*scenarios[, 2*(D - 1) + 4]
  scenarios[, 2*(D - 1) + 5] <- n*scenarios[, 2*(D - 1) + 5]
  information <- informationCO(n, D, L, sigma.e)
  I.1         <- information$I.1
  I           <- information$I
  scenarios   <- scenarioProbabilities(D, L, delta, e, f, I, Lambda, scenarios)
  type.I      <- sum(scenarios[, 2*(D - 1) + 2]*scenarios[, 2*(D - 1) + 6])
  type.II     <- 1 - sum(scenarios[, 2*(D - 1) + 3]*scenarios[, 2*(D - 1) + 7])
  if (L > 1){
    EN.H0     <- sum(scenarios[, 2*(D - 1) + 1]*
                       scenarios[, 2*(D - 1) + 6]*scenarios[, 2*(D - 1) + 4])
    EN.H1     <- sum(scenarios[, 2*(D - 1) + 1]*
                       scenarios[, 2*(D - 1) + 7]*scenarios[, 2*(D - 1) + 4])
    maxN      <- max(scenarios[, 2*(D - 1) + 4])
    EO.H0     <- sum(scenarios[, 2*(D - 1) + 1]*
                       scenarios[, 2*(D - 1) + 6]*scenarios[, 2*(D - 1) + 5])
    EO.H1     <- sum(scenarios[, 2*(D - 1) + 1]*
                       scenarios[, 2*(D - 1) + 7]*scenarios[, 2*(D - 1) + 5])
    maxO      <- max(scenarios[, 2*(D - 1) + 5])
  } else {
    EN.H0     <- n
    EN.H1     <- n
    maxN      <- n
    EO.H0     <- n*D
    EO.H1     <- n*D
    maxO      <- n*D
  }
  design.performance        <- c(type.I, 1 - type.II, EN.H0, EN.H1,
                                 maxN, EO.H0, EO.H1, maxO)
  names(design.performance) <- c("FWER", "Power",
                                 "E(N|H0)", "E(N|H1)", "maxN",
                                 "E(O|H0)", "E(O|H1)", "maxO")

  if (summary == TRUE){
    print("Outputting...")
  }

  output        <- list(alpha = alpha, beta = beta, D = D, Delta = Delta,
                        delta = delta, design.performance = design.performance,
                        exact.design.performance = exact.design.performance,
                        e = e, f = f, I = I, I.1 = I.1, I.exact = I.exact,
                        I.1.exact = I.1.exact, L = L, Lambda = Lambda,
                        n = n, n.exact = n.exact, scenarios = scenarios,
                        scenarios.exact = scenarios.exact,
                        sequence.type = sequence.type, sigma.e = sigma.e)
  class(output) <- "gsCO"
  return(output)
}
