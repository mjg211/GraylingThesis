max.fwer.gs.co <- function(D = 4, L = 3, sigma.e = sqrt(6.51),
                           sigma.b = sqrt(10.12), alpha = 0.05,
                           beta = 0.2, delta = 1.25, sequences,
                           e = c(8.53, 3.18, 1.92), f = c(-0.30, 0.59, 1.92),
                           n.or.m = "m", n.m = 3, max.calls = 1000,
                           lower = NULL, upper = NULL, init = NULL,
                           seed = Sys.time(), summary = TRUE){

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
  if (sigma.b <= 0){
    stop("Between person standard deviation sigma.b must be strictly positive.")
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
  if (!(sequences %in% c("latin", "williams")) & (!is.list(sequences))){
    stop("sequences must be set to \"latin\" or \"williams\", or be a list of the form described in the help documentation.")
  }
  if (is.list(sequences)){
    for (d in 2:D){
      if (is.null(sequences[[d]])){
        stop("sequences must be set to \"latin\" or \"williams\", or be a list of the form described in the help documentation.")
      } else {
        if (sum(unique(sequences[[d]]) == 0:(d - 1)) != d){
          stop("sequences must be set to \"latin\" or \"williams\", or be a list of the form described in the help documentation.")
        }
      }
    }
  }
  if (!(is.vector(e)) | length(e) != L){
    stop("e must be a vector of length L.")
  }
  if (!(is.vector(f)) | length(f) != L){
    stop("f must be a vector of length L.")
  }
  if (!(n.or.m %in% c("n", "m"))){
    stop("n.or.m must be set to \"n\" or \"m\".")
  }
  if ((n.m%%1 != 0) | (n.m < 1)){
    stop("n.m must be a whole number greater than or equal to 1.")
  }
  if ((max.calls%%1 != 0) | (max.calls < 1)){
    stop("max.calls must be a whole number greater than or equal to 1.")
  }
  if (!is.null(lower)){
    if (!(is.vector(lower)) | length(lower) != D - 1){
      stop("lower must be a vector of length D - 1.")
    }
  }
  if (!is.null(upper)){
    if (!(is.vector(upper)) | length(upper) != D - 1){
      stop("upper must be a vector of length D - 1.")
    }
  }
  if (!is.null(init)){
    if (!(is.vector(init)) | length(init) != D - 1){
      stop("init must be a vector of length D - 1.")
    }
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }

  ##### FUNCTION INITIALISATION ###############################################

  latinSquare <- function(D = 4, labels = 0:(D - 1)){
    design           <- matrix(0, D, D)
    design[1, ]      <- labels
    for (k in 2:D){
      design[k, ]    <- c(design[k - 1, 2:D], design[k - 1, 1])
    }
    rownames(design) <- paste("k =", seq_len(D))
    colnames(design) <- paste("j =", seq_len(D))
    return(design)
  }

  williamsSquare <- function(D = 4, labels = 0:(D - 1)){
    if (D == 2){
      prelim.design        <- matrix(c(0, 1, 1, 0), 2, 2)
    } else {
      k.1                  <- c(0, rep(1, D - 1))
      factor               <- ifelse((3:D)%%2 == 0, -1, 1)
      for (j in 3:D){
        k.1[j]             <- k.1[j - 1] + factor[j - 2]*(D + 1 - j)
      }
      prelim.design        <- rbind(k.1, matrix(rep(0, D*(D - 1)), ncol = D))
      for (k in 2:D){
        prelim.design[k, ] <- (prelim.design[k - 1, ] + 1)%%D
      }
      if (D%%2 == 1){
        prelim.design      <- rbind(prelim.design,
                                    t(apply(prelim.design, 1, rev)))
      }
    }
    design <- matrix(0, nrow = nrow(prelim.design), ncol = ncol(prelim.design))
    for (i in 0:(D - 1)){
      design[which(prelim.design == i)] <- labels[i + 1]
    }
    rownames(design) <- paste("k =", seq_len(nrow(design)))
    colnames(design) <- paste("j =", seq_len(D))
    return(design)
  }

  informationLambdaFinder <- function(D, L, drop.scenarios,
                                      XT.Sigma.inv.X.matrices, n.or.m){
    information <- list()
    lambdas     <- list()
    Lambda      <- matrix(0, L*(D - 1), L*(D - 1))
    I           <- numeric(L*(D - 1))
    for (i in 1:nrow(drop.scenarios)){
      cumul.XT  <- matrix(0, 2*(D - 1) + 1, 2*(D - 1) + 1)
      for (j in 1:L){
        if (drop.scenarios[i, j] > 1){
          if (n.or.m == "n"){
            cumul.XT <- cumul.XT + XT.Sigma.inv.X.matrices[[drop.scenarios[i, j]]]/
              nrow(sequences[[drop.scenarios[i, j]]])
          } else {
            cumul.XT <- cumul.XT + XT.Sigma.inv.X.matrices[[drop.scenarios[i, j]]]
          }
        }
        inv.cumul.XT                         <- ginv(cumul.XT)[(D + 1):(2*(D - 1) + 1),
                                                               (D + 1):(2*(D - 1) + 1)]
        I[(1 + (D - 1)*(j - 1)):((D - 1)*j)] <- 1/diag(inv.cumul.XT)
        Lambda[(1 + (D - 1)*(j - 1)):((D - 1)*j),
               (1 + (D - 1)*(j - 1)):((D - 1)*j)] <- sqrt(diag(I[(1 + (D - 1)*(j - 1)):((D - 1)*j)]))%*%
          inv.cumul.XT%*%
          sqrt(diag(I[(1 + (D - 1)*(j - 1)):((D - 1)*j)]))
        if (j > 1){
          for (k in 1:(j - 1)){
            Lambda[(1 + (D - 1)*(j - 1)):((D - 1)*j),
                   (1 + (D - 1)*(k - 1)):((D - 1)*k)] <- sqrt(diag(I[(1 + (D - 1)*(k - 1)):((D - 1)*k)]))%*%
              inv.cumul.XT%*%
              sqrt(diag(I[(1 + (D - 1)*(j - 1)):((D - 1)*j)]))
            Lambda[(1 + (D - 1)*(k - 1)):((D - 1)*k),
                   (1 + (D - 1)*(j - 1)):((D - 1)*j)] <- t(Lambda[(1 + (D - 1)*(j - 1)):((D - 1)*j),
                                                                  (1 + (D - 1)*(k - 1)):((D - 1)*k)])
          }
        }
      }
      information[[i]] <- I
      lambdas[[i]]     <- Lambda
    }
    output <- list(information = information, lambdas = lambdas)
    return(output)
  }

  designMatrices <- function(n, D, sequences, period.effect = TRUE,
                             sigma.e, sigma.b, remaining = 0:(D - 1),
                             s = 1, seed = as.integer(Sys.time())){

    if (is.character(sequences)){
      if (sequences == "williams"){
        sequences <- williamsSquare(length(remaining))
      } else if (sequences == "latin"){
        sequences <- latinSquare(length(remaining))
      }
    }

    p <- ncol(sequences)
    k <- nrow(sequences)

    X       <- matrix(0, nrow = n*p, ncol = 1 +
                        (1 + as.integer(period.effect))*(D - 1))
    X[, 1]  <- 1
    Z       <- matrix(0, nrow = n*p, ncol = n)
    G       <- diag(sigma.b^2, nrow = n, ncol = n)
    R       <- diag(sigma.e^2, nrow = n*p, ncol = n*p)
    subject <- rep(0, n*p)

    for (i in 1:n){
      Z[(1 + p*(i - 1)):(p*i), i]    <- 1
      subject[(1 + p*(i - 1)):(p*i)] <- s + i - 1
    }

    if (period.effect == TRUE){
      period.block                                  <- matrix(0, nrow = p,
                                                              ncol = D - 1)
      period.block[2:p, 1:(p - 1)]                  <- diag(p - 1)
      period.complete                               <- matrix(0, nrow = n*p,
                                                              ncol = D - 1)
      for (i in 1:n){
        period.complete[(1 + p*(i - 1)):(p*i), ] <- period.block
      }
      X[, 2:D] <- period.complete
    }

    to   <- n
    while (floor(to/p) != ceiling(to/p)){
      to <- to - 1
    }
    for (i in 1:(to/k)){
      for (j in 1:k){
        for (l in 1:p){
          if (sequences[j, l] != 0){
            X[k*p*(i - 1) + p*(j - 1) + l,
              1 + as.integer(period.effect)*(D - 1) + sequences[j, l]] <- 1
          }
        }
      }
    }
    if (to != n){
      rand.seq <- sample.int(k, n - to)
      for (j in 1:length(rand.seq)){
        for (l in 1:p){
          if (sequences[rand.seq[j], l] != 0){
            X[to*p + p*(j - 1) + l,
              1 + as.integer(period.effect)*(D - 1) + sequences[rand.seq[j], l]] <- 1
          }
        }
      }
    }

    Sigma          <- Z%*%G%*%t(Z) + R
    XT.Sigma.inv.X <- t(X)%*%ginv(Sigma)%*%X
    Covariance     <- ginv(XT.Sigma.inv.X)
    Lambda <- Covariance[(2 + as.integer(period.effect)*(D - 1)):(D + as.integer(period.effect)*(D - 1)),
                         (2 + as.integer(period.effect)*(D - 1)):(D + as.integer(period.effect)*(D - 1))]
    output <- list(Covariance = Covariance, D = D, G = G, Lambda = Lambda, n = n, period.effect = period.effect,
                   R = R, remaining = remaining, s = s, seed = seed, sequences = sequences,
                   Sigma = Sigma, sigma.e = sigma.e, sigma.b = sigma.b,
                   X = X, XT.Sigma.inv.X = XT.Sigma.inv.X, Z = Z)
    return(output)
  }

  fwer <- function(parameters, D, L, e, f, information, lambdas,
                   scenarios, drop.scenarios){
    tau <- parameters
    if (all(tau > 0)){
      return(0)
    } else {
      local.scenarios <- numeric(nrow(scenarios))
      for (i in 1:nrow(scenarios)){
        if (any(scenarios[i, D:(2*(D - 1))][which(tau <= 0)] == 1)){
          local.scenarios[i] <- 1
        }
      }
      local.scenarios <- scenarios[which(local.scenarios == 1), ]
      ind.lower <- numeric(L)
      ind.upper <- numeric(L)
      lower     <- numeric((D - 1)*L)
      upper     <- numeric((D - 1)*L)
      for (i in 1:nrow(local.scenarios)){
        sorted.scenario <- sort.int(local.scenarios[i, 1:(D - 1)], decreasing = TRUE, index.return = TRUE)
        local.scenarios[i, 1:(2*(D - 1))] <- c(sorted.scenario$x,
                                               local.scenarios[i, D:(2*(D - 1))][sorted.scenario$ix])
        tau.scenario <- tau[sorted.scenario$ix]
        for (j in 1:(D - 1)){
          if ((local.scenarios[i, j] == 1) & (local.scenarios[i, D - 1 + j] == 0)){
            ind.lower <- rep(-Inf, L)
            ind.upper <- c(f[1], rep(Inf, L - 1))
          } else if ((local.scenarios[i, j] == 1) & (local.scenarios[i, D - 1 + j] == 1)){
            ind.lower <- c(e[1], rep(-Inf, L - 1))
            ind.upper <- rep(Inf, L)
          } else if ((local.scenarios[i, j] %in% 2:(L - 1)) & (local.scenarios[i, D - 1 + j] == 0)){
            ind.lower <- c(f[1:(local.scenarios[i, j] - 1)], -Inf,
                           rep(-Inf, L - local.scenarios[i, j]))
            ind.upper <- c(e[1:(local.scenarios[i, j] - 1)], f[local.scenarios[i, j]],
                           rep(Inf, L - local.scenarios[i, j]))
          } else if ((local.scenarios[i, j] %in% 2:(L - 1)) & (local.scenarios[i, D - 1 + j] == 1)){
            ind.lower <- c(f[1:(local.scenarios[i, j] - 1)], e[local.scenarios[i, j]],
                           rep(-Inf, L - local.scenarios[i, j]))
            ind.upper <- c(e[1:(local.scenarios[i, j] - 1)], Inf,
                           rep(Inf, L - local.scenarios[i, j]))
          } else if ((local.scenarios[i, j] == L) & (local.scenarios[i, D - 1 + j] == 0)){
            ind.lower <- c(f[1:(local.scenarios[i, j] - 1)], -Inf)
            ind.upper <- e
          } else if ((local.scenarios[i, j] == L) & (local.scenarios[i, D - 1 + j] == 1)){
            ind.lower <- f
            ind.upper <- c(e[1:(local.scenarios[i, j] - 1)], Inf)
          }
          lower[seq(from = j, to = j + (D - 1)*(L - 1), length.out = L)] <- ind.lower
          upper[seq(from = j, to = j + (D - 1)*(L - 1), length.out = L)] <- ind.upper
        }
        I      <- information[[local.scenarios[i, 2*(D - 1) + 1]]]
        Lambda <- lambdas[[local.scenarios[i, 2*(D - 1) + 1]]]
        local.scenarios[i, 2*(D - 1) + 2] <- pmvnorm(lower = lower,
                                                     upper = upper,
                                                     mean = rep(tau.scenario, L)*sqrt(I),
                                                     sigma = Lambda)[1]
      }
      return(-sum(local.scenarios[, 2*(D - 1) + 2]))
    }
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  set.seed(seed)

  if (is.null(lower)){
    lower   <- rep(-20, D - 1)
  }
  if (is.null(upper)){
    upper   <- rep(20, D - 1)
  }
  if (is.null(init)){
    init      <- numeric(D - 1)
    for (i in 1:length(lower)){
      init[i] <- runif(1, min = lower[i], max = upper[i])
    }
  }

  if (is.character(sequences)){
    if (sequences == "williams"){
      sequences <- list()
      for (d in 2:D){
        sequences[[d]] <- williamsSquare(d)
      }
    } else if (sequences == "latin"){
      sequences <- list()
      for (d in 2:D){
        sequences[[d]] <- latinSquare(d)
      }
    }
  }

  all.drop.scenarios           <- permutations(D, L, repeats.allowed = TRUE)
  sorted.drop.scenarios        <- all.drop.scenarios
  for (i in 1:nrow(all.drop.scenarios)){
    sorted.drop.scenarios[i, ] <- sort(all.drop.scenarios[i, ],
                                       decreasing = TRUE)
  }
  drop.scenarios               <- unique(sorted.drop.scenarios)
  stage.1.D                    <- numeric(nrow(drop.scenarios))
  for (i in 1:nrow(drop.scenarios)){
    if (drop.scenarios[i, 1] == D){
      stage.1.D[i]             <- 1
    }
  }
  drop.scenarios               <- drop.scenarios[as.logical(stage.1.D), ]
  all.omega <- permutations(n = L, r = D - 1, repeats.allowed = TRUE)
  colnames(all.omega) <- paste("omega_", 1:(D - 1), sep = "")
  all.psi   <- permutations(n = 2, r = D - 1, v = 0:1, repeats.allowed = TRUE)
  colnames(all.psi) <- paste("psi_", 1:(D - 1), sep = "")
  all.scenarios <- matrix(0, nrow = L^(D - 1)*2^(D - 1),
                          ncol = 2*(D - 1))
  for (i in 1:(L^(D - 1))){
    all.scenarios[(1 + (i - 1)*2^(D - 1)):(i*2^(D - 1)),
                  ] <- cbind(matrix(all.omega[i, ], nrow = 2^(D - 1),
                                    ncol = D - 1, byrow = TRUE),
                             all.psi)
  }
  scenarios           <- cbind(all.scenarios, matrix(0, nrow = nrow(all.scenarios),
                                                     ncol = 2))
  colnames(scenarios) <- c(paste("omega_", 1:(D - 1), sep = ""),
                           paste("psi_", 1:(D - 1), sep = ""),
                           "drop.sc", "P(tau)")
  for (i in 1:nrow(scenarios)){
    each.stage <- numeric(L)
    for (j in 1:L){
      each.stage[j] <- length(which(scenarios[i, 1:(D - 1)] >= j)) + 1
    }
    for (j in 1:nrow(drop.scenarios)){
      if (sum(drop.scenarios[j, ] == each.stage) == L){
        scenarios[i, 2*(D - 1) + 1] <- j
        break
      }
    }
  }

  XT.Sigma.inv.X.matrices        <- list()
  for (d in 2:D){
    XT.Sigma.inv.X.matrices[[d]] <- designMatrices(D = D,
                                                   sequences = sequences[[d]],
                                                   n = nrow(sequences[[d]]),
                                                   sigma.e = sigma.e,
                                                   sigma.b = sigma.b)$XT.Sigma.inv.X
  }

  information.lambdas <- informationLambdaFinder(D, L, drop.scenarios,
                                                 XT.Sigma.inv.X.matrices,
                                                 n.or.m)
  information         <- information.lambdas$information
  lambdas             <- information.lambdas$lambdas
  for (i in 1:nrow(drop.scenarios)){
    information[[i]] <- n.m*information[[i]]
  }

  if (summary == TRUE){
    print("Searching for maximal FWER. Output from GenSA may follow...")
  }

  max.fwer <- GenSA(par = init, fn = fwer, lower = lower,
                    upper = upper, control = list(max.call = max.calls,
                                                  verbose = summary), D = D,
                    L = L, e = e, f = f, scenarios = scenarios,
                    drop.scenarios = drop.scenarios, information = information,
                    lambdas = lambdas)

  if (summary == TRUE){
    print("Outputting...")
  }

  output <- list(alpha = alpha, beta = beta, D = D, delta = delta, e = e, f = f,
                 fwer = -max.fwer$value, init = init, L = L, lower = lower,
                 max.calls = max.calls, max.fwer.search = max.fwer, n.m = n.m,
                 n.or.m = n.or.m, seed = seed, sequences = sequences,
                 sigma.b = sigma.b, sigma.e = sigma.e, tau = max.fwer$par,
                 upper = upper)
  return(output)
}
