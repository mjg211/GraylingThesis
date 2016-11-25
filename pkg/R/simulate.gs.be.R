simulate.gs.be <- function(D = 3, L = 3, sigma.e = sqrt(log(0.3^2 + 1)),
                           sigma.b = sqrt(2*sigma.e^2), n = 12,
                           a = c(-0.413, 0.623, 1.899), B = log(1.25),
                           theta = 0, sequence.type = "williams",
                           REML = TRUE, adjust = TRUE, mu.0 = 0,
                           pi = rep(0, D - 1), parallel = TRUE, cpus = 8,
                           replicates = 10000, summary = TRUE){

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
    stop("Between person standard deviation sigma.e must be strictly positive.")
  }
  if ((n%%1 != 0) | (n < 1)){
    stop("n must be a whole number greater than or equal to 1.")
  }
  if (!is.vector(a) | (length(a) != L)){
    stop("a must be a vector of length L.")
  }
  if (B <= 0){
    stop("BE margin B must be strictly positive.")
  }
  if (!is.vector(theta)){
    stop("theta must be a vector.")
  }
  if (!(sequence.type %in% c("latin", "williams"))){
    stop("sequence.type must be set to \"latin\" or \"williams\".")
  }
  if ((replicates%%1 != 0) | (replicates < 1)){
    stop("replicates must be a whole number greater than or equal to 1.")
  }
  if (!is.logical(REML)){
    stop("REML must be logical.")
  }
  if (!is.logical(adjust)){
    stop("adjust must be logical.")
  }
  if (!is.numeric(mu.0)){
    stop("mu.0 must be numeric.")
  }
  if (!is.vector(pi) | (length(pi) != D - 1)){
    stop("pi must be a vector of length D - 1.")
  }
  if (sigma.b <= 0){
    stop("Between person standard deviation sigma.b must be strictly positive.")
  }
  if (!is.logical(parallel)){
    stop("parallel must be set to TRUE or FALSE.")
  }
  if ((cpus%%1 != 0) | (cpus < 1)){
    stop("cpus must be a whole number greater than or equal to 1.")
  }
  if ((replicates%%1 != 0) | (replicates < 1)){
    stop("replicates must be a whole number greater than or equal to 1.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }

  ##### FUNCTION INITIALISATION ###############################################

  individualGSBE <- function(D, L, n, a, mu.0, tau, sigma.e, sigma.b, pi,
                             sequences, Sigmas, REML, adjust){
    drugs.rem  <- 0:(D - 1)
    decision.m <- numeric(D - 1)
    decision.p <- numeric(D - 1)
    each.stage <- numeric(L)
    response   <- NULL
    period     <- NULL
    treatment  <- NULL
    individual <- NULL
    pi         <- c(0, pi)
    tau        <- c(0, tau)
    for (l in 1:L){
      each.stage[l]   <- length(drugs.rem)
      sequences.l     <- sequences[[length(drugs.rem)]]
      switches        <- list()
      for (d in 1:length(drugs.rem)){
        switches[[d]] <- which(sequences.l == d - 1)
      }
      for (d in 1:length(drugs.rem)){
        sequences.l[switches[[d]]] <- drugs.rem[d]
      }

      period         <- c(period, rep(1:length(drugs.rem), n))
      treatment.l    <- rep(as.vector(t(sequences.l)), n/nrow(sequences.l))
      treatment      <- c(treatment, treatment.l)
      individual.l   <- NULL
      for (i in 1:n){
        individual.l <- c(individual.l, rep(i + (l - 1)*n, length(drugs.rem)))
      }
      individual     <- c(individual, individual.l)

      mean.l         <- numeric(n*length(drugs.rem))
      for (i in 1:n){
        for (j in 1:length(drugs.rem)){
          mean.l[(i - 1)*length(drugs.rem) + j] <- mu.0 + pi[j] +
            tau[treatment.l[(i - 1)*length(drugs.rem) + j] + 1]
        }
      }
      response.l     <- rmvnorm(1, mean = mean.l, sigma = Sigmas[[length(drugs.rem)]])
      response       <- c(response, response.l)

      df.analysis        <- data.frame(Response = response,
                                       Period = factor(period, levels = 1:D),
                                       Treatment = factor(treatment, levels = 0:(D - 1)),
                                       Individual = factor(individual, unique(individual)))
      interim.analysis   <- lmer(Response ~ Period + Treatment + (1|Individual),
                                 data = df.analysis, REML = REML)
      beta.hat           <- fixef(interim.analysis)
      tau.hat            <- beta.hat[(D + 1):(2*D - 1)]
      I.hat              <- 1/diag(as.matrix(vcov(interim.analysis)))[(D + 1):(2*D - 1)]
      Z.hat.minus        <- (tau.hat + B)*sqrt(I.hat)
      Z.hat.plus         <- (tau.hat - B)*sqrt(I.hat)
      if (adjust == TRUE){
        DF.l <- n*sum(each.stage) - l*n  - 2*(D - 1)
        a[l] <- qt(pnorm(a[l], lower.tail = FALSE), lower.tail = FALSE, df = DF.l)
      }
      new.drugs.rem     <- NULL
      for (d in drugs.rem[-1]){
        if (Z.hat.minus[d] > a[l] & Z.hat.plus[d] < -a[l]){
          new.drugs.rem <- c(new.drugs.rem, d)
        }
      }
      if (l == L){
        for (d in drugs.rem[-1]){
          if (Z.hat.minus[d] > a[l]){
            decision.m[d] <- 1
          }
          if (Z.hat.plus[d] < -a[l]){
            decision.p[d] <- 1
          }
        }
      }
      drugs.rem         <- c(0, new.drugs.rem)
      if (length(drugs.rem) == 1){
        break
      }
    }
    N   <- n*max(which(each.stage > 0))
    O   <- n*sum(each.stage)
    R.1.minus <- as.numeric(decision.m[1] == 1)
    R.A.minus <- as.numeric(any(decision.m == 1))
    R.1       <- as.numeric(decision.m[1] == 1 & decision.p[1] == 1)
    R.A       <- as.numeric(any(decision.m == 1 & decision.p == 1))
    return(c(R.1.minus, R.A.minus, R.1, R.A, N, O))
  }

  wrapper <- function(rep){
    result <- individualGSBE(D, L, n, a, mu.0, rep(theta[i], D - 1),
                             sigma.e, sigma.b, pi, sequences, Sigmas, REML,
                             adjust)
    return(result)
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  sequences              <- list()
  for (d in 2:D){
    if (sequence.type == "williams"){
      sequences[[d]]     <- williams(d) - 1
    } else {
      sequences.d        <- matrix(0, d, d)
      sequences.d[1, ]   <- 0:(d - 1)
      for (i in 2:d){
        sequences.d[i, ] <- c(sequences[i - 1, 2:d], sequences[i - 1, 1])
      }
      sequences[[d]]     <- sequences.d
    }
  }

  Sigmas        <- list()
  for (d in 2:D){
    Sigmas.d    <- matrix(0, d*n, d*n)
    Sigma.ind   <- matrix(sigma.b^2, d, d) + diag(sigma.e^2, d, d)
    for (i in 1:n){
      Sigmas.d[(1 + (i - 1)*d):(i*d), (1 + (i - 1)*d):(i*d)] <- Sigma.ind
    }
    Sigmas[[d]] <- Sigmas.d
  }

  if (summary == TRUE){
    print("Determining performance across theta...")
  }

  performance <- matrix(0, nrow = length(theta), ncol = 6)
  sink("NULL")
  for (i in 1:length(theta)){
    suppressMessages(sfInit(parallel = parallel, cpus = cpus))
    suppressMessages(sfLibrary(mvtnorm))
    suppressMessages(sfLibrary(lme4))
    suppressMessages(sfLibrary(stats))
    sfExport("D", "L", "n", "a", "mu.0", "theta", "B",
             "sigma.e", "sigma.b", "pi", "sequences",
             "Sigmas", "REML", "adjust", "i", "individualGSBE")
    results <- sfLapply(1:replicates, wrapper)
    suppressMessages(sfStop())
    for (j in 1:replicates){
      performance[i, ] <- performance[i, ] + results[[j]]
    }
    performance[i, ]   <- performance[i, ]/replicates
  }
  sink()
  performance           <- cbind(theta, performance)
  colnames(performance) <- c("theta",
                             "P(Reject H_01- | theta)",
                             "P(Reject H_0d- for some d | theta)",
                             "P(Reject H_01 | theta)",
                             "P(Reject H_0d for some d | theta)",
                             "E(N | theta)",
                             "E(O | theta)")

  if (summary == TRUE){
    print("Outputting...")
  }

  output <- list(a = a, adjust = adjust, B = B, cpus = cpus, D = D, L = L,
                 mu.0 = mu.0, n = n, parallel = parallel,
                 performance = performance, pi = pi, REML = REML,
                 replicates = replicates, sequence.type = sequence.type,
                 sigma.b = sigma.b, sigma.e = sigma.e, theta = theta)
  return(output)
}
