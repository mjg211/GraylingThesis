ssre.be <- function(method = "A", p = 0.5, D = 4,
                    sigma.e.tilde = sqrt(log(0.3^2 + 1)),
                    sigma.e = sqrt(log(0.3^2 + 1)),
                    sigma.b = sqrt(2*sigma.e^2), alpha = 0.05, beta = 0.2,
                    B = log(1.25), delta = 0, theta = 0, pi = rep(0, D - 1),
                    n.max = 48, current = FALSE, REML = TRUE,
                    t.boundaries = TRUE, parallel = TRUE, cpus = 8,
                    replicates = 10000, summary = TRUE){

  ##### ERROR CHECKING ########################################################

  if (!(method %in% c("A", "B", "C"))){
    stop("method must be set to \"A\", \"B\" or \"C\".")
  }
  if ((p <= 0) | (p >= 1)){
    stop("p must be strictly between 0 and 1.")
  }
  if ((D%%1 != 0) | (D < 2)){
    stop("D must be a whole number greater than or equal to 2.")
  }
  if (sigma.e.tilde <= 0){
    stop("Assumed within person standard deviation sigma.e.tilde must be strictly positive.")
  }
  if (sigma.e <= 0){
    stop("Within person standard deviation sigma.e must be strictly positive.")
  }
  if (sigma.b <= 0){
    stop("Between person standard deviation sigma.e must be strictly positive.")
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
  if (!is.numeric(theta)){
    stop("theta must be numeric.")
  }
  if (!is.vector(pi) | (length(pi) != D - 1)){
    stop("pi must be a vector of length D - 1.")
  }
  if ((n.max%%1 != 0) | (n.max < 2)){
    stop("n.max must be a whole number greater than or equal to 2.")
  }
  if (!is.logical(current)){
    stop("current must be logical.")
  }
  if (!is.logical(REML)){
    stop("REML must be logical.")
  }
  if (!is.logical(t.boundaries)){
    stop("t.boundaries must be logical.")
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

  findAlpha <- function(r, alpha, sigma){
    return(((1 - pmvnorm(lower = rep(-r, 2), upper = rep(r, 2),
                         mean = rep(0, 2), sigma = sigma)[1]) - alpha)^2)
  }

  singleTrial <- function(rep, n.be){
    periods.1   <- periods[[1]]
    subjects.1  <- subjects[[1]]
    responses.1 <- as.numeric(rnorm(n.be*D)%*%chol.Sigma.mats[[1]]) +
      means.1
    df.analysis <- data.frame(Period = factor(periods.1, levels = 1:D),
                              Subject = factor(subjects.1, levels = unique(subjects.1)),
                              Treatment = factor(treatments.1, levels = 0:(D - 1)),
                              Response = responses.1)
    interim.analysis <- lmer(Response ~ Period + Treatment + (1 | Subject),
                             data = df.analysis, REML = REML)
    est.sigma.b        <- sqrt(VarCorr(interim.analysis)[[1]][1])
    est.sigma.e        <- sigma(interim.analysis)

    drugs.rem <- 0:(D - 1)
    decision  <- numeric(D - 1)
    Under.P   <- 0
    Over.P    <- 0
    if (method == "A"){
      if (t.boundaries == TRUE){
        r <- qmvt(1 - alpha, df = n.be*D - n.be - 2*(D - 1),
                  sigma = matrix(0.5, D - 1, D - 1) +
                    diag(0.5, D - 1, D - 1))$quantile
      } else {
        r   <- qmvnorm(1 - alpha,
                       sigma = matrix(0.5, D - 1, D - 1) +
                         diag(0.5, D - 1, D - 1))$quantile
      }
      power <- pmvnorm(lower = c(r, -Inf), upper = c(Inf, -r),
                       mean = c(delta + B, delta - B)*sqrt(n.be/(2*est.sigma.e^2)),
                       sigma = matrix(1, 2, 2))[1]
      if (power >= 1 - beta){
        Over.P             <- 1
        beta.hat           <- fixef(interim.analysis)
        tau.hat            <- beta.hat[(D + 1):(2*D - 1)]
        I.hat              <- 1/diag(as.matrix(vcov(interim.analysis)))[(D + 1):(2*D - 1)]
        Z.hat.plus         <- (tau.hat + B)*sqrt(I.hat)
        Z.hat.minus        <- (tau.hat - B)*sqrt(I.hat)
        for (d in drugs.rem[-1]){
          if (Z.hat.plus[d] >= r & Z.hat.minus[d] <= -r){
            decision[d] <- 1
          }
        }
        return(c(decision, as.numeric(any(decision == 1)), n.be, n.be*D, Under.P, Over.P,
                 est.sigma.e, est.sigma.b, est.sigma.e, est.sigma.b, r))
      } else {
        if (delta == 0){
          reest.n.be <- ceiling(2*(qmvnorm(1 - alpha,
                                           sigma = matrix(0.5, D - 1, D - 1) +
                                             diag(0.5, D - 1, D - 1))$quantile +
                                     qnorm(1 - beta/2))^2*est.sigma.e^2/B^2)
        } else {
          reest.n.be <- ceiling(2*(qmvnorm(1 - alpha,
                                           sigma = matrix(0.5, D - 1, D - 1) +
                                             diag(0.5, D - 1, D - 1))$quantile +
                                     qnorm(1 - beta))^2*est.sigma.e^2/(delta - B)^2)
        }
        while (reest.n.be %% D != 0){
          reest.n.be <- reest.n.be + 1
        }
        add.n <- reest.n.be - n.be
        if (add.n + n.be > n.max){
          add.n <- which(1:(n.max - n.be) %% length(drugs.rem) == 0)
          add.n <- add.n[length(add.n)]
          Under.P <- 1
        }
        if (add.n == 0){
          add.n <- length(drugs.rem)
        }
      }
    } else if (method == "B"){
      beta.hat           <- fixef(interim.analysis)
      tau.hat            <- beta.hat[(D + 1):(2*D - 1)]
      I.hat              <- 1/diag(as.matrix(vcov(interim.analysis)))[(D + 1):(2*D - 1)]
      Z.hat.plus         <- (tau.hat + B)*sqrt(I.hat)
      Z.hat.minus        <- (tau.hat - B)*sqrt(I.hat)
      if (t.boundaries == TRUE){
        r <- qmvt(1 - alpha.BC, df = n.be*D - n.be - 2*(D - 1),
                  sigma = matrix(0.5, D - 1, D - 1) +
                    diag(0.5, D - 1, D - 1))$quantile
      } else {
        r   <- qmvnorm(1 - alpha.BC,
                       sigma = matrix(0.5, D - 1, D - 1) +
                         diag(0.5, D - 1, D - 1))$quantile
      }
      new.drugs.rem <- NULL
      for (d in drugs.rem[-1]){
        if (Z.hat.plus[d] >= r & Z.hat.minus[d] <= -r){
          decision[d] <- 1
        } else {
          new.drugs.rem <- c(new.drugs.rem, d)
        }
      }
      drugs.rem <- c(0, new.drugs.rem)
      if (length(drugs.rem) == 1){
        return(c(decision, as.numeric(any(decision == 1)), n.be, n.be*D, Under.P, Over.P,
                 est.sigma.e, est.sigma.b, est.sigma.e, est.sigma.b, r))
      } else {
        power <- pmvnorm(lower = c(r, -Inf), upper = c(Inf, -r),
                         mean = c(delta + B, delta - B)*sqrt(n.be/(2*est.sigma.e^2)),
                         sigma = matrix(1, 2, 2))[1]
        if (power >= 1 - beta){
          Over.P <- 1
          return(c(decision, as.numeric(any(decision == 1)), n.be, n.be*D, Under.P, Over.P,
                   est.sigma.e, est.sigma.b, est.sigma.e, est.sigma.b, r))
        } else {
          if (delta == 0){
            reest.n.be <- ceiling(2*(qmvnorm(1 - alpha.BC,
                                             sigma = matrix(0.5, D - 1, D - 1) +
                                               diag(0.5, D - 1, D - 1))$quantile +
                                       qnorm(1 - beta/2))^2*est.sigma.e^2/B^2)
          } else {
            reest.n.be <- ceiling(2*(qmvnorm(1 - alpha.BC,
                                             sigma = matrix(0.5, D - 1, D - 1) +
                                               diag(0.5, D - 1, D - 1))$quantile +
                                       qnorm(1 - beta))^2*est.sigma.e^2/(delta - B)^2)
          }
          add.n <- reest.n.be - n.be
          if (add.n <= 0){
            Over.P <- 1
            return(c(decision, as.numeric(any(decision == 1)), n.be, n.be*D, Under.P, Over.P,
                     est.sigma.e, est.sigma.b, est.sigma.e, est.sigma.b, r))
          } else {
            while (add.n %% length(drugs.rem) != 0){
              add.n <- add.n + 1
            }
            if (add.n + n.be > n.max){
              add.n <- which(1:(n.max - n.be) %% length(drugs.rem) == 0)
              add.n <- add.n[length(add.n)]
              Under.P <- 1
            }
            if (add.n == 0){
              add.n <- length(drugs.rem)
            }
          }
        }
      }
    } else if (method == "C"){
      if (t.boundaries == TRUE){
        r <- qmvt(1 - alpha, df = n.be*D - n.be - 2*(D - 1),
                  sigma = matrix(0.5, D - 1, D - 1) +
                    diag(0.5, D - 1, D - 1))$quantile
      } else {
        r   <- qmvnorm(1 - alpha,
                       sigma = matrix(0.5, D - 1, D - 1) +
                         diag(0.5, D - 1, D - 1))$quantile
      }
      power <- pmvnorm(lower = c(r, -Inf), upper = c(Inf, -r),
                       mean = c(delta + B, delta - B)*sqrt(n.be/(2*est.sigma.e^2)),
                       sigma = matrix(1, 2, 2))
      if (power >= 1 - beta){
        Over.P             <- 1
        beta.hat           <- fixef(interim.analysis)
        tau.hat            <- beta.hat[(D + 1):(2*D - 1)]
        I.hat              <- 1/diag(as.matrix(vcov(interim.analysis)))[(D + 1):(2*D - 1)]
        Z.hat.plus         <- (tau.hat + B)*sqrt(I.hat)
        Z.hat.minus        <- (tau.hat - B)*sqrt(I.hat)
        for (d in drugs.rem[-1]){
          if (Z.hat.plus[d] >= r & Z.hat.minus[d] <= -r){
            decision[d] <- 1
          }
        }
        return(c(decision, as.numeric(any(decision == 1)), n.be, n.be*D, Under.P, Over.P,
                 est.sigma.e, est.sigma.b, est.sigma.e, est.sigma.b, r))
      } else {
        if (t.boundaries == TRUE){
          r <- qmvt(1 - alpha.BC, df = n.be*D - n.be - 2*(D - 1),
                    sigma = matrix(0.5, D - 1, D - 1) +
                      diag(0.5, D - 1, D - 1))$quantile
        } else {
          r   <- qmvnorm(1 - alpha.BC,
                         sigma = matrix(0.5, D - 1, D - 1) +
                           diag(0.5, D - 1, D - 1))$quantile
        }
        beta.hat           <- fixef(interim.analysis)
        tau.hat            <- beta.hat[(D + 1):(2*D - 1)]
        I.hat              <- 1/diag(as.matrix(vcov(interim.analysis)))[(D + 1):(2*D - 1)]
        Z.hat.plus         <- (tau.hat + B)*sqrt(I.hat)
        Z.hat.minus        <- (tau.hat - B)*sqrt(I.hat)
        new.drugs.rem <- NULL
        for (d in drugs.rem[-1]){
          if (Z.hat.plus[d] >= r & Z.hat.minus[d] <= -r){
            decision[d] <- 1
          } else {
            new.drugs.rem <- c(new.drugs.rem, d)
          }
        }
        drugs.rem <- c(0, new.drugs.rem)
        if (length(drugs.rem) == 1){
          return(c(decision, as.numeric(any(decision == 1)), n.be, n.be*D, Under.P, Over.P,
                   est.sigma.e, est.sigma.b, est.sigma.e, est.sigma.b, r))
        } else {
          if (delta == 0){
            reest.n.be <- ceiling(2*(qmvnorm(1 - alpha.BC,
                                             sigma = matrix(0.5, D - 1, D - 1) +
                                               diag(0.5, D - 1, D - 1))$quantile +
                                       qnorm(1 - beta/2))^2*est.sigma.e^2/B^2)
          } else {
            reest.n.be <- ceiling(2*(qmvnorm(1 - alpha.BC,
                                             sigma = matrix(0.5, D - 1, D - 1) +
                                               diag(0.5, D - 1, D - 1))$quantile +
                                       qnorm(1 - beta))^2*est.sigma.e^2/(delta - B)^2)
          }
          add.n <- reest.n.be - n.be
          if (add.n <= 0){
            Over.P <- 1
            return(c(decision, as.numeric(any(decision == 1)), n.be, n.be*D, Under.P, Over.P,
                     est.sigma.e, est.sigma.b, est.sigma.e, est.sigma.b, r))
          } else {
            while (add.n %% length(drugs.rem) != 0){
              add.n <- add.n + 1
            }
            if (add.n + n.be > n.max){
              add.n <- which(1:(n.max - n.be) %% length(drugs.rem) == 0)
              add.n <- add.n[length(add.n)]
              Under.P <- 1
            }
            if (add.n == 0){
              add.n <- length(drugs.rem)
            }
          }
        }
      }
    }
    R <- length(drugs.rem)
    latin.square        <- matrix(0, R, R)
    latin.square[1, ]   <- drugs.rem
    for (i in 2:R){
      latin.square[i, ] <- c(latin.square[i - 1, 2:R],
                             latin.square[i - 1, 1])
    }
    treatments.2 <- numeric(R^2)
    means.2      <- numeric(R^2)
    for (s in 1:R){
      treatment.s      <- numeric(R)
      mean.s           <- numeric(R)
      for (p in 1:R){
        treatment.s[p] <- latin.square[s, p]
        if (p > 1){
          mean.s[p]    <- mean.s[p] + pi.j[p - 1]
        }
        if (treatment.s[p] != 0){
          mean.s[p]    <- mean.s[p] + tau[treatment.s[p]]
        }
      }
      treatments.2[(1 + (s - 1)*R):(s*R)] <- treatment.s
      means.2[(1 + (s - 1)*R):(s*R)]      <- mean.s
    }
    periods.2  <- periods[[R]][[add.n]]
    subjects.2 <- subjects[[R]][[add.n]]
    treatments.2  <- rep(treatments.2, add.n/R)
    means.2       <- rep(means.2, add.n/R)
    responses.2  <- as.numeric(rnorm(add.n*R)%*%chol.Sigma.mats[[R]][[add.n]]) +
                        means.2
    df.analysis        <- data.frame(Response = c(responses.1, responses.2),
                                     Period = factor(c(periods.1, periods.2), levels = 1:D),
                                     Treatment = factor(c(treatments.1, treatments.2),
                                                        levels = 0:(D - 1)),
                                     Subject = factor(c(subjects.1, subjects.2),
                                                      levels = unique(c(subjects.1, subjects.2))))
    final.analysis <- lmer(Response ~ Period + Treatment + (1 | Subject),
                           data = df.analysis, REML = REML)
    f.est.sigma.b        <- sqrt(VarCorr(final.analysis)[[1]][1])
    f.est.sigma.e        <- sigma(final.analysis)
    beta.hat           <- fixef(final.analysis)
    tau.hat            <- beta.hat[(D + 1):(2*D - 1)]
    I.hat              <- 1/diag(as.matrix(vcov(final.analysis)))[(D + 1):(2*D - 1)]
    Z.hat.plus         <- (tau.hat + B)*sqrt(I.hat)
    Z.hat.minus        <- (tau.hat - B)*sqrt(I.hat)
    if (method == "A"){
      if (t.boundaries == TRUE){
        r <- qmvt(1 - alpha, df = n.be*D + add.n*R - n.be - add.n - 2*(D - 1),
                  sigma = matrix(0.5, D - 1, D - 1) +
                    diag(0.5, D - 1, D - 1))$quantile
      } else {
        r   <- qmvnorm(1 - alpha,
                       sigma = matrix(0.5, D - 1, D - 1) +
                         diag(0.5, D - 1, D - 1))$quantile
      }
    } else if (method == "B"){
      if (current == TRUE){
        if (t.boundaries == TRUE){
          r <- qmvt(1 - alpha.BC, df = n.be*D + add.n*R - n.be - add.n - 2*(D - 1),
                    sigma = matrix(0.5, length(drugs.rem) - 1, length(drugs.rem) - 1) +
                      diag(0.5, length(drugs.rem) - 1, length(drugs.rem) - 1))$quantile
        } else {
          r   <- qmvnorm(1 - alpha.BC,
                         sigma = matrix(0.5, length(drugs.rem) - 1, length(drugs.rem) - 1) +
                           diag(0.5, length(drugs.rem) - 1, length(drugs.rem) - 1))$quantile
        }
      } else {
        if (t.boundaries == TRUE){
          r <- qmvt(1 - alpha.BC, df = n.be*D + add.n*R - n.be - add.n - 2*(D - 1),
                    sigma = matrix(0.5, D - 1, D - 1) +
                      diag(0.5, D - 1, D - 1))$quantile
        } else {
          r   <- qmvnorm(1 - alpha.BC,
                         sigma = matrix(0.5, D - 1, D - 1) +
                           diag(0.5, D - 1, D - 1))$quantile
        }
      }
    } else if (method == "C"){
      if (current == TRUE){
        if (t.boundaries == TRUE){
          r <- qmvt(1 - alpha.BC, df = n.be*D + add.n*R - n.be - add.n - 2*(D - 1),
                    sigma = matrix(0.5, length(drugs.rem) - 1, length(drugs.rem) - 1) +
                      diag(0.5, length(drugs.rem) - 1, length(drugs.rem) - 1))$quantile
        } else {
          r   <- qmvnorm(1 - alpha.BC,
                         sigma = matrix(0.5, length(drugs.rem) - 1, length(drugs.rem) - 1) +
                           diag(0.5, length(drugs.rem) - 1, length(drugs.rem) - 1))$quantile
        }
      } else {
        if (t.boundaries == TRUE){
          r <- qmvt(1 - alpha.BC, df = n.be*D + add.n*R - n.be - add.n - 2*(D - 1),
                    sigma = matrix(0.5, D - 1, D - 1) +
                      diag(0.5, D - 1, D - 1))$quantile
        } else {
          r   <- qmvnorm(1 - alpha.BC,
                         sigma = matrix(0.5, D - 1, D - 1) +
                           diag(0.5, D - 1, D - 1))$quantile
        }
      }
    }
    for (d in drugs.rem[-1]){
      if (Z.hat.plus[d] >= r & Z.hat.minus[d] <= -r){
        decision[d] <- 1
      }
    }
    return(c(decision, as.numeric(any(decision == 1)), n.be + add.n,
             n.be*D + add.n*length(drugs.rem), Under.P, Over.P,
             est.sigma.e, est.sigma.b, f.est.sigma.e, f.est.sigma.b, r))
  }

  wrapper <- function(rep){
    result <- singleTrial(rep, n.be)
    return(result)
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
  while (n.be %% D != 0){
    n.be <- n.be + 1
  }
  n.be   <- ceiling(p*n.be)
  while (n.be %% D != 0){
    n.be <- n.be + 1
  }
  latin.square        <- matrix(0, D, D)
  latin.square[1, ]   <- 0:(D - 1)
  for (i in 2:D){
    latin.square[i, ] <- c(latin.square[i - 1, 2:D],
                           latin.square[i - 1, 1])
  }
  tau              <- rep(theta, D - 1)
  chol.Sigma.mats  <- list()
  periods          <- list()
  subjects         <- list()
  treatment        <- numeric(D^2)
  mean             <- numeric(D^2)
  for (s in 1:D){
    treatment.s      <- numeric(D)
    mean.s           <- numeric(D)
    for (j in 1:D){
      treatment.s[j] <- latin.square[s, j]
      if (j > 1){
        mean.s[j]    <- mean.s[j] + pi.j[j - 1]
      }
      if (treatment.s[j] != 0){
        mean.s[j]    <- mean.s[j] + tau[treatment.s[j]]
      }
    }
    treatment[(1 + (s - 1)*D):(s*D)] <- treatment.s
    mean[(1 + (s - 1)*D):(s*D)]      <- mean.s
  }
  subject                          <- numeric(n.be*D)
  for (i in 1:n.be){
    subject[(1 + (i - 1)*D):(i*D)] <- rep(i, D)
  }
  periods[[1]]  <- rep(1:D, n.be)
  subjects[[1]] <- subject
  treatments.1  <- rep(treatment, n.be/D)
  means.1       <- rep(mean, n.be/D)
  Sigma.mat     <- matrix(0, n.be*D, n.be*D)
  Sigma.block   <- matrix(sigma.b.true^2, D, D) + diag(sigma.e.true^2, D, D)
  for (i in 1:n.be){
    Sigma.mat[(1 + (i - 1)*D):(i*D), (1 + (i - 1)*D):(i*D)] <- Sigma.block
  }
  chol.Sigma.mats[[1]] <- chol(Sigma.mat)
  for (d in 2:D){
    poss.n               <- which(1:n.max %% d == 0)
    chol.Sigma.mats[[d]] <- list()
    periods[[d]]         <- list()
    subjects[[d]]        <- list()
    Sigma.mat            <- matrix(0, poss.n[length(poss.n)]*d, poss.n[length(poss.n)]*d)
    Sigma.block          <- matrix(sigma.b.true^2, d, d) + diag(sigma.e.true^2, d, d)
    subject              <- numeric(poss.n[length(poss.n)]*d)
    for (i in 1:poss.n[length(poss.n)]){
      Sigma.mat[(1 + (i - 1)*d):(i*d), (1 + (i - 1)*d):(i*d)] <- Sigma.block
      subject[(1 + (i - 1)*d):(i*d)]                          <- rep(i, d)
    }
    periods[[d]][[poss.n[length(poss.n)]]]         <- rep(1:d, poss.n[length(poss.n)])
    subjects[[d]][[poss.n[length(poss.n)]]]        <- subject + n.be
    chol.Sigma.mats[[d]][[poss.n[length(poss.n)]]] <- chol(Sigma.mat)
    for (n in poss.n){
      periods[[d]][[n]]         <- periods[[d]][[poss.n[length(poss.n)]]][1:(n*d)]
      subjects[[d]][[n]]        <- subjects[[d]][[poss.n[length(poss.n)]]][1:(n*d)]
      chol.Sigma.mats[[d]][[n]] <- chol.Sigma.mats[[d]][[poss.n[length(poss.n)]]][1:(n*d), 1:(n*d)]
    }
  }
  if (method != "A"){
    r        <- optim(par = qnorm(1 - alpha), fn = findAlpha, method = "Brent",
                      lower = 0, upper = qnorm(1 - alpha/2), alpha = 2*alpha,
                      sigma = matrix(c(1, sqrt(1/2), sqrt(1/2), 1), 2, 2))$par
    alpha.BC <- pnorm(r, lower.tail = FALSE)
    alpha.BC <- f*alpha.BC
  } else {
    alpha.BC <- NA
  }

  if (summary == TRUE){
    print("Determining performance of SSRE procedure...")
  }

  sink("NULL")
  suppressMessages(sfInit(parallel = parallel, cpus = cpus))
  suppressMessages(sfLibrary(lme4))
  suppressMessages(sfLibrary(stats))
  suppressMessages(sfLibrary(MASS))
  suppressMessages(sfLibrary(mvtnorm))
  sfExport("D", "n.be", "alpha", "alpha.BC", "beta", "delta", "B", "theta",
           "pi.j", "n.max", "REML", "t.boundaries",
           "periods", "treatments.1", "subjects", "means.1",
           "chol.Sigma.mats", "method")
  sfExport("singleTrial")
  results           <- sfLapply(1:replicates, wrapper)
  suppressMessages(sfStop())
  sink()
  results           <- matrix(unlist(results), nrow = replicates, ncol = 10 + D - 1, byrow = TRUE)
  colnames(results) <- c(paste("P(Reject H_0", 1:(D - 1), sep = ""), "P(Reject H_0d for some d)",
                         "N", "O", "Under.P", "Over.P", "int.est.sigma.e", "int.est.sigma.b",
                         "final.est.sigma.e", "final.est.sigma.b", "r")
  av.results        <- c(method, D, n.be, alpha, beta, delta, B, theta,
                         sigma.e, sigma.e.true,
                         sigma.b.true, pi.j, p, n.max, t.boundaries,
                         REML, replicates, colMeans(results))
  names(av.results) <- c("method", "D", "n.be", "alpha", "beta", "delta",
                         "B", "theta", "sigma.e",
                         "sigma.e.true", "sigma.b.true",
                         paste("pi_", 2:D, sep = ""), "p", "n.max", "t.boundaries",
                         "REML", "replicates", colnames(results))

  if (summary == TRUE){
    print("Outputting...")
  }

  output            <- list(all.results = results, alpha = alpha, alpha.BC = alpha.BC,
                            av.results = av.results, B = B, beta = beta, cpus = cpus,
                            current = current, D = D, delta = delta, method = method,
                            n.max = n.max, p = p, parallel = parallel, pi = pi,
                            REML = REML, replicates = replicates, sigma.b = sigma.b,
                            sigma.e = sigma.e, sigma.e.tilde = sigma.e.tilde,
                            t.boundaries = t.boundaries, theta = theta)
  return(output)
}
