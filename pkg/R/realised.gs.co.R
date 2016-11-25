realised.gs.co <- function(trial.type = "1A", all.data = TRUE, D = 4,
                           sigma.e = sqrt(6.51), sigma.b = sqrt(10.12),
                           sequence.type = "williams", lambda = 1.1, w = 1,
                           p = 6, n.vec = c(36, 72, 108),
                           f = c(0.35, 1.37, 2.08), e = c(2.90, 2.35, 2.08),
                           mu.0 = 10.65, tau = rep(0, D - 1),
                           pi = c(-0.77, -0.96, -0.55), REML = TRUE,
                           parallel = TRUE, cpus = 8, replicates = 10000,
                           seed = Sys.time(), summary = TRUE){

  ##### ERROR CHECKING ########################################################

  if (!(trial.type %in% c("1A", "1B", "2"))){
    stop("trial.type must be set to \"1A\", \"1B\" or \"2\".")
  }
  if ((D%%1 != 0) | (D < 2)){
    stop("D must be a whole number greater than or equal to 2.")
  }
  if (sigma.e <= 0){
    stop("Within person standard deviation sigma.e must be strictly positive.")
  }
  if (sigma.b <= 0){
    stop("Between person standard deviation sigma.e must be strictly positive.")
  }
  if (!(sequence.type %in% c("latin", "williams"))){
    stop("sequence.type must be set to \"latin\" or \"williams\".")
  }
  if (lambda <= 0){
    stop("lambda must be strictly positive.")
  }
  if ((w%%1 != 0) | (w < 1)){
    stop("w must be a whole number greater than or equal to 1.")
  }
  if ((p%%1 != 0) | (p < 1)){
    stop("p must be a whole number greater than or equal to 1.")
  }
  if (any(n.vec%%1 != 0) | any(n.vec < 1)){
    stop("n.vec must be a vector of increasing integer values.")
  }
  for (l in 2:length(n.vec)){
    if (n.vec[l] <= n.vec[l - 1]){
      stop("n.vec must be a vector of increasing integer values.")
    }
  }
  if (!is.vector(e) | (length(e) != length(n.vec))){
    stop("e must be a vector of length length(n.vec).")
  }
  if (!is.vector(f) | (length(f) != length(n.vec))){
    stop("f must be a vector of length length(n.vec).")
  }
  if (!is.numeric(mu.0)){
    stop("mu.0 must be numeric.")
  }
  if (!is.vector(tau) | (length(tau) != D - 1)){
    stop("tau must be a vector of length D - 1.")
  }
  if (!is.vector(pi) | (length(pi) != D - 1)){
    stop("pi must be a vector of length D - 1.")
  }
  if (!is.logical(REML)){
    stop("REML must be logical.")
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

  sequences <- function(D, drugs.rem, sequence.type){
    if (sequence.type == "latin"){
      sequences        <- matrix(0, D, D)
      sequences[1, ]   <- 0:(D - 1)
      for (d in 2:D){
        sequences[d, ] <- c(sequences[d - 1, 2:D], sequences[d - 1, 1])
      }
    } else {
      sequences        <- williams(D) - 1
    }
    switches           <- list()
    for (d in 1:length(drugs.rem)){
      switches[[d]]    <- which(sequences == d - 1)
    }
    for (d in 1:length(drugs.rem)){
      sequences[switches[[d]]] <- drugs.rem[d]
    }
    return(sequences)
  }

  individualX <- function(D, R, sequence){
    X          <- matrix(0, R, 2*D - 1)
    X[, 1]     <- 1
    X[, 2:R]   <- rbind(rep(0, R - 1), diag(1, R - 1, R - 1))
    for (r in 1:R){
      if (sequence[r] != 0){
        X[r, D + sequence[r]] <- 1
      }
    }
    return(X)
  }

  designMatrices <- function(patient.info, current.status){
    Y <- NULL
    X <- NULL
    Z <- NULL
    included.patients <- which(current.status == "C")
    for (i in included.patients){
      Y <- c(Y, patient.info[[i]]$Y)
      X <- rbind(X, patient.info[[i]]$X)
      Z <- rbind(cbind(Z, rep(1, nrow(Z))), cbind(matrix(0, patient.info[[i]]$R, ncol(Z)),
                                                  rep(1, patient.info[[i]]$R)))
    }
    design.matrices   <- list()
    design.matrices$X <- X
    design.matrices$Y <- Y
    design.matrices$Z <- Z
    return(design.matrices)
  }

  analysisData <- function(patient.info, all.data, n, patients.old){
    periods   <- NULL
    subjects   <- NULL
    treatments <- NULL
    responses  <- NULL
    if (all.data == FALSE){
      to.include <- 1:n
    } else {
      to.include <- NULL
      for (i in 1:patients.old){
        if (!is.null(patient.info[[i]]$Y)){
          to.include <- c(to.include, i)
        }
      }
    }
    for (i in to.include){
      periods    <- c(periods, 1:length(patient.info[[i]]$Y))
      subjects   <- c(subjects, rep(i, length(patient.info[[i]]$Y)))
      treatments <- c(treatments, patient.info[[i]]$sequence[1:length(patient.info[[i]]$Y)])
      responses  <- c(responses, patient.info[[i]]$Y)
    }
    df.analysis <- data.frame(period = factor(periods, levels = 1:D),
                              subject = factor(subjects, levels = unique(subjects)),
                              treatment = factor(treatments, levels = 0:(D - 1)),
                              response = responses)
    return(df.analysis)
  }

  singleTrial <- function(rep){
    drugs.rem           <- 0:(D - 1)
    main.sequences      <- list()
    main.sequences[[1]] <- sequences(D, drugs.rem, sequence.type)
    curr.sequences      <- main.sequences[[1]]
    beta           <- c(mu.0, pi, tau)
    R              <- D
    C              <- 0
    L              <- length(n.vec)
    patients.old   <- 0
    observations   <- 0
    waitlist.time  <- 0
    patient.info   <- list()
    current.status <- NULL
    week           <- 1
    stage.dec      <- numeric(D - 1)
    decision       <- numeric(D - 1)
    beta.hat <- matrix(0, L, 2*D - 1)
    I.hat    <- matrix(0, L, D - 1)
    for (l in 1:L){
      if (trial.type == "1B"){
        if (l > 1){
          if (any(stage.dec == l - 1)){
            for (i in 1:patients.old){
              per.i <- patient.info[[i]]$curr.period
              num.per.i <- patient.info[[i]]$R
              if (per.i < num.per.i){
                if (!all(patient.info[[i]]$sequence[(per.i + 1):num.per.i] %in% drugs.rem)){
                  to.replace <- which(!(patient.info[[i]]$sequence[(per.i + 1):num.per.i] %in% drugs.rem))
                  for (j in to.replace){
                    patient.info[[i]]$sequence[per.i + j] <- round(runif(1, -0.5, R + 0.5))
                  }
                  patient.info[[i]]$X <- individualX(D, patient.info[[i]]$R,
                                                     patient.info[[i]]$sequence)
                }
              }
            }
          }
        }
      }
      while (C < n.vec[l]){
        update.patients.WL <- which(current.status[0:patients.old] == "WL")
        for (i in update.patients.WL){
          waitlist.time <- waitlist.time + 1
          if (sum(current.status != "WL") < n.vec[l]){
            if (is.vector(curr.sequences)){
              patient.info[[i]]$sequence <- curr.sequences
              curr.sequences <- main.sequences[[l]]
            } else {
              rand.seq <- ceiling(runif(1, 0, nrow(curr.sequences)))
              patient.info[[i]]$sequence <- curr.sequences[rand.seq, ]
              curr.sequences <- curr.sequences[-rand.seq, ]
            }
            patient.info[[i]]$X <- individualX(D, R, patient.info[[i]]$sequence)
            patient.info[[i]]$current.status <- "D"
            patient.info[[i]]$curr.period <- 1
            patient.info[[i]]$weeks.status <- 0
            current.status[i]                <- "D"
          } else {
            patient.info[[i]]$weeks.status <- patient.info[[i]]$weeks.status + 1
          }
        }

        rec.week     <- rpois(1, lambda)
        patients.new <- min(patients.old + rec.week, n.vec[L])
        if (patients.new != patients.old){
          for (i in 1:(patients.new - patients.old)){
            patient.info[[patients.old + i]] <- list()
            patient.info[[patients.old + i]]$week.rec <- week
            patient.info[[patients.old + i]]$R <- R
            if (trial.type == "1A" | trial.type == "1B"){
              patient.info[[patients.old + i]]$current.status <- "D"
              patient.info[[patients.old + i]]$curr.period <- 1
            } else if (trial.type == "2"){
              if (patients.old + i > n.vec[l]){
                patient.info[[patients.old + i]]$current.status <- "WL"
              } else {
                patient.info[[patients.old + i]]$current.status <- "D"
                patient.info[[patients.old + i]]$curr.period <- 1
              }
            }
            if (patient.info[[patients.old + i]]$current.status != "WL"){
              if (is.vector(curr.sequences)){
                patient.info[[patients.old + i]]$sequence <- curr.sequences
                curr.sequences <- main.sequences[[l]]
              } else {
                rand.seq <- ceiling(runif(1, 0, nrow(curr.sequences)))
                patient.info[[patients.old + i]]$sequence <- curr.sequences[rand.seq, ]
                curr.sequences <- curr.sequences[-rand.seq, ]
              }
              patient.info[[patients.old + i]]$X <- individualX(D, R, patient.info[[patients.old + i]]$sequence)
            }
            patient.info[[patients.old + i]]$weeks.status   <- 0
            patient.info[[patients.old + i]]$s              <- rnorm(n = 1, mean = 0, sd = sigma.b)
            current.status <- c(current.status, patient.info[[patients.old + i]]$current.status)
          }
        }
        update.patients.C  <- which(current.status[0:patients.old] == "C")
        update.patients.D  <- which(current.status[0:patients.old] == "D")
        update.patients.WO <- which(current.status[0:patients.old] == "WO")
        for (i in update.patients.C){
          patient.info[[i]]$weeks.status   <- patient.info[[i]]$weeks.status + 1
        }
        for (i in update.patients.D){
          patient.info[[i]]$weeks.status <- patient.info[[i]]$weeks.status + 1
          if (patient.info[[i]]$weeks.status == p){
            patient.info[[i]]$Y <- c(patient.info[[i]]$Y,
                                     rnorm(n = 1, mean = patient.info[[i]]$X[patient.info[[i]]$curr.period, ]%*%beta + patient.info[[i]]$s,
                                           sd = sigma.e))
            observations <- observations + 1
            if (patient.info[[i]]$curr.period == length(patient.info[[i]]$sequence)){
              patient.info[[i]]$current.status <- "C"
              patient.info[[i]]$weeks.status   <- 0
              current.status[i]                <- "C"
            } else {
              patient.info[[i]]$current.status <- "WO"
              patient.info[[i]]$weeks.status   <- 0
              current.status[i]                <- "WO"
            }
          }
        }
        for (i in update.patients.WO){
          patient.info[[i]]$weeks.status <- patient.info[[i]]$weeks.status + 1
          if (patient.info[[i]]$weeks.status == w){
            patient.info[[i]]$current.status <- "D"
            patient.info[[i]]$weeks.status <- 0
            patient.info[[i]]$curr.period <- patient.info[[i]]$curr.period + 1
            current.status[i]                <- "D"
          }
        }
        C            <- length(which(current.status == "C"))
        week         <- week + 1
        patients.old <- patients.new
        print(patients.old)
      }

      analysis.df <- analysisData(patient.info, all.data, n.vec[l], patients.old)
      interim.analysis <- lmer(response ~ period + treatment + (1|subject),
                               data = analysis.df, REML = REML)
      beta.hat[l, ] <- fixef(interim.analysis)
      tau.hat       <- beta.hat[l, (D + 1):(2*(D - 1) + 1)]
      I.hat[l, ]    <- 1/diag(vcov(interim.analysis))[(D + 1):(2*(D - 1) + 1)]
      Z.hat         <- tau.hat*sqrt(I.hat[l, ])
      new.drugs.rem <- NULL
      for (i in drugs.rem[-1]){
        if (Z.hat[i] <= f[l]){
          decision[i]   <- 0
          stage.dec[i]  <- l
        } else if (Z.hat[i] > e[l]){
          decision[i]   <- 1
          stage.dec[i]  <- l
        } else if ((Z.hat[i] > f[l]) && (Z.hat[i] <= e[l])){
          new.drugs.rem <- c(new.drugs.rem, i)
        }
      }
      drugs.rem <- c(0, new.drugs.rem)
      R <- length(drugs.rem)
      if (R == 1){
        break
      }
      if (l < L){
        main.sequences[[l + 1]] <- sequences(R, drugs.rem, sequence.type)
        curr.sequences <- main.sequences[[l + 1]]
      }
    }

    patients.active <- length(which(current.status[0:patients.new] != "WL"))
    output <-c(waitlist.time, waitlist.time/patients.new,
               patients.new, patients.active, any(decision == 1), decision[1] == 1,
               observations, week - 1)
    return(output)
  }

  wrapper <- function(i){
    result.i <- singleTrial(i)
    return(result.i)
  }

  ##### MAIN COMPUTATIONS #####################################################

  set.seed(seed)

  if (summary == TRUE){
    print("Determining average performance...")
  }

  sink("NULL")
  suppressMessages(sfInit(parallel = parallel, cpus = cpus))
  suppressMessages(sfLibrary(mvtnorm))
  suppressMessages(sfLibrary(lme4))
  suppressMessages(sfLibrary(stats))
  suppressMessages(sfLibrary(crossdes))
  sfExport("trial.type", "all.data", "D", "sequence.type", "n.vec",
           "sigma.e", "sigma.b", "f", "e", "mu.0", "tau", "pi",
           "lambda", "w", "p", "REML", "singleTrial", "sequences",
           "individualX", "designMatrices",
           "analysisData")
  results           <- sfLapply(1:replicates, wrapper)
  suppressMessages(sfStop())
  sink()
  av.results        <- numeric(8)
  for (j in 1:replicates){
    av.results      <- av.results + results[[j]]
  }
  results           <- matrix(unlist(results), nrow = replicates, ncol = 8, byrow = TRUE)
  av.results        <- av.results/replicates
  names(av.results) <- c("total.waitlist.time", "average.waitlist.time",
                          "total.patients", "active.patients", "power.any",
                          "power.1", "observations", "trial.length")
  colnames(results) <- names(av.results)

  if (summary == TRUE){
    print("Outputting...")
  }

  output <- list(all.data = all.data, av.results = av.results, cpus = cpus,
                 D = D, e = e, f = f, lambda = lambda, mu.0 = mu.0,
                 n.vec = n.vec, p = p, parallel = parallel, pi = pi, REML = REML,
                 replicates = replicates, results = results, seed = seed,
                 sequence.type = sequence.type, sigma.b = sigma.b, sigma.e = sigma.e,
                 tau = tau, trial.type = trial.type, w = w)
 return(output)
}
