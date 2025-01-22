library(coda)
library(lavaan)

source("SDM Sampler_HMC.R")
source("Tools.R")

options(scipen = 999)

if (file.exists(filename)) {
  state <- readRDS(filename)
  which.reps <- (state$summary.lavaan[1]+1):rep # Resume calculations
} else {
  state <- NULL
  which.reps <- 1:rep
}

if(!exists("N")) {
  N <- 500
}
K <- 30

dtype <- 2 # 1: MVN-distributed data, 2: binary data
J <- 50000 # Number of samples to estimate observed-data likelihood

Nt.M0 <- 1 # Only general factor
Nt.M1 <- length(tau.sim) # Full model

# Model M0 (one general factor)
u.M0 <- matrix(0, nrow = Nt.M0, ncol = K)
u.M0[1,] <- 1 # f1

# Model M1 (one general factor + second factor with three component-specific dimensions)
u.M1 <- matrix(0, nrow = Nt.M1, ncol = K)
u.M1[1, 1:(K / 3)] <- 1 # f3.1
u.M1[2, (K / 3 + 1):(2/3 *K)] <- 1 # f3.2
u.M1[3, (2/3*K + 1):K] <- 1 # f3.3
u.M1[Nt.M1,] <- 1 # f1

# Covariance structure (data is generated under M1)
if (dtype == 1) {
  sig2k.sim <- runif(K, 1, 1.5)
  Sigma.sim <- diag(sig2k.sim)
  for (tt in 1:Nt.M1) {
    Sigma.sim <- Sigma.sim + u.M1[tt,] %*% t(u.M1[tt,]) * tau.sim[tt]
  }
} else if (dtype == 2) {
  Sigma.sim <- diag(K)
  for (tt in 1:Nt.M1) {
    Sigma.sim <- Sigma.sim + u.M1[tt,] %*% t(u.M1[tt,]) * tau.sim[tt]
  }
  # Correlation matrix
  Sigma.sim <- Sigma.sim - diag(K)
  sig2k.sim <- 1 - diag(Sigma.sim)
  Sigma.sig2k <- diag(sig2k.sim)
  Sigma.sim <- Sigma.sim + Sigma.sig2k
}

# Number of chains
nChains <- 2

# Configure TS-IG BCSM
XG.TSIG <- 100
burnin.TSIG <- floor(XG.TSIG / 10)

# Configure HMC BCSM
XG.HMC <- 10000
burnin.HMC <- floor(XG.HMC / 10)
runGibbs <- TRUE
runHMC <- TRUE
pre.est <- "mean" # Use means from TS-IG BCSM as set of starting values for HMC
if (dtype == 1) {
  # M0
  epsilon.M0 <- 0.8
  epsilon.scale.M0 <- 1
  epsilon.sd.M0 <- 0#epsilon.M0 * 0.2
  L.M0 <- 25
  L.sd.M0 <- L.M0 * 0.2
  # M1
  epsilon.M1 <- 0.25
  epsilon.scale.M1 <- 1
  epsilon.sd.M1 <- 0#epsilon.M1 * 0.2
  L.M1 <- 25
  L.sd.M1 <- L.M1 * 0.2
} else if (dtype == 2) {
  # M0
  epsilon.M0 <- 0.2
  epsilon.scale.M0 <- 1
  epsilon.sd.M0 <- 0#epsilon.M0 * 0.02
  L.M0 <- 6
  L.sd.M0 <- 1 #L.M0 * 0.5
  # M1
  epsilon.M1 <- 0.055
  if (tau.sim[1] <= -0.02)
    epsilon.M1 <- 0.04
  epsilon.scale.M1 <- 1
  epsilon.sd.M1 <- 0#epsilon.M1 * 0.02
  L.M1 <- 6
  L.sd.M1 <- 1 #L.M1 * 0.5
}

# Prior specification
if (dtype == 1) {
  pr.tau.lower <- -2
  pr.tau.upper <- 2
} else if (dtype == 2) {
  pr.tau.lower <- -1
  pr.tau.upper <- 1
}
if (dtype == 1) {
  pr.sig2k.upper <- 10
} else if (dtype == 2) {
  pr.sig2k.upper <- 1
}
pr.iinter.mu0 <- 0
pr.iinter.var0 <- 1000
a0 <- 5 # shape (for Gibbs-sampler)
b0 <- 5 # rate (for Gibbs-sampler)

# Thresholds for frequentist fit measure indices 
thr.RMSEA1 <- 0.05
thr.RMSEA2 <- 0.10
thr.SRMR <- 0.08
thr.TLI <- 0.95

# Number of replications
#rep <- 5

## Data structures ##
if (is.null(state)) {
  # Bayes factor
  logBF <- matrix(NA, nrow = rep, ncol = 1) # M1 vs M0
  logBF.ESS <- matrix(NA, nrow = rep, ncol = 1) # M1 vs M0, given corrected ESS
  
  # Posterior model probabilities 
  PMP <- matrix(NA, nrow = rep, ncol = 2) # M1 vs M0, M0 vs M1
  PMP.ESS <- matrix(NA, nrow = rep, ncol = 2) # M1 vs M0, M0 vs M1, given corrected ESS
  
  # Marginal log likelihoods
  marllik.M0 <- matrix(NA, nrow = rep, ncol = 1)
  marllik.M1 <- matrix(NA, nrow = rep, ncol = 1)
  marllik.M0.ESS <- matrix(NA, nrow = rep, ncol = 1) # Given corrected ESS
  marllik.M1.ESS <- matrix(NA, nrow = rep, ncol = 1) # Given corrected ESS
  
  # Posterior means and SDs
  est.M0 <- matrix(NA, nrow = rep, ncol = 2*Nt.M0) # First half are the posterior means of tau, other half the SDs 
  est.M1 <- matrix(NA, nrow = rep, ncol = 2*Nt.M1) # First half are the posterior means of tau, other half the SDs 
  
  # Lavaan
  lavaan.est.M0 <- matrix(NA, nrow = rep, ncol = 2*Nt.M0) # First half are the estimates of tau, other half the SEs 
  lavaan.est.M1 <- matrix(NA, nrow = rep, ncol = 2*Nt.M1) # First half are the estimates of tau, other half the SEs 
  lavaan.FM.M0 <- matrix(NA, nrow = rep, ncol = 7) # RMSEA, RMSEA.Scaled, SRMR, CFI, CFI.Scaled, TLI, TLI.Scaled
  lavaan.FM.M1 <- matrix(NA, nrow = rep, ncol = 7) # RMSEA, RMSEA.Scaled, SRMR, CFI, CFI.Scaled, TLI, TLI.Scaled
  lavaan.FM.COMP <- matrix(NA, nrow = rep, ncol = 4) # CFI, CFI.Scaled, TLI, TLI.scaled
} else {
  # Bayes factor
  logBF <- state$results.BCSM$logBF # M1 vs M0
  logBF.ESS <- state$results.BCSM$logBF.ESS # M1 vs M0, given corrected ESS
  
  # Posterior model probabilities 
  PMP <- state$results.BCSM$PMP # M1 vs M0, M0 vs M1
  PMP.ESS <- state$results.BCSM$PMP.ESS # M1 vs M0, M0 vs M1, given corrected ESS
  
  # Marginal log likelihoods
  marllik.M0 <- state$results.BCSM$marllik.M0
  marllik.M1 <- state$results.BCSM$marllik.M1
  marllik.M0.ESS <- state$results.BCSM$marllik.M0.ESS # Given corrected ESS
  marllik.M1.ESS <- state$results.BCSM$marllik.M1.ESS # Given corrected ESS
  
  # Posterior means and SDs
  est.M0 <- state$results.BCSM$est.M0 # First half are the posterior means of tau, other half the SDs 
  est.M1 <- state$results.BCSM$est.M1 # First half are the posterior means of tau, other half the SDs 
  
  # Lavaan
  lavaan.est.M0 <- state$results.lavaan$lavaan.est.M0  # First half are the estimates of tau, other half the SEs 
  lavaan.est.M1 <- state$results.lavaan$lavaan.est.M1 # First half are the estimates of tau, other half the SEs 
  lavaan.FM.M0 <- state$results.lavaan$lavaan.FM.M0 # RMSEA, RMSEA.Scaled, SRMR, CFI, CFI.Scaled, TLI, TLI.Scaled
  lavaan.FM.M1 <- state$results.lavaan$lavaan.FM.M1 # RMSEA, RMSEA.Scaled, SRMR, CFI, CFI.Scaled, TLI, TLI.Scaled
  lavaan.FM.COMP <- state$results.lavaan$lavaan.FM.COMP # CFI, CFI.Scaled, TLI, TLI.scaled
}


# Summary
summary.BCSM <- NULL
summary.lavaan <- NULL

start_time <- Sys.time()
for (rr in which.reps) {
  # Simulate item intercepts
  iinter.sim <- rnorm(K, 0, 1)
  if (dtype == 2) {
    iinter.sim <- iinter.sim - mean(iinter.sim)
  }
  
  # Mean structure
  mu.sim <- iinter.sim
  
  # Generate data
  done <- FALSE
  while (!done) {
    Z <- MASS::mvrnorm(n = N, mu = mu.sim, Sigma = Sigma.sim, empirical = FALSE)
    if (dtype == 1) {
      Y <- Z
    } else if (dtype == 2) {
      Y <- matrix(0, nrow = N, ncol = K)
      Y[which(Z > 0)] <- 1
    }
    done <- TRUE
    for (kk in 1:K) {
      if (all(Y[, kk] == 0) | all(Y[, kk] == 1)) { 
        done <- FALSE
      } 
    }
  }
  colnames(Y) <- colnames(Z) <- paste0("x", 1:ncol(Y))
  
  # Lavaan model specification
  m.M0 <- paste("f1 =~", paste(colnames(Y), collapse = "+1*"))
  m.M1 <- paste("f1 =~", paste(colnames(Y), collapse = "+1*"),
                "\nf3.1 =~", paste(colnames(Y[, as.logical(u.M1[1, ])]), collapse = "+1*"), 
                "\nf3.2 =~", paste(colnames(Y[, as.logical(u.M1[2, ])]), collapse = "+1*"),
                "\nf3.3 =~", paste(colnames(Y[, as.logical(u.M1[3, ])]), collapse = "+1*"))

  # Fit M0
  cat("### Fitting M0 ###", "\n")
  Samples.TSIG.M0 <- numeric(0)
  Samples.HMC.M0 <- numeric(0)
  for (cc in 1:nChains) {
    cat("*** Chain:", cc, "\n")
    out.M0 <- SamplerSDM_HMC(Y = Y, u = u.M0, dtype = dtype, epsilon = epsilon.M0, epsilon.scale = epsilon.scale.M0, epsilon.sd = epsilon.sd.M0, enforceLB.tau = TRUE
                            , L = L.M0, L.sd = L.sd.M0, save.Z = FALSE, save.Z.Gibbs = FALSE, save.DistPar = FALSE
                            , pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, a0 = a0, b0 = b0
                            , pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0
                            , XG = XG.HMC, XG.Gibbs = XG.TSIG, Burnin.Gibbs = burnin.TSIG, runHMC = runHMC, runGibbs = runGibbs, Gibbs.trmethod = "ICDF", pre.est = pre.est, Gibbs.Z.start = Z
                            , ii.print = 5000, ii.print.Gibbs = 5000, verbose = FALSE)
    Samples.TSIG.M0 <- rbind(Samples.TSIG.M0, cbind(out.M0$Samples.tau.Gibbs[(burnin.TSIG + 1):XG.TSIG], out.M0$Samples.sig2k.Gibbs[(burnin.TSIG + 1):XG.TSIG,], out.M0$Samples.iinter.Gibbs[(burnin.TSIG + 1):XG.TSIG,]))
    Samples.HMC.M0 <- rbind(Samples.HMC.M0, cbind(out.M0$Samples.tau[(burnin.HMC + 1):XG.HMC], out.M0$Samples.sig2k[(burnin.HMC + 1):XG.HMC,], out.M0$Samples.iinter[(burnin.HMC + 1):XG.HMC,]))
  }
  colnames(Samples.TSIG.M0) <- colnames(Samples.HMC.M0) <- c("tau", paste("sig2_", 1:K, sep = ""), paste("iinter_", 1:K, sep = ""))
  est.M0[rr, 1:Nt.M0] <- colMeans(matrix(Samples.HMC.M0[,1:Nt.M0], nrow = nrow(Samples.HMC.M0), ncol = Nt.M0))
  est.M0[rr, (Nt.M0+1):(2*Nt.M0)] <- colSDs(matrix(Samples.HMC.M0[,1:Nt.M0], nrow = nrow(Samples.HMC.M0), ncol = Nt.M0))

  # Fit M1
  cat("### Fitting M1 ###", "\n")
  Samples.TSIG.M1 <- numeric(0)
  Samples.HMC.M1 <- numeric(0)
  for (cc in 1:nChains) {
    cat("*** Chain:", cc, "\n")
    out.M1 <- SamplerSDM_HMC(Y = Y, u = u.M1, dtype = dtype, epsilon = epsilon.M1, epsilon.scale = epsilon.scale.M1, epsilon.sd = epsilon.sd.M1, enforceLB.tau = TRUE
                            , L = L.M1, L.sd = L.sd.M1, save.Z = FALSE, save.Z.Gibbs = FALSE, save.DistPar = FALSE
                            , pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, a0 = a0, b0 = b0
                            , pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0
                            , XG = XG.HMC, XG.Gibbs = XG.TSIG, Burnin.Gibbs = burnin.TSIG, runHMC = runHMC, runGibbs = runGibbs, Gibbs.trmethod = "ICDF", Gibbs.tau.start = tau.sim, Gibbs.Z.start = Z, pre.est = pre.est
                            , ii.print = 5000, ii.print.Gibbs = 5000, verbose = FALSE)
    Samples.TSIG.M1 <- rbind(Samples.TSIG.M1, cbind(out.M1$Samples.tau.Gibbs[(burnin.TSIG + 1):XG.TSIG,], out.M1$Samples.sig2k.Gibbs[(burnin.TSIG + 1):XG.TSIG,], out.M1$Samples.iinter.Gibbs[(burnin.TSIG + 1):XG.TSIG,]))
    Samples.HMC.M1 <- rbind(Samples.HMC.M1, cbind(out.M1$Samples.tau[(burnin.HMC + 1):XG.HMC,], out.M1$Samples.sig2k[(burnin.HMC + 1):XG.HMC,], out.M1$Samples.iinter[(burnin.HMC + 1):XG.HMC,]))
  }
  colnames(Samples.TSIG.M1) <- colnames(Samples.HMC.M1) <- c(paste("tau_", 1:Nt.M1, sep = ""), paste("sig2_", 1:K, sep = ""), paste("iinter_", 1:K, sep = ""))
  est.M1[rr, 1:Nt.M1] <- colMeans(matrix(Samples.HMC.M1[,1:Nt.M1], nrow = nrow(Samples.HMC.M1), ncol = Nt.M1))
  est.M1[rr, (Nt.M1+1):(2*Nt.M1)] <- colSDs(matrix(Samples.HMC.M1[,1:Nt.M1], nrow = nrow(Samples.HMC.M1), ncol = Nt.M1))

  # Lavaan
  fit.M0 <- lavaan::cfa(m.M0, data = Y, ordered = TRUE, orthogonal = TRUE, parameterization = "delta")
  lavaan.est.M0[rr, 1:Nt.M0] <- diag(lavaan::lavInspect(fit.M0, what = "est")$psi)
  lavaan.est.M0[rr, (Nt.M0+1):(2*Nt.M0)] <- diag(lavaan::lavInspect(fit.M0, what = "se")$psi)
  lavaan.FM.M0[rr, ] <- lavaan::fitmeasures(fit.M0, c("rmsea", "rmsea.scaled", "srmr", "cfi", "cfi.scaled", "tli", "tli.scaled")) 
  colnames(lavaan.FM.M0) <- c("rmsea", "rmsea.scaled", "srmr", "cfi", "cfi.scaled", "tli", "tli.scaled")
  
  fit.M1 <- lavaan::cfa(m.M1, data = Y, ordered = TRUE, orthogonal = TRUE, parameterization = "delta")
  lavaan.est.M1[rr, 1:Nt.M1] <- diag(lavaan::lavInspect(fit.M1, what = "est")$psi)
  lavaan.est.M1[rr, (Nt.M1+1):(2*Nt.M1)] <- diag(lavaan::lavInspect(fit.M1, what = "se")$psi)
  lavaan.FM.M1[rr, ] <- lavaan::fitmeasures(fit.M1, c("rmsea", "rmsea.scaled", "srmr", "cfi", "cfi.scaled", "tli", "tli.scaled")) 
  colnames(lavaan.FM.M1) <- c("rmsea", "rmsea.scaled", "srmr", "cfi", "cfi.scaled", "tli", "tli.scaled")
  lavaan.FM.COMP[rr, ] <- lavaan::fitmeasures(fit.M1, c("cfi", "cfi.scaled", "tli", "tli.scaled"), baseline.model = fit.M0) 
  colnames(lavaan.FM.COMP) <- c("cfi", "cfi.scaled", "tli", "tli.scaled")
  # rbind(lavaan.FM.M0[rr, ], lavaan.FM.M1[rr, ])
  # rbind(c(lavaan.est.M0[rr, 1:Nt.M0], NA, NA), lavaan.est.M1[rr, 1:Nt.M1])
  # lavaan::fitmeasures(fit.M1, c("cfi", "cfi.scaled", "tli", "tli.scaled"), baseline.model = fit.M0) 
  
    
  #### Bayes factor (approximated through BICs) ####
  cat("### Computing marginal likelihood estimates ###", "\n")

  # Estimate marginal likelihood under M0
  cat("*** M0", "\n")

  # Posterior samples
  iinter.HMC <- Samples.HMC.M0[, (Nt.M0 + K + 1):ncol(Samples.HMC.M0)]
  sig2k.HMC <- Samples.HMC.M0[, (Nt.M0 + 1):(Nt.M0 + K)]
  tau.HMC <- as.matrix(Samples.HMC.M0[, 1:(Nt.M0)], nrow = nrow(Samples.HMC.M0), ncol = Nt.M0)
  
  # Posterior mean estimates
  iinter.HMC.est <- colMeans(iinter.HMC)
  sig2k.HMC.est <- colMeans(sig2k.HMC)
  tau.HMC.est <- colMeans(tau.HMC)

  # BIC
  Sigma.M0 <- diag(sig2k.HMC.est)
  for (tt in 1:Nt.M0) {
    Sigma.M0 <- Sigma.M0 + u.M0[tt,] %*% t(u.M0[tt,]) * tau.HMC.est[tt]
  }
  Np.M0 <- Nt.M0 + K # Nt covariance parameter(s), K item intercepts
  NN.M0 <- N * K # Effective sample size, safe choice
  NN.M0.ESS <- sum(diag(Sigma.M0)) / sum(Sigma.M0) * NN.M0 # Effective sample size, corrected for correlations
  Z <-  MASS::mvrnorm(N, iinter.HMC.est, Sigma.M0, empirical = TRUE)
  B11.list <- vector(mode = "list", length = K)
  B12.list <- vector(mode = "list", length = K)
  B21.list <- vector(mode = "list", length = K)
  B22.inv.list <- vector(mode = "list", length = K)
  mu.part.list <- vector(mode = "list", length = K)
  sig2.list <- vector(mode = "list", length = K)
  for (kk in 1:K) {
    condition.on <- (1:K)[-kk]
    B11.list[[kk]] <- Sigma.M0[kk, kk]
    B12.list[[kk]] <- Sigma.M0[kk, condition.on]
    B21.list[[kk]] <- Sigma.M0[condition.on, kk]
    B22.inv.list[[kk]] <- solve(Sigma.M0[condition.on, condition.on])
    mu.part.list[[kk]] <- B12.list[[kk]] %*% B22.inv.list[[kk]]
    sig2.list[[kk]] <- B11.list[[kk]] - B12.list[[kk]] %*% B22.inv.list[[kk]] %*% B21.list[[kk]]
  }
  ll <- numeric(N)
  for (jj in 1:J) {
    for (kk in 1:K) {
      condition.on <- (1:K)[-kk] # condition on Z[-kk]
      B11 <- B11.list[[kk]]
      B12 <- B12.list[[kk]]
      B21 <- B21.list[[kk]]
      B22.inv <- B22.inv.list[[kk]]
      mu <- iinter.HMC.est[kk] + mu.part.list[[kk]] %*% (matrix(t(Z[, condition.on]), nrow = K - 1, ncol = N) - matrix(iinter.HMC.est[condition.on], nrow = K - 1, ncol = N, byrow = FALSE))
      sig2 <-sig2.list[[kk]] 
      Z[Y[, kk] == 0, kk] <- extraDistr::rtnorm(n = length(mu[1, Y[, kk] == 0]), mean = mu[1, Y[, kk] == 0], sd = sqrt(sig2[1, 1]), a = -Inf, b = 0)
      Z[Y[, kk] == 1, kk] <- extraDistr::rtnorm(n = length(mu[1, Y[, kk] == 1]), mean = mu[1, Y[, kk] == 1], sd = sqrt(sig2[1, 1]), a = 0, b = Inf)
      Z[, kk] <- Z[, kk] / sd(Z[, kk]) # Delta parameterization
    }
    ll <- ll + Rfast::dmvnorm(Z, mu = iinter.HMC.est, sigma = Sigma.M0, logged = FALSE) # Estimate likelihood of observed data
  }
  ll <- ll / J
  ll <- log(ll)
  ll <- sum(ll)
  BIC.M0 <- Np.M0 * log(NN.M0) - 2 * ll
  BIC.M0.ESS <- Np.M0 * log(NN.M0.ESS) - 2 * ll
  mar.M0 <- -BIC.M0 / 2
  mar.M0.ESS <- -BIC.M0.ESS / 2


  # Estimate marginal likelihood under M1
  cat("*** M1", "\n")

  # Posterior samples
  iinter.HMC <- Samples.HMC.M1[, (Nt.M1 + K + 1):ncol(Samples.HMC.M1)]
  sig2k.HMC <- Samples.HMC.M1[, (Nt.M1 + 1):(Nt.M1 + K)]
  tau.HMC <- as.matrix(Samples.HMC.M1[, 1:(Nt.M1)], nrow = nrow(Samples.HMC.M1), ncol = Nt.M1)
  
  # Posterior mean estimates
  iinter.HMC.est <- colMeans(iinter.HMC)
  sig2k.HMC.est <- colMeans(sig2k.HMC)
  tau.HMC.est <- colMeans(tau.HMC)

  # BIC
  Sigma.M1 <- diag(sig2k.HMC.est)
  for (tt in 1:Nt.M1) {
    Sigma.M1 <- Sigma.M1 + u.M1[tt,] %*% t(u.M1[tt,]) * tau.HMC.est[tt]
  }
  Np.M1 <- Nt.M1 + K # Nt covariance parameter(s), K item intercepts
  NN.M1 <- N * K # Effective sample size, safe choice
  NN.M1.ESS <- sum(diag(Sigma.M1)) / sum(Sigma.M1) * NN.M1 # Effective sample size, corrected for correlations
  Z <-  MASS::mvrnorm(N, iinter.HMC.est, Sigma.M1, empirical = TRUE)
  B11.list <- vector(mode = "list", length = K)
  B12.list <- vector(mode = "list", length = K)
  B21.list <- vector(mode = "list", length = K)
  B22.inv.list <- vector(mode = "list", length = K)
  mu.part.list <- vector(mode = "list", length = K)
  sig2.list <- vector(mode = "list", length = K)
  for (kk in 1:K) {
    condition.on <- (1:K)[-kk]
    B11.list[[kk]] <- Sigma.M1[kk, kk]
    B12.list[[kk]] <- Sigma.M1[kk, condition.on]
    B21.list[[kk]] <- Sigma.M1[condition.on, kk]
    B22.inv.list[[kk]] <- solve(Sigma.M1[condition.on, condition.on])
    mu.part.list[[kk]] <- B12.list[[kk]] %*% B22.inv.list[[kk]]
    sig2.list[[kk]] <- B11.list[[kk]] - B12.list[[kk]] %*% B22.inv.list[[kk]] %*% B21.list[[kk]]
  }
  ll <- numeric(N)
  for (jj in 1:J) {
    for (kk in 1:K) {
      condition.on <- (1:K)[-kk] # condition on Z[-kk]
      B11 <- B11.list[[kk]]
      B12 <- B12.list[[kk]]
      B21 <- B21.list[[kk]]
      B22.inv <- B22.inv.list[[kk]]
      mu <- iinter.HMC.est[kk] + mu.part.list[[kk]] %*% (matrix(t(Z[, condition.on]), nrow = K - 1, ncol = N) - matrix(iinter.HMC.est[condition.on], nrow = K - 1, ncol = N, byrow = FALSE))
      sig2 <-sig2.list[[kk]] 
      Z[Y[, kk] == 0, kk] <- extraDistr::rtnorm(n = length(mu[1, Y[, kk] == 0]), mean = mu[1, Y[, kk] == 0], sd = sqrt(sig2[1, 1]), a = -Inf, b = 0)
      Z[Y[, kk] == 1, kk] <- extraDistr::rtnorm(n = length(mu[1, Y[, kk] == 1]), mean = mu[1, Y[, kk] == 1], sd = sqrt(sig2[1, 1]), a = 0, b = Inf)
      Z[, kk] <- Z[, kk] / sd(Z[, kk]) # Delta parameterization
    }
    ll <- ll + Rfast::dmvnorm(Z, mu = iinter.HMC.est, sigma = Sigma.M1, logged = FALSE) # Estimate likelihood of observed data
  }
  ll <- ll / J
  ll <- log(ll)
  ll <- sum(ll)
  BIC.M1 <- Np.M1 * log(NN.M1) - 2 * ll
  BIC.M1.ESS <- Np.M1 * log(NN.M1.ESS) - 2 * ll
  mar.M1 <- -BIC.M1 / 2
  mar.M1.ESS <- -BIC.M1.ESS / 2

  marllik.M0[rr,] <- mar.M0
  marllik.M1[rr,] <- mar.M1
  marllik.M0.ESS[rr,] <- mar.M0.ESS
  marllik.M1.ESS[rr,] <- mar.M1.ESS
  logBF[rr,] <- mar.M1 - mar.M0
  logBF.ESS[rr,] <- mar.M1.ESS - mar.M0.ESS
  PMP[rr,] <- bridgesampling::post_prob(mar.M1, mar.M0)
  PMP.ESS[rr,] <- bridgesampling::post_prob(mar.M1.ESS, mar.M0.ESS)

  results.BCSM <- list(est.M0 = est.M0, est.M1 = est.M1, marllik.M0 = marllik.M0, marllik.M1 = marllik.M1, marllik.M0.ESS = marllik.M0.ESS, marllik.M1.ESS = marllik.M1.ESS, logBF = logBF, logBF.ESS = logBF.ESS, PMP = PMP, PMP.ESS = PMP.ESS)
  results.lavaan <- list(lavaan.est.M0 = lavaan.est.M0, lavaan.est.M1 = lavaan.est.M1, lavaan.FM.M0 = lavaan.FM.M0, lavaan.FM.M1 = lavaan.FM.M1, lavaan.FM.COMP = lavaan.FM.COMP)

  summary.BCSM <- c(rr, tau.sim[1], mean(logBF.ESS[1:rr, ]), mean(logBF[1:rr, ]), sd(logBF.ESS[1:rr, ]), sd(logBF[1:rr, ]), length(which(PMP.ESS[1:rr, 1] > 0.5))/rr*100, length(which(PMP.ESS[1:rr, 2] > 0.5))/rr*100, length(which(PMP[1:rr, 1] > 0.5))/rr*100, length(which(PMP[1:rr, 2] > 0.5))/rr*100)
  summary.lavaan <- c(rr, tau.sim[1], colMeans(FM.threshold(lavaan.FM.M0[1:rr, 1], lavaan.FM.M1[1:rr, 1], thr.RMSEA1)$select)*100 # RMSEA thr = .05
                      , colMeans(FM.threshold(lavaan.FM.M0[1:rr, 1], lavaan.FM.M1[1:rr, 1], thr.RMSEA2)$select)*100 # RMSEA thr = .10
                      , colMeans(FM.delta(lavaan.FM.M0[1:rr, 1], lavaan.FM.M1[1:rr, 1])$select)*100 # RMSEA delta
                      , colMeans(FM.threshold(lavaan.FM.M0[1:rr, 3], lavaan.FM.M1[1:rr, 3], thr.SRMR)$select)*100 # SRMR thr = .08
                      , colMeans(FM.delta(lavaan.FM.M0[1:rr, 3], lavaan.FM.M1[1:rr, 3])$select)*100 # SRMR delta 
                      , colMeans(FM.threshold(1-lavaan.FM.M0[1:rr, 6], 1-lavaan.FM.M1[1:rr, 6], 1-thr.TLI)$select)*100 # TLI thr = .95
                      , colMeans(FM.delta(1-lavaan.FM.M0[1:rr, 6], 1-lavaan.FM.M1[1:rr, 6])$select)*100 # TLI delta
                      )
 names(summary.BCSM) <- c("Replication", "tau.sim", "mean(logBF.ESS)", "mean(logBF)", "sd(logBF.ESS)", "sd(logBF)", "Select.ESS% M1", "Select.ESS% M0", "Select% M1", "Select% M0")
 names(summary.lavaan) <- c("Replication", "tau.sim", "RMSEA.thr1% M0", "RMSEA.thr1% M1", "RMSEA.thr2% M0", "RMSEA.thr2% M1", "RMSEA.delta% M0", "RMSEA.delta% M1"
                            , "SRMR.thr% M0", "SRMR.thr% M1", "SRMR.delta% M0", "SRMR.delta% M1"
                            , "TLI.thr% M0", "TLI.thr% M1", "TLI.delta% M0", "TLI.delta% M1")
 
 print(round(summary.BCSM, 2))
 print(round(summary.lavaan, 2))

 saveRDS(list(tau.sim = tau.sim, results.BCSM = results.BCSM, results.lavaan = results.lavaan, summary.BCSM = summary.BCSM, summary.lavaan = summary.lavaan), file = filename)
 
 # cat("Replication", rr, "complete", "|", "logBF:", logBF[rr,], "|", "PMP:", PMP[rr,], "|", "logBF.ESS:", logBF.ESS[rr,], "|", "PMP.ESS:", PMP.ESS[rr,], "\n")
}
end_time <- Sys.time()
end_time - start_time


