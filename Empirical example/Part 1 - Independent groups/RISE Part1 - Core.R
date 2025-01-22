library(coda)

source("SDM Sampler_HMC.R")
source("../RISE Sample.R") # update u, create Y6, Y7, Y8 and Y9
options(scipen = 999)

N <- nrow(Y6)
K <- ncol(u)
Nt <- nrow(u)

dtype <- 2 # 1: MVN-distributed data, 2: binary data
J <- 500000 # Number of samples to estimate observed-data likelihood

# Number of chains
nChains <- 2

# Configure TS-IG BCSM
XG.TSIG <- 100
burnin.TSIG <- floor(XG.TSIG / 10)

# Configure HMC BCSM
XG.HMC <- 30000
burnin.HMC <- floor(XG.HMC / 10)
runGibbs <- TRUE
runHMC <- TRUE
pre.est <- "mean" # Use random sample from TS-IG BCSM as set of starting values for HMC

# Prior specification
pr.tau.lower <- -1
pr.tau.upper <- 1
pr.iinter.mu0 <- 0
pr.iinter.var0 <- 1000
a0 <- 15 # shape (for Gibbs sampler)
b0 <- 15 # rate (for Gibbs sampler)

# Number of replications (one per group)
rep <- 4

## Data structures ##
# Marginal log likelihoods
marllik <- state$marllik #matrix(NA, nrow = rep, ncol = 1)
marllik.ESS <- state$marllik.ESS #matrix(NA, nrow = rep, ncol = 1) # Given corrected ESS

# Posterior means and SDs
est <- state$est# matrix(NA, nrow = rep, ncol = 2*Nt) # First half are the posterior means of tau, other half the SDs 

start_time <- Sys.time()
for (rr in which.reps) {
  # Select data
  if (rr == 1) {
    Y <- Y6
  } else if (rr == 2) {
    Y <- Y7
  } else if (rr == 3) {
    Y <- Y8
  } else if (rr == 4) {
    Y <- Y9
  }

  # Fit BCSM
  cat("### Fitting Dependence Model ###", "\n")
  Samples.HMC <- numeric(0)
  for (cc in 1:nChains) {
    cat("*** Chain:", cc, "\n")
    out <- SamplerSDM_HMC(Y = Y, u = u, dtype = dtype, ftype = ftype, epsilon = epsilon, epsilon.scale = epsilon.scale, epsilon.sd = epsilon.sd, enforceLB.tau = TRUE
                            , L = L, L.sd = L.sd, save.Z = FALSE, save.Z.Gibbs = FALSE, save.DistPar = FALSE
                            , pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = NULL, a0 = a0, b0 = b0
                            , pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0
                            , XG = XG.HMC, XG.Gibbs = XG.TSIG, Burnin.Gibbs = burnin.TSIG, runHMC = runHMC, runGibbs = runGibbs, Gibbs.trmethod = "RS", pre.est = pre.est
                            , ii.print = 100, ii.print.Gibbs = 1000, verbose = FALSE)
    Samples.HMC <- rbind(Samples.HMC, cbind(out$Samples.tau[(burnin.HMC + 1):XG.HMC, ], out$Samples.sig2k[(burnin.HMC + 1):XG.HMC,], out$Samples.iinter[(burnin.HMC + 1):XG.HMC,]))
  }
  if (Nt == 1) {
    colnames(Samples.HMC) <- c("tau", paste("sig2_", 1:K, sep = ""), paste("iinter_", 1:K, sep = ""))
  } else {
    colnames(Samples.HMC) <- c(paste("tau_", 1:Nt, sep = ""), paste("sig2_", 1:K, sep = ""), paste("iinter_", 1:K, sep = ""))
  }
  est[rr, 1:Nt] <- colMeans(matrix(Samples.HMC[,1:Nt], nrow = nrow(Samples.HMC), ncol = Nt))
  est[rr, (Nt+1):(2*Nt)] <- colSDs(matrix(Samples.HMC[,1:Nt], nrow = nrow(Samples.HMC), ncol = Nt))
    
  #### Bayes factor (approximated through BICs) ####
  cat("### Computing marginal likelihood estimates ###", "\n")

  # Estimate marginal likelihood
  
  # Posterior samples
  iinter.HMC <- Samples.HMC[, (Nt + K + 1):ncol(Samples.HMC)]
  sig2k.HMC <- Samples.HMC[, (Nt + 1):(Nt + K)]
  tau.HMC <- as.matrix(Samples.HMC[, 1:(Nt)], nrow = nrow(Samples.HMC), ncol = Nt)
  
  # Posterior mean estimates
  iinter.HMC.est <- colMeans(iinter.HMC)
  sig2k.HMC.est <- colMeans(sig2k.HMC)
  tau.HMC.est <- colMeans(tau.HMC)

  # BIC
  Sigma <- diag(sig2k.HMC.est)
  for (tt in 1:Nt) {
    Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * tau.HMC.est[tt]
  }
  Np <- Nt + K # Nt covariance parameter(s), K item intercepts
  NN <- N * K # Effective sample size, safe choice
  NN.ESS <- sum(diag(Sigma)) / sum(Sigma) * NN # Effective sample size, corrected for correlations
  Z <- MASS::mvrnorm(N, iinter.HMC.est, Sigma, empirical = TRUE)
  B11.list <- vector(mode = "list", length = K)
  B12.list <- vector(mode = "list", length = K)
  B21.list <- vector(mode = "list", length = K)
  B22.inv.list <- vector(mode = "list", length = K)
  mu.part.list <- vector(mode = "list", length = K)
  sig2.list <- vector(mode = "list", length = K)
  for (kk in 1:K) {
    condition.on <- (1:K)[-kk]
    B11.list[[kk]] <- Sigma[kk, kk]
    B12.list[[kk]] <- Sigma[kk, condition.on]
    B21.list[[kk]] <- Sigma[condition.on, kk]
    B22.inv.list[[kk]] <- solve(Sigma[condition.on, condition.on])
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
    ll <- ll + Rfast::dmvnorm(Z, mu = iinter.HMC.est, sigma = Sigma, logged = FALSE) # Estimate likelihood of observed data
    if (jj %% 10000 == 0) {
      ll.tmp <- ll / jj
      ll.tmp <- log(ll.tmp)
      ll.tmp <- sum(ll.tmp)
      print(ll.tmp - (1/2)*Np*log(NN))
    }
  }
  ll <- ll / J
  ll <- log(ll)
  ll <- sum(ll)
  BIC <- Np * log(NN) - 2 * ll
  BIC.ESS <- Np * log(NN.ESS) - 2 * ll
  marllik[rr,] <- -BIC / 2
  marllik.ESS[rr,] <- -BIC.ESS / 2

  saveRDS(list(est = est, marllik = marllik, marllik.ESS = marllik.ESS), file = filename)
  cat("Replication", rr, "complete", "|", "est:", round(est[rr,],3), "|", "marllik:", marllik[rr,], "|", "marllik.ESS:", marllik.ESS[rr,],  "\n")
}
end_time <- Sys.time()
end_time - start_time


