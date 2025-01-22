library(coda)

source("SDM Sampler_HMC_VS.R")
source("../RISE Sample.R") # update u, create Y6, Y7, Y8 and Y9
options(scipen = 999)

N <- nrow(Y6)
K <- ncol(u)
Nt <- nrow(u)
Ng <- 4 # Number of groups

dtype <- 2 # 1: MVN-distributed data, 2: binary data
J <- 500000 # Number of samples to estimate observed-data likelihood

# Number of chains
nChains <- 2

# Configure TS-IG BCSM
XG.TSIG <- 20
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
pr.sig2k.upper <- 1
pr.iinter.mu0 <- 0
pr.iinter.var0 <- 1000
a0 <- 15 # shape (for Gibbs sampler)
b0 <- 15 # rate (for Gibbs sampler)


## Data structures ##
# Marginal log likelihoods
marllik <- numeric(Ng)
marllik.ESS <- numeric(Ng) # Given corrected ESS

gdata <- list()
for (gg in 1:Ng) {
  gdata[[gg]] <- list()
  
  # Posterior samples
  gdata[[gg]]$Samples.tau <- numeric(0)
  gdata[[gg]]$Samples.sig2k <- numeric(0)
  gdata[[gg]]$Samples.ginter <- numeric(0)
  gdata[[gg]]$Samples.tinter <- numeric(0)
  gdata[[gg]]$Samples.iinter <- numeric(0)
  
  # Posterior means and SDs
  # First half are the posterior means of tau, other half the SDs 
  gdata[[gg]]$est.tau <- matrix(NA, nrow = 1, ncol = 2*Nt) 
  gdata[[gg]]$est.sig2k <- matrix(NA, nrow = 1, ncol = 2*K) 
  gdata[[gg]]$est.ginter <- matrix(NA, nrow = 1, ncol = 2) 
  gdata[[gg]]$est.tinter <- matrix(NA, nrow = 1, ncol = 2*Nt) 
  gdata[[gg]]$est.iinter <- matrix(NA, nrow = 1, ncol = 2*K) 
}

Y <- list()
Y[[1]] <- as.matrix(Y6) # Reference group
Y[[2]] <- as.matrix(Y7)
Y[[3]] <- as.matrix(Y8)
Y[[4]] <- as.matrix(Y9)

# Fit BDM
cat("### Fitting Dependence Model ###", "\n")
for (cc in 1:nChains) {
  cat("*** Chain:", cc, "\n")
  out <- SamplerSDM_HMC_VS(Y = Y, u = u, dtype = dtype, ftype = ftype, est.means = TRUE, epsilon = epsilon, epsilon.scale = epsilon.scale, epsilon.sd = epsilon.sd, enforceLB.tau = TRUE
                            , L = L, L.sd = L.sd, save.Z = FALSE, save.Z.Gibbs = FALSE, save.DistPar = FALSE
                            , pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, a0 = a0, b0 = b0
                            , pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0
                            , XG = XG.HMC, XG.Gibbs = XG.TSIG, Burnin.Gibbs = burnin.TSIG, runHMC = runHMC, runGibbs = runGibbs, Gibbs.trmethod = "ICDF"
                            , ii.print = 10, ii.print.Gibbs = 100, verbose = FALSE)
  for (gg in 1:Ng) {
    gdata[[gg]]$Samples.tau <- rbind(gdata[[gg]]$Samples.tau, out$gdata[[gg]]$Samples.tau[burnin.HMC:XG.HMC,])
    gdata[[gg]]$Samples.sig2k <- rbind(gdata[[gg]]$Samples.sig2k, out$gdata[[gg]]$Samples.sig2k[burnin.HMC:XG.HMC,])
    gdata[[gg]]$Samples.ginter <- rbind(gdata[[gg]]$Samples.ginter, out$gdata[[gg]]$Samples.ginter[burnin.HMC:XG.HMC,])
    gdata[[gg]]$Samples.tinter <- rbind(gdata[[gg]]$Samples.tinter, out$gdata[[gg]]$Samples.tinter[burnin.HMC:XG.HMC,])
    gdata[[gg]]$Samples.iinter <- rbind(gdata[[gg]]$Samples.iinter, out$gdata[[gg]]$Samples.iinter[burnin.HMC:XG.HMC,])
  }
}


# Compute posterior means and SDs 
for (gg in 1:Ng) {
  if (Nt == 1) {
    gdata[[gg]]$est.tau[1] <- mean(gdata[[gg]]$Samples.tau)
    gdata[[gg]]$est.tau[2] <- sd(gdata[[gg]]$Samples.tau)
  }
  else {
    gdata[[gg]]$est.tau[1:Nt] <- colMeans(gdata[[gg]]$Samples.tau)
    gdata[[gg]]$est.tau[(Nt+1):(2*Nt)] <- colSDs(gdata[[gg]]$Samples.tau)
  }
  
  gdata[[gg]]$est.sig2k[1:K] <- colMeans(gdata[[gg]]$Samples.sig2k)
  gdata[[gg]]$est.sig2k[(K+1):(2*K)] <- colSDs(gdata[[gg]]$Samples.sig2k)
  
  gdata[[gg]]$est.ginter[1] <- mean(gdata[[gg]]$Samples.ginter)
  gdata[[gg]]$est.ginter[2] <- sd(gdata[[gg]]$Samples.ginter)
  
  if (Nt == 1) {
    gdata[[gg]]$est.tinter[1] <- mean(gdata[[gg]]$Samples.tinter)
    gdata[[gg]]$est.tinter[2] <- sd(gdata[[gg]]$Samples.tinter)
  }
  else {
    gdata[[gg]]$est.tinter[1:Nt] <- colMeans(gdata[[gg]]$Samples.tinter)
    gdata[[gg]]$est.tinter[(Nt+1):(2*Nt)] <- colSDs(gdata[[gg]]$Samples.tinter)
  }
  
  gdata[[gg]]$est.iinter[1:K] <- colMeans(gdata[[gg]]$Samples.iinter)
  gdata[[gg]]$est.iinter[(K+1):(2*K)] <- colSDs(gdata[[gg]]$Samples.iinter)
}

#### Bayes factor (approximated through BICs) ####
cat("### Computing marginal likelihood estimates ###", "\n")

NN <- numeric(Ng)
NN.ESS <- numeric(Ng)
if (ftype == "BF") {
  Np <- Ng*Nt + Ng*K + (Ng-1) + (Ng-1)*(Nt-1) # Ng*Nt covariance parameter(s), Ng*K item intercepts, Ng-1 group intercepts, (Ng-1)*(Nt-1) subtest intercepts 
}  else if (ftype == "MD") {
  Np <- Ng*Nt + Ng*K + (Ng-1) + (Ng-1)*Nt # Ng*Nt covariance parameter(s), Ng*K item intercepts, Ng-1 group intercepts, (Ng-1)*Nt subtest intercepts 
}  else if (ftype == "UD") {
  Np <- Ng*Nt + Ng*K + (Ng-1) # Ng*Nt covariance parameter(s), Ng*K item intercepts, Ng-1 group intercepts 
}
# Estimate marginal likelihood
for (gg in 1:Ng) {
  # Posterior mean estimates
  tau.HMC.est <- gdata[[gg]]$est.tau[1:Nt] 
  sig2k.HMC.est <- gdata[[gg]]$est.sig2k[1:K] 
  ginter.HMC.est <- gdata[[gg]]$est.ginter[1]
  tinter.HMC.est <- gdata[[gg]]$est.tinter[1:Nt]
  iinter.HMC.est <- gdata[[gg]]$est.iinter[1:K]

  # BIC
  Sigma <- diag(sig2k.HMC.est)
  for (tt in 1:Nt) {
    Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * tau.HMC.est[tt]
  }
  NN[gg] <- out$gdata[[gg]]$N * K # Effective sample size, safe choice
  NN.ESS[gg] <- sum(diag(Sigma)) / sum(Sigma) * NN[gg] # Effective sample size, corrected for correlations
  
  mu.est <- iinter.HMC.est
  if (gg != 1) {
    mu.est <- mu.est + ginter.HMC.est
    for (tt in 1:Nt) {
      if (is.null(out$which.f1) || tt != out$which.f1) {
        mu.est <- mu.est + tinter.HMC.est[tt] * u[tt,] 
      }
    }
  }
  
  Z <- MASS::mvrnorm(out$gdata[[gg]]$N, mu.est, Sigma, empirical = TRUE)
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
  ll <- numeric(out$gdata[[gg]]$N)
  for (jj in 1:J) {
    for (kk in 1:K) {
      condition.on <- (1:K)[-kk] # condition on Z[-kk]
      B11 <- B11.list[[kk]]
      B12 <- B12.list[[kk]]
      B21 <- B21.list[[kk]]
      B22.inv <- B22.inv.list[[kk]]
      mu <- mu.est[kk] + mu.part.list[[kk]] %*% (matrix(t(Z[, condition.on]), nrow = K - 1, ncol = out$gdata[[gg]]$N) - matrix(mu.est[condition.on], nrow = K - 1, ncol = out$gdata[[gg]]$N, byrow = FALSE))
      sig2 <- sig2.list[[kk]] 
      Z[Y[[gg]][, kk] == 0, kk] <- extraDistr::rtnorm(n = length(mu[1, Y[[gg]][, kk] == 0]), mean = mu[1, Y[[gg]][, kk] == 0], sd = sqrt(sig2[1, 1]), a = -Inf, b = 0)
      Z[Y[[gg]][, kk] == 1, kk] <- extraDistr::rtnorm(n = length(mu[1, Y[[gg]][, kk] == 1]), mean = mu[1, Y[[gg]][, kk] == 1], sd = sqrt(sig2[1, 1]), a = 0, b = Inf)
      Z[, kk] <- Z[, kk] / sd(Z[, kk]) # Delta parameterization
    }
    ll <- ll + Rfast::dmvnorm(Z, mu = mu.est, sigma = Sigma, logged = FALSE) # Estimate likelihood of observed data
  }
  ll <- ll / J
  ll <- log(ll)
  ll <- sum(ll)
  #BIC <- Np * log(NN) - 2 * ll
  #BIC.ESS <- Np * log(NN.ESS) - 2 * ll
  marllik[gg] <- ll #-BIC / 2
  marllik.ESS[gg] <- ll # -BIC.ESS / 2
  
}
print(sum(marllik) - mean((1/2)*Np*log(NN)))
saveRDS(list(gdata = gdata, ll = marllik, marllik = sum(marllik) - mean((1/2)*Np*log(NN)), marllik.ESS = sum(marllik.ESS) - mean((1/2)*Np*log(NN.ESS))), file = filename)



