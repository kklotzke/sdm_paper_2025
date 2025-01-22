# Implementation of the Hamiltonian Metropolis Hastings algorithm for Structured Dependence Modelling
# The model is a SDM with a UD dependence structure, a MD dependence structure or a bi-factor dependence structure
#  
# Group data is vertically scaled through the item intercepts
# Dependence structures are estimated freely per group using same classification matrix
#
# Konrad Klotzke

# dtype 1: MVN-distributed data, 2: binary data 
# pre.est 'mean': posterior mean of Gibbs sampler, 'last': last sample of Gibbs sampler, "random": random sample of Gibbs sampler as set of starting values for HMC
SamplerSDM_HMC_VS <- function(Y.list, u, dtype = 1, ftype = c("UD", "MD", "BF"), est.means = TRUE, epsilon, epsilon.scale, epsilon.sd, L, L.sd, M, enforceLB.tau = TRUE, pr.tau.lower = -2, pr.tau.upper = 2, pr.sig2k.upper = 10, pr.iinter.mu0 = 0, pr.iinter.var0 = 10, a0 = 0.0001, b0 = 0.0001, XG, XG.Gibbs, Burnin.Gibbs, runHMC = TRUE, runGibbs = TRUE, Gibbs.trmethod = "ICDF", Gibbs.tau.start = NULL, Gibbs.Z.start = NULL, pre.est = "mean", pre.iinter = NULL, pre.sig2k = NULL, pre.tau = NULL, pre.Z = NULL, save.Z = FALSE, save.Z.Gibbs = FALSE, save.DistPar = FALSE, ii.print = 100, ii.print.Gibbs = 100, verbose = FALSE) {

  Y <- Y.list[[1]] # Gibbs-sampler uses first group data to generate starting values for HMC
  N <- nrow(Y)
  K <- ncol(Y) # Balanced design
  Nt <- nrow(u)

  Ng <- length(Y.list)
  gdata <- list()
  for (gg in 1:Ng) {
    gdata[[gg]] <- list()
    gdata[[gg]]$N <- nrow(Y.list[[gg]])
    #gdata[[gg]]$K <- ncol(Y.list[[gg]])
    if (save.Z) {
      gdata[[gg]]$Z.vec <- matrix(NA, nrow = XG, ncol = gdata[[gg]]$N * K)
    }
    if (dtype == 1) {
      gdata[[gg]]$Z <- Y.list[[gg]]
      gdata[[gg]]$SSWk <- numeric(K)
      for (kk in 1:K) {
        gdata[[gg]]$SSWk[kk] <- sum((gdata[[gg]]$Z[, kk] - mean(gdata[[gg]]$Z[, kk])) ^ 2)
      }
      if (save.DistPar) {
        gdata[[gg]]$DistPar.sig2k <- vector(mode = "list", length = XG)
      }
    }
    if (missing(epsilon)) {
      gdata[[gg]]$epsilon <- 0.1
    }
    else {
      gdata[[gg]]$epsilon <- epsilon[[gg]]
    }
    if (missing(epsilon.scale)) {
      gdata[[gg]]$epsilon.scale <- rep(1, Nt)
    }
    else {
      gdata[[gg]]$epsilon.scale <- epsilon.scale[[gg]]
    }
    if (missing(epsilon.sd)) {
      gdata[[gg]]$epsilon.sd <- 0#epsilon / 10
    }
    else {
      gdata[[gg]]$epsilon.sd <- epsilon.sd[[gg]]
    }
    if (missing(L)) {
      gdata[[gg]]$L <- 6
    }
    else {
      gdata[[gg]]$L <- L[[gg]]
    }
    if (missing(L.sd)) {
        gdata[[gg]]$L.sd <- L[[gg]] / 10
    }
    else {
      gdata[[gg]]$L.sd <- L.sd[[gg]]
    }
    if (length(gdata[[gg]]$epsilon.scale) == 1) {
      gdata[[gg]]$epsilon.scale <- rep(gdata[[gg]]$epsilon.scale, Nt)
    }
    if (length(gdata[[gg]]$epsilon.scale) != Nt) {
      cat("ERROR: length(epsilon.scale) must be equal to nrow(u).", "\n")
      return(NULL)
    }
  }
  
  if (ftype == "BF") {
    which.f1 <- which(rowSums(u) == ncol(u))
  }
  else if (ftype == "UD") {
    which.f1 <- 1
  }
  else {
    which.f1 <- NULL
  }

  if (is.null(Gibbs.tau.start)) {
    Gibbs.tau.start <- rep(0.01, Nt)
  }

  if (is.null(Gibbs.Z.start)) {
    Gibbs.Z.start <- matrix(0, nrow = N, ncol = K)
  }

  if (dtype == 1) {
    Z <- Y
  }
  else if (dtype == 2) {
    pr.tau.lower <- max(-1, pr.tau.lower)
    pr.tau.upper <- min(1, pr.tau.upper)
    Z <- Gibbs.Z.start
  }

  if (length(pr.tau.lower) == 1) {
    pr.tau.lower <- rep(pr.tau.lower, Nt)
  }
  if (length(pr.tau.upper) == 1) {
    pr.tau.upper <- rep(pr.tau.upper, Nt)
  }

  if (save.Z.Gibbs) {
    Z.vec.Gibbs <- matrix(NA, nrow = XG.Gibbs, ncol = N * K)
  }
  else {
    Z.vec.Gibbs <- NULL
  }

  if (save.DistPar) {
    if (dtype == 1) {
      DistPar.sig2k.Gibbs <- vector(mode = "list", length = XG.Gibbs)
    }
    else {
      DistPar.sig2k.Gibbs <- NULL
    }
    DistPar.tau.Gibbs <- vector(mode = "list", length = XG.Gibbs)
  }
  else {
    DistPar.sig2k.Gibbs <- NULL
    DistPar.tau.Gibbs <- NULL
  }

  ## Run Gibbs-sampler to obtain starting values
  if (runGibbs) {

    # Sum of squares for covariance parameters
    SSB <- numeric(Nt)

    if (dtype == 1) {
      # Sum of squares for measurement error variance parameters
      SSWk <- numeric(K)
      for (kk in 1:K) {
        SSWk[kk] <- sum((Z[, kk] - mean(Z[, kk])) ^ 2)
      }
    }

    Samples.iinter.Gibbs <- matrix(NA, nrow = XG.Gibbs, ncol = K) # Item intercepts
    Samples.sig2k.Gibbs <- matrix(NA, nrow = XG.Gibbs, ncol = K) # Measurement error variance parameters
    Samples.tau.Gibbs <- matrix(NA, nrow = XG.Gibbs, ncol = Nt)
    Samples.iinter.Gibbs[1,] <- 0
    Samples.tau.Gibbs[1,] <- Gibbs.tau.start
    if (dtype == 1) {
      Samples.sig2k.Gibbs[1,] <- 1
    }
    else if (dtype == 2) {
      Samples.sig2k.Gibbs[1,] <- numeric(K)
      for (kk in 1:K) {
        shift <- 0
        for (tt in 1:Nt) {
          shift <- shift + Samples.tau.Gibbs[1, tt] * u[tt, kk] ^ 2
        }
        Samples.sig2k.Gibbs[1, kk] <- 1 - shift
      }
    }
    tr.lb <- tr.ub <- numeric(Nt)
    tr.ub <- pr.tau.upper

    for (ii in 2:XG.Gibbs) {
      iinter <- Samples.iinter.Gibbs[ii - 1,]
      sig2k <- Samples.sig2k.Gibbs[ii - 1,]
      tau <- Samples.tau.Gibbs[ii - 1,]

      # Mean structure
      B.sum <- iinter

      # Covariance structure
      Sigma <- diag(sig2k)
      for (tt in 1:Nt) {
        Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * tau[tt]
      }

      # Latent responses
      if (dtype == 2) {
        Z <- Sample.Z(Z = Z, Y = Y, B.sum = B.sum, Sigma = Sigma, par = "delta")
      }
      Z.star <- Z - matrix(B.sum, nrow = N, ncol = K, byrow = TRUE) # matrix(colMeans(Z), N, K, byrow = TRUE)#

      if (save.Z.Gibbs) {
        Z.vec.Gibbs[ii,] <- as.vector(Z)
      }

      if (save.DistPar) {
        if (dtype == 1) {
          DistPar.sig2k.Gibbs[[ii]] <- matrix(NA, nrow = K, ncol = 6) # shift, scale, shape, rate, a, b
        }
        DistPar.tau.Gibbs[[ii]] <- matrix(NA, nrow = Nt, ncol = 6) # shift, scale, shape, rate, a, b
      }

      # Update SSB
      for (tt in 1:(Nt)) {
        SSB[tt] <- sum((rowMeans(Z.star[, as.logical(u[tt,])]) - mean(Z.star[, as.logical(u[tt,])])) ^ 2)
      }

      Sigma <- diag(sig2k)
      for (tt in 1:Nt) {
        scale <- sum(u[tt,]) ^ 2 / sum(as.logical(u[tt,])) ^ 2 # Rescale for non-binary loadings
        if ((ftype == "BF" && tt != which.f1) || (ftype == "MD")) {
          vec <- which(as.logical(u[tt,]))
          tr.lb[tt] <- (-sum(sig2k[as.logical(u[tt,])]) / sum(as.logical(u[tt,])) ^ 2) / scale # Blockdiagonal structure
          if (dtype == 2) {
            if (ftype == "BF") {
              tr.ub[tt] <- min(1 - tau[which.f1], pr.tau.upper[tt]) # ensure that var(Z_k) = 1 with all positive sig2k (delta parameterization)
              shift <- -tr.lb[tt] + tau[which.f1]
            }
            else if (ftype == "MD") {
              tr.ub[tt] <- min(1, pr.tau.upper[tt])
              shift <- -tr.lb[tt]
            }
          }
        }
        else {
          tr.lb[tt] <- -1 / (t(u[tt,]) %*% solve(Sigma) %*% u[tt,])[1, 1]
          if (dtype == 2) {
            if (ftype != "BF") {
              tr.ub[tt] <- pr.tau.upper[tt]
            }
            else {
              tr.ub[tt] <- min(1 - max(tau[(1:Nt - 1)[-which.f1]]), pr.tau.upper[tt]) # ensure that var(Z_k) = 1 with all positive sig2k (delta parameterization)
            }
          }
          tmp <- 0
          for (xx in (1:Nt)[-tt]) {
            tmp <- tmp + tau[xx] * sum(u[xx,] * as.logical(u[tt,])) ^ 2
          }
          shift <- (sum(sig2k[as.logical(u[tt,])]) + tmp) / sum(as.logical(u[tt,])) ^ 2
        }
        shape <- a0 + N / 2
        rate <- b0 + SSB[tt] / 2
        if (save.DistPar) {
          DistPar.tau.Gibbs[[ii]][tt, 1] <- shift
          DistPar.tau.Gibbs[[ii]][tt, 2] <- scale
          DistPar.tau.Gibbs[[ii]][tt, 3] <- a0
          DistPar.tau.Gibbs[[ii]][tt, 4] <- b0
          DistPar.tau.Gibbs[[ii]][tt, 5] <- shift + max(-shift, scale * tr.lb[tt], scale * pr.tau.lower[tt])
          DistPar.tau.Gibbs[[ii]][tt, 6] <- shift + min(scale * tr.ub[tt], scale * pr.tau.upper[tt])
        }
        if (Gibbs.trmethod == "ICDF") {
          # Note: CD.L and CD.U are based on scale*tau, boundaries are rescaled
          CD.L <- invgamma::pinvgamma(shift + max(-shift, scale * tr.lb[tt], scale * pr.tau.lower[tt]), shape = shape, rate = rate, lower.tail = TRUE)
          CD.U <- invgamma::pinvgamma(shift + min(scale * tr.ub[tt], scale * pr.tau.upper[tt]), shape = shape, rate = rate, lower.tail = TRUE)
          Samples.tau.Gibbs[ii, tt] <- tau[tt] <- (invgamma::qinvgamma(runif(1) * (CD.U - CD.L) + CD.L, shape = shape, rate = rate) - shift) / scale
        }
        if ((Gibbs.trmethod == "RS") || !is.finite(tau[tt])) {
          tdraws <- (invgamma::rinvgamma(100000, shape = shape, rate = rate) - shift) / scale
          tdraws <- tdraws[is.finite(tdraws)]
          tdraws <- tdraws[(tdraws > max(tr.lb[tt], pr.tau.lower[tt])) & (tdraws < min(tr.ub[tt], pr.tau.upper[tt]))] # Truncate
          rr <- 0
          while (length(tdraws) == 0 & rr <= 5) {
            # Sample more if needed
            tdraws <- (invgamma::rinvgamma(100000, shape = shape, rate = rate) - shift) / scale
            tdraws <- x <- tdraws[is.finite(tdraws)]
            tdraws <- tdraws[(tdraws > max(tr.lb[tt], pr.tau.lower[tt])) & (tdraws < min(tr.ub[tt], pr.tau.upper[tt]))] # Truncate
            rr <- rr + 1
          }
          if (length(tdraws) == 0) {
            if (all(x <= max(tr.lb[tt], pr.tau.lower[tt]))) {
              tdraws <- runif(1, max(tr.lb[tt], pr.tau.lower[tt]), max(tr.lb[tt], pr.tau.lower[tt]) + 0.01) # Sampler is at edge of lower limit of parameter space
            }
            else if (all(x >= min(tr.ub[tt], pr.tau.upper[tt]))) {
              tdraws <- runif(1, min(tr.ub[tt], pr.tau.upper[tt]) - 0.01, min(tr.ub[tt], pr.tau.upper[tt])) # Sampler is at the edge of upper limit of parameter space
            }
            #tdraws <- Gibbs.tau.start[tt] # DEBUG, REMOVE LATER
          }
          Samples.tau.Gibbs[ii, tt] <- tau[tt] <- tdraws[1]
        }
        Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * tau[tt]
      }

      # Item intercepts
      Samples.iinter.Gibbs[ii,] <- Sample.IInter(Z = Z, Sigma1 = Sigma, mu0 = pr.iinter.mu0, var0 = pr.iinter.var0)
      if (dtype == 2) {
        Samples.iinter.Gibbs[ii,] <- Samples.iinter.Gibbs[ii,] - mean(Samples.iinter.Gibbs[ii,])
      }

      if (any(!is.finite(tau)))
        browser()

      # Measurement error variance parameters
      for (kk in 1:K) {
        shift <- 0
        for (tt in 1:Nt) {
          shift <- shift + tau[tt] * u[tt, kk] ^ 2
        }
        if (dtype == 1) {
          # Sample measurement error variance parameter from proper TS-IG posterior
          lower <- shift + max(-shift, 0) # sig2k > 0
          upper <- shift + pr.sig2k.upper
          if (save.DistPar) {
            DistPar.sig2k.Gibbs[[ii]][kk, 1] <- shift
            DistPar.sig2k.Gibbs[[ii]][kk, 2] <- 1
            DistPar.sig2k.Gibbs[[ii]][kk, 3] <- a0
            DistPar.sig2k.Gibbs[[ii]][kk, 4] <- b0
            DistPar.sig2k.Gibbs[[ii]][kk, 5] <- lower
            DistPar.sig2k.Gibbs[[ii]][kk, 6] <- upper
          }
          CD.L <- invgamma::pinvgamma(lower, shape = a0 + N / 2, rate = b0 + SSWk[kk] / 2, lower.tail = TRUE) # Lower bound sig2k = 0
          CD.U <- invgamma::pinvgamma(upper, shape = a0 + N / 2, rate = b0 + SSWk[kk] / 2, lower.tail = TRUE) # Upper bound sig2k = pr.sig2k.upper
          sig2k[kk] <- invgamma::qinvgamma(runif(1) * (CD.U - CD.L) + CD.L, shape = a0 + N / 2, rate = b0 + SSWk[kk] / 2) - shift
        }
        else if (dtype == 2) {
          # Compute measurement error variance parameter
          sig2k[kk] <- 1 - shift
        }
        Samples.sig2k.Gibbs[ii, kk] <- sig2k[kk]
      }

      if (ii %% ii.print.Gibbs == 0) {
        cat("-- Gibbs-Sampling Iteration ", ii, "-- ", "\n")
        if (Nt == 1) {
          print(round(mean(Samples.tau.Gibbs[2:ii,]), 3))
        }
        else {
          print(round(colMeans(Samples.tau.Gibbs[2:ii,]), 3))
        }
      }
    }

    if (Nt == 1) {
      slimit <- sqrt(var(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs, 1]))
    }
    else {
      slimit <- sqrt(min(eigen(cov(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs,]))$value))
    }
    if ((pre.est == "mean") || !runHMC) {
      iinter.est <- colMeans(Samples.iinter.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs,])
      sig2k.est <- colMeans(Samples.sig2k.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs,])
      if (Nt == 1) {
        tau.est <- mean(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs, 1])
      }
      else {
        tau.est <- colMeans(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs,])
      }
    }
    else if (pre.est == "last") {
      iinter.est <- Samples.iinter.Gibbs[XG.Gibbs,]
      sig2k.est <- Samples.sig2k.Gibbs[XG.Gibbs,]
      tau.est <- Samples.tau.Gibbs[XG.Gibbs,]
    }
    else if (pre.est == "random") {
      iinter.est <- Samples.iinter.Gibbs[sample((Burnin.Gibbs + 1):XG.Gibbs, 1),]
      sig2k.est <- Samples.sig2k.Gibbs[sample((Burnin.Gibbs + 1):XG.Gibbs, 1),]
      tau.est <- Samples.tau.Gibbs[sample((Burnin.Gibbs + 1):XG.Gibbs, 1),]
    }

    print(tau.est)
  }
  else {
    slimit <- NULL
    Samples.tau.Gibbs <- NULL

    if (is.null(pre.iinter))
      iinter.est <- rep(0, K)
    else
      iinter.est <- pre.iinter

    if (is.null(pre.sig2k))
      sig2k.est <- rep(0.2, K)
    else
      sig2k.est <- pre.sig2k

    if (is.null(pre.tau))
      tau.est <- rep(0.1, Nt)
    else
      tau.est <- pre.tau

    if (!is.null(pre.Z))
      Z <- pre.Z
  }

  ## Hamiltonian Monte Carlo
  if (runHMC) {

    # Initialize data structures and HMC parameters
    for (gg in 1:Ng) {
      gdata[[gg]]$mh.accept <- 0
      gdata[[gg]]$Samples.ginter <- matrix(NA, nrow = XG, ncol = 1) # Group means
      gdata[[gg]]$Samples.tinter <- matrix(NA, nrow = XG, ncol = Nt) # Factor means
      gdata[[gg]]$Samples.iinter <- matrix(NA, nrow = XG, ncol = K) # Item intercepts
      gdata[[gg]]$Samples.sig2k <- matrix(NA, nrow = XG, ncol = K) # Measurement error variance parameters
      gdata[[gg]]$Samples.tau <- matrix(NA, nrow = XG, ncol = Nt)
      gdata[[gg]]$Samples.p <- matrix(NA, nrow = XG, ncol = Nt)
      gdata[[gg]]$Samples.ginter[1,] <- 0
      gdata[[gg]]$Samples.tinter[1,] <- 0
      gdata[[gg]]$Samples.iinter[1,] <- iinter.est
      gdata[[gg]]$Samples.tau[1,] <- tau.est
      gdata[[gg]]$sig2k.shift <- numeric(K) # Will replace this later with better method to sample measurement error variances
      gdata[[gg]]$Z <- Z
      if (any(!is.finite(gdata[[gg]]$Z)))
        gdata[[gg]]$Z <- matrix(0, nrow = gdata[[gg]]$N, ncol = K)
      if (dtype == 1) {
        gdata[[gg]]$Samples.sig2k[1,] <- sig2k.est
      }
      else if (dtype == 2) {
        gdata[[gg]]$Samples.sig2k[1,] <- numeric(K)
        for (kk in 1:K) {
          shift <- 0
          for (tt in 1:Nt) {
            shift <- shift + gdata[[gg]]$Samples.tau[1, tt] * u[tt, kk] ^ 2
          }
          gdata[[gg]]$sig2k.shift[kk] <- shift
          gdata[[gg]]$Samples.sig2k[1, kk] <- 1 - shift
        }
      }
      if (missing(M)) {
        #M <- diag(Nt)
        if (Nt == 1) {
          x <- var(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs, 1])
          gdata[[gg]]$M <- matrix(1 / x, nrow = 1, ncol = 1)
        }
        else {
          x <- diag((cov(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs,])))
          gdata[[gg]]$M <- diag(1 / x)
        }
        print(1 / x)
      }
      gdata[[gg]]$Samples.p[1,] <- MASS::mvrnorm(n = 1, mu = rep(0, Nt), Sigma = gdata[[gg]]$M) # Initial momentum
    }

    for (ii in 2:XG) {
      ## Estimate covariance structure
      for (gg in 1:Ng) {
        L.mean <- gdata[[gg]]$L
        epsilon.mean <- gdata[[gg]]$epsilon
        
        ginter <- gdata[[gg]]$Samples.ginter[ii - 1,]
        tinter <- gdata[[gg]]$Samples.tinter[ii - 1,]
        iinter <- gdata[[gg]]$Samples.iinter[ii - 1,]
        sig2k <- gdata[[gg]]$Samples.sig2k[ii - 1,]
        free.means <- est.means && (gg != 1)

        # Current position parameter(s)
        tau <- gdata[[gg]]$Samples.tau[ii - 1,]

        # Mean structure 
        B.sum <- iinter + ginter
        for (tt in 1:Nt) {
          if (is.null(which.f1) || tt != which.f1) {
            B.sum <- B.sum + u[tt,] * tinter[tt]
          }
        }

        # Covariance structure
        Sigma <- diag(sig2k)
        for (tt in 1:Nt) {
          Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * tau[tt]
        }

        # Latent responses
        if (dtype == 2) {
          gdata[[gg]]$Z <- Sample.Z(Z = gdata[[gg]]$Z, Y = Y.list[[gg]], B.sum = B.sum, Sigma = Sigma, par = "delta")
        }

        if (dtype == 2) {
          sig2k <- NULL # Computed during evaluation 
        }

        if (save.Z) {
          gdata[[gg]]$Z.vec[ii,] <- as.vector(gdata[[gg]]$Z)
        }

        if (save.DistPar && (dtype == 1)) {
          DistPar.sig2k[[ii]] <- matrix(NA, nrow = K, ncol = 6) # shift, scale, shape, rate, a, b
        }

        # Sample momentum parameter(s)
        p <- MASS::mvrnorm(n = 1, mu = rep(0, Nt), Sigma = gdata[[gg]]$M)

        # Vary number of steps for better parameter space exploration / to prevent getting stuck in a loop
        L <- round(runif(1, L.mean - gdata[[gg]]$L.sd, L.mean + gdata[[gg]]$L.sd))

        # Vary step size for better parameter space exploration / to prevent getting stuck in a loop
        epsilon <- runif(1, epsilon.mean - gdata[[gg]]$epsilon.sd, epsilon.mean + gdata[[gg]]$epsilon.sd)

        # Position parameter(s) to update
        tau.star <- tau

        # Momentum parameter(s) to update
        p.star <- p

        # Half step for momentum parameter(s)
        G <- getGrad(x = tau.star, Z = gdata[[gg]]$Z, sig2k = sig2k, ginter = ginter, tinter = tinter, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = gdata[[gg]]$sig2k.shift, which.f1 = which.f1, free.means = free.means)
        p.star <- p.star - gdata[[gg]]$epsilon.scale * epsilon * G / 2

        ll <- 1
        while (ll <= L) {
          if ((ll == 1) | ((ll > 1) & !all(G == 0))) {
            if (any(!is.finite(tau.star)))
              browser()

            # Full step for position parameter(s)
            if (!enforceLB.tau) {
              tau.star <- tau.star + gdata[[gg]]$epsilon.scale * epsilon * p.star / diag(gdata[[gg]]$M)
            }
            else {
              for (tt in 1:Nt) {
                # Full step for position parameter(s) 
                tau.star[tt] <- tau.star[tt] + gdata[[gg]]$epsilon.scale[tt] * epsilon * p.star[tt] / (diag(gdata[[gg]]$M)[tt])
                lb <- getLB(x = tau.star, sig2k = sig2k, u = u, which.tt = tt)
                while (tau.star[tt] < lb) {
                  # Enforce lower bound
                  tau.star[tt] <- lb + (lb - tau.star[tt])
                  p.star[tt] <- -p.star[tt]
                }
              }
            }

            # Full step for momentum parameter(s), unless last iteration
            if (ll != L) {
              G <- getGrad(x = tau.star, Z = gdata[[gg]]$Z, sig2k = sig2k, ginter = ginter, tinter = tinter, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = gdata[[gg]]$sig2k.shift, which.f1 = which.f1, free.means = free.means)
              p.star <- p.star - gdata[[gg]]$epsilon.scale * epsilon * G
            }
            if (any(!is.finite(p.star)))
              browser()

            # cat("#########################################\n")
            # print(round(tau.star, 2))
            # print(round(p.star, 2))
            # print(round(G, 2))
            # cat("#########################################\n")
            # browser()
            # if (any(G == -1000) | any(G == 1000))
            #   print("inf")
            # if (any(G == 0)) {
            #   print(which(G==0))
            #   #browser()
            # }
          }
          ll <- ll + 1
        }

        # Half step for momentum parameter(s)
        G <- getGrad(x = tau.star, Z = gdata[[gg]]$Z, sig2k = sig2k, ginter = ginter, tinter = tinter, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = gdata[[gg]]$sig2k.shift, which.f1 = which.f1, free.means = free.means)
        p.star <- p.star - gdata[[gg]]$epsilon.scale * epsilon * G / 2


        # Make proposal distribution symmetric
        p.star <- -p.star

        # cat("tau.star:", round(tau.star, 3), "|", "p.star:", round(p.star, 3), "\n")

        # Metropolis step
        U.min1 <- -lpost(x = tau, Z = gdata[[gg]]$Z, sig2k = sig2k, ginter = ginter, tinter = tinter, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = gdata[[gg]]$sig2k.shift, which.f1 = which.f1, free.means = free.means)
        K.min1 <- sum((p ^ 2) / (2 * diag(gdata[[gg]]$M))) # (t(p) %*% solve(M) %*% p)[1, 1] / 2
        H.min1 <- U.min1 + K.min1
        U.cand <- -lpost(x = tau.star, Z = gdata[[gg]]$Z, sig2k = sig2k, ginter = ginter, tinter = tinter, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = gdata[[gg]]$sig2k.shift, which.f1 = which.f1, free.means = free.means)
        K.cand <- sum((p.star ^ 2) / (2 * diag(gdata[[gg]]$M))) #  (t(p.star) %*% solve(M) %*% p.star)[1, 1] / 2
        H.cand <- U.cand + K.cand

        ar <- exp(H.min1 - H.cand)
        if (is.na(ar)) {
          ar <- 0
        }
        ar <- min(1, ar)
        un <- runif(1, 0, 1)
        if (un > ar) {
          # Reject
          gdata[[gg]]$Samples.tau[ii,] <- tau
          gdata[[gg]]$Samples.p[ii,] <- p
          if (verbose) {
            cat("ii:", ii, "|", "Udiff:", round(U.cand - U.min1, 2), "|", "Kdiff:"
                , round(K.cand - K.min1, 2)
                , "|", "REJECT", "\n")
            print(round(tau.star, 3))
          }
        }
        else {
          # Accept
          gdata[[gg]]$Samples.tau[ii,] <- tau <- tau.star
          gdata[[gg]]$Samples.p[ii,] <- p.star
          gdata[[gg]]$mh.accept <- gdata[[gg]]$mh.accept + 1
          if (verbose) {
            cat("ii:", ii, "|", "Udiff:", round(U.cand - U.min1, 2), "|", "Kdiff:"
                  , round(K.cand - K.min1, 2)
                  , "|", "ACCEPT", "\n")
            print(round(tau.star, 3))
          }
        }

        # Measurement error variance parameters 
        for (kk in 1:K) {
          shift <- 0
          for (tt in 1:Nt) {
            shift <- shift + tau[tt] * u[tt, kk] ^ 2
          }
          if (dtype == 1) {
            # Sample measurement error variance parameter from proper TS-IG posterior
            lower <- shift + max(-shift, 0) # sig2k > 0
            upper <- shift + pr.sig2k.upper
            if (save.DistPar) {
              gdata[[gg]]$DistPar.sig2k[[ii]][kk, 1] <- shift
              gdata[[gg]]$DistPar.sig2k[[ii]][kk, 2] <- 1
              gdata[[gg]]$DistPar.sig2k[[ii]][kk, 3] <- a0
              gdata[[gg]]$DistPar.sig2k[[ii]][kk, 4] <- b0
              gdata[[gg]]$DistPar.sig2k[[ii]][kk, 5] <- lower
              gdata[[gg]]$DistPar.sig2k[[ii]][kk, 6] <- upper
            }
            CD.L <- invgamma::pinvgamma(lower, shape = a0 + gdata[[gg]]$N / 2, rate = b0 + gdata[[gg]]$SSWk[kk] / 2, lower.tail = TRUE) # Lower bound sig2k = 0
            CD.U <- invgamma::pinvgamma(upper, shape = a0 + gdata[[gg]]$N / 2, rate = b0 + gdata[[gg]]$SSWk[kk] / 2, lower.tail = TRUE) # Upper bound sig2k = pr.sig2k.upper
            sig2k[kk] <- invgamma::qinvgamma(runif(1) * (CD.U - CD.L) + CD.L, shape = a0 + gdata[[gg]]$N / 2, rate = b0 + gdata[[gg]]$SSWk[kk] / 2) - shift
          }
          else if (dtype == 2) {
            # Compute measurement error variance parameter
            sig2k[kk] <- 1 - shift
          }
          gdata[[gg]]$sig2k.shift[kk] <- shift
          gdata[[gg]]$Samples.sig2k[ii, kk] <- sig2k[kk]
        }

        # Update covariance structure
        Sigma <- diag(sig2k)
        for (tt in 1:Nt) {
          Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * tau[tt]
        }
        gdata[[gg]]$Sigma <- Sigma
      }

      if (est.means) {
        ## Estimate mean structure
        # Group means 
        gdata[[1]]$Samples.ginter[ii, 1] <- 0 # Reference group
        ginter <- numeric(Ng)
        for (gg in 2:Ng) {
          B.sum <- gdata[[gg]]$Samples.iinter[ii - 1,]
          for (tt in 1:Nt) {
            if (is.null(which.f1) || tt != which.f1) {
              B.sum <- B.sum + u[tt,] * gdata[[gg]]$Samples.tinter[ii - 1, tt]
            }
          }
          gdata[[gg]]$Samples.ginter[ii, 1] <- ginter[gg] <- Sample.TInter(Z = gdata[[gg]]$Z - matrix(B.sum, nrow = gdata[[gg]]$N, ncol = K, byrow = TRUE), u = matrix(1, nrow = 1, ncol = K), Sigma1 = gdata[[gg]]$Sigma, mu0 = pr.iinter.mu0, var0 = pr.iinter.var0)
        }

        # Factor means
        gdata[[1]]$Samples.tinter[ii,] <- 0 # Reference group
        tinter <- matrix(NA, nrow = Nt, ncol = Ng)
        if (Nt == 1 && !is.null(which.f1)) {
          tinter[1,] <- ginter # Factor mean is equal to group mean if all items load on single factor 
          for (gg in 2:Ng) {
            gdata[[gg]]$Samples.tinter[ii, 1] <- tinter[1, gg] 
          }
        }
        else if (Nt > 1) {
          for (gg in 2:Ng) {
            B.sum <- gdata[[gg]]$Samples.iinter[ii - 1,] + gdata[[gg]]$Samples.ginter[ii, 1]
            if (!is.null(which.f1)) {
              tinter[-which.f1, gg] <- Sample.TInter(Z = gdata[[gg]]$Z - matrix(B.sum, nrow = gdata[[gg]]$N, ncol = K, byrow = TRUE), u = u[-which.f1,], Sigma1 = gdata[[gg]]$Sigma, mu0 = pr.iinter.mu0, var0 = pr.iinter.var0)
              tinter[which.f1, gg] <- gdata[[gg]]$Samples.ginter[ii, 1] # Factor mean is equal to group mean if all items load on single factor
            }
            else {
              tinter[, gg] <- Sample.TInter(Z = gdata[[gg]]$Z - matrix(B.sum, nrow = gdata[[gg]]$N, ncol = K, byrow = TRUE), u = u, Sigma1 = gdata[[gg]]$Sigma, mu0 = pr.iinter.mu0, var0 = pr.iinter.var0)
            }
          }
          for (gg in 2:Ng) {
            if(!is.null(which.f1)) {
              gdata[[gg]]$Samples.tinter[ii, -which.f1] <- tinter[-which.f1, gg] - mean(tinter[-which.f1, gg]) # Identify model
              gdata[[gg]]$Samples.tinter[ii, which.f1] <- tinter[which.f1, gg]
            }
            else {
              gdata[[gg]]$Samples.tinter[ii, ] <- tinter[, gg] - mean(tinter[, gg]) # Identify model
            }
          }
        }
      }
      else {
        ginter <- numeric(Ng)
        tinter <- matrix(0, nrow = Nt, ncol = Ng)
        for (gg in 1:Ng) {
          gdata[[gg]]$Samples.ginter[ii, 1] <- 0
          gdata[[gg]]$Samples.tinter[ii,] <- 0
        }
      }
    
      # Item intercepts 
      Samples.iinter.gmean <- numeric(K)
      for (gg in 1:Ng) {
        B.sum <- numeric(K)
        if (gg != 1) {
          for (tt in 1:Nt) {
            B.sum <- B.sum + u[tt,] * gdata[[gg]]$Samples.tinter[ii, tt]
          }
          if (is.null(which.f1)) {
            # Otherwise factor mean of f1 = group mean
            B.sum <- B.sum + gdata[[gg]]$Samples.ginter[ii, 1]
          }
        }
        Samples.iinter.gmean <- Samples.iinter.gmean + Sample.IInter(Z = gdata[[gg]]$Z - matrix(B.sum, nrow = gdata[[gg]]$N, ncol = K, byrow = TRUE), Sigma1 = gdata[[gg]]$Sigma, mu0 = pr.iinter.mu0, var0 = pr.iinter.var0)
      }
      Samples.iinter.gmean <- Samples.iinter.gmean / Ng
      #Samples.iinter.gmean <- Samples.iinter.gmean - mean(Samples.iinter.gmean) # Identify model
      for (gg in 1:Ng) {
        gdata[[gg]]$Samples.iinter[ii,] <- Samples.iinter.gmean
      }

      if (ii %% ii.print == 0) {
        cat("-------- HMC Iteration ", ii, "-------- ", "\n")
        for (gg in 1:Ng) {
          cat("##", "Group", gg, "##", "\n")
          cat("ESS:", round(coda::effectiveSize(coda::as.mcmc(gdata[[gg]]$Samples.tau[2:ii,])), 0), "\n")
          cat("mh.accept:", gdata[[gg]]$mh.accept, "\n")
          if (Nt == 1) {
            cat("tau:", round(mean(gdata[[gg]]$Samples.tau[2:ii]), 3), "\n")
          }
          else {
            cat("tau:", round(colMeans(gdata[[gg]]$Samples.tau[2:ii,]), 3), "\n")
          }
          if (Nt == 1) {
            cat("ginter:", round(mean(gdata[[gg]]$Samples.ginter[2:ii]), 3), "\n")
          }
          else {
            cat("tinter:", round(colMeans(gdata[[gg]]$Samples.tinter[2:ii,]), 3), "\n")
          }
          cat("###################", "\n")
        }
      }
    }
    # End iteration 
  }
  else {
    mh.accept <- NULL
    Samples.iinter <- NULL
    Samples.sig2k <- NULL
    Samples.tau <- NULL
    Samples.p <- NULL
    M <- NULL
  }

  return(list("gdata" = gdata, "est.means" = est.means, "which.f1" = which.f1
  , "Samples.tau.Gibbs" = Samples.tau.Gibbs, "Samples.iinter.Gibbs" = Samples.iinter.Gibbs, "Samples.sig2k.Gibbs" = Samples.sig2k.Gibbs, "runGibbs" = runGibbs, "pre.est" = pre.est
  , "u" = u, "epsilon" = epsilon, "epsilon.scale" = epsilon.scale, "epsilon.sd" = "epsilon.sd", "L" = L, "L.sd" = L.sd, "enforceLB.tau" = enforceLB.tau, "dtype" = dtype, "ftype" = ftype, "slimit" = slimit
  , "pr.tau.lower" = pr.tau.lower, "pr.tau.upper" = pr.tau.upper, "pr.sig2k.upper" = pr.sig2k.upper, "pr.iinter.mu0" = pr.iinter.mu0, "pr.iinter.var0" = pr.iinter.var0, "a0" = a0, "b0" = b0
  , "Z.vec.Gibbs" = Z.vec.Gibbs, "DistPar.sig2k.Gibbs" = DistPar.sig2k.Gibbs, "DistPar.tau.Gibbs" = DistPar.tau.Gibbs))
}

# Log posterior 
lpost <- function(x, Z, sig2k = NULL, ginter, tinter, iinter, u, pr.tau.lower = -1, pr.tau.upper = 1, pr.sig2k.upper = NULL, pr.iinter.mu0, pr.iinter.var0, a0 = NULL, b0 = NULL, sig2k.shift = NULL, neg = FALSE
, pr.tau = "Uniform", tau.shift = NULL, tau.scale = NULL, which.f1 = NULL, free.means) {
  K <- ncol(Z)
  Nt <- length(x)
  if (length(pr.tau.lower) == 1)
    pr.tau.lower <- rep(pr.tau.lower, Nt)
  if (length(pr.tau.upper) == 1)
    pr.tau.upper <- rep(pr.tau.upper, Nt)

  # Compute K measurement error variance parameters (delta parameterization)
  if (is.null(sig2k)) {
    dtype <- 2
    sig2k <- numeric(K)
    for (kk in 1:K) {
      shift <- 0
      for (tt in 1:Nt) {
        shift <- shift + x[tt] * u[tt, kk] ^ 2
      }
      sig2k[kk] <- 1 - shift
    }
  } else {
    dtype <- 1
  }
  # if (ginter != 0) {
  #   browser()
  # }
  mu <- iinter
  if (free.means) {
    mu <- mu + ginter
    for (tt in 1:Nt) {
      if (is.null(which.f1) || tt != which.f1) {
        mu <- mu + tinter[tt] * u[tt,] 
      }
    }
  }
  Sigma <- diag(sig2k)
  for (tt in 1:Nt) {
    Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * x[tt]
  }
  # Enforce PD
  if (all(is.finite(Sigma)) && all(eigen(Sigma)$values > 0)) {
    llik <- sum(Rfast::dmvnorm(Z, mu = mu, sigma = Sigma, logged = TRUE)) # Likelihood data
    if (pr.tau == "Uniform") {
      llik <- llik + sum(dunif(x, min = pr.tau.lower, max = pr.tau.upper, log = TRUE)) # Uniform prior covariances
    } else if (pr.tau == "TS-IG") {
      llik <- llik + sum(dinvgamma_tr(x, shift = tau.shift, scale = tau.scale, shape = a0, rate = b0, a = pr.tau.lower, b = pr.tau.upper, log = TRUE))
    }
    if (free.means) {
      llik <- llik + sum(dnorm(ginter, mean = pr.iinter.mu0, sd = sqrt(pr.iinter.var0), log = TRUE)) # Normal prior group means
      if (is.null(which.f1)) {
        llik <- llik + sum(dnorm(tinter, mean = pr.iinter.mu0, sd = sqrt(pr.iinter.var0), log = TRUE)) # Normal prior factor means 
      }
      else if (!is.null(which.f1) && Nt > 1) {
        llik <- llik + sum(dnorm(tinter[-which.f1], mean = pr.iinter.mu0, sd = sqrt(pr.iinter.var0), log = TRUE)) # Normal prior factor means, exclude f1 mean
      }
    }
    llik <- llik + sum(dnorm(iinter, mean = pr.iinter.mu0, sd = sqrt(pr.iinter.var0), log = TRUE)) # Normal prior item intercepts 
    if (dtype == 1) {
      llik <- llik + sum(dinvgamma_tr(sig2k, shift = sig2k.shift, shape = a0, rate = b0, a = sig2k.shift + pmax(-sig2k.shift, 0), b = sig2k.shift + pr.sig2k.upper, log = TRUE)) # TS-IG prior measurement error variances
    }
  }
  else {
    llik <- -Inf # Not PD
  }
  if (neg)
    llik <- -1 * llik
  return(llik)
}

# Compute gradient of log posterior with respect to tau
getGrad <- function(x, Z, sig2k, ginter, tinter, iinter, u, pr.tau.lower, pr.tau.upper, pr.sig2k.upper = NULL, pr.iinter.mu0, pr.iinter.var0, a0 = NULL, b0 = NULL, sig2k.shift = NULL, which.f1 = NULL, free.means) {
  G <- as.vector(rootSolve::gradient(f = lpost, x = x, Z = Z, sig2k = sig2k, ginter = ginter, tinter = tinter, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = sig2k.shift, neg = TRUE, which.f1 = which.f1, free.means = free.means))
  # if (any(!is.finite(G))) {
  #   print(G)
  # }
  G[which(is.na(G))] <- 0 # x is outside of prior bounds
  G[which(G == -Inf)] <- -1000
  G[which(G == Inf)] <- 1000
  return(G)
}

# Get lower bound of tau[tt]
getLB <- function(x, sig2k, u, which.tt) {
  Nt <- nrow(u)
  K <- ncol(u)

  # Compute K measurement error variance parameters (delta parameterization)
  if (is.null(sig2k)) {
    sig2k <- numeric(K)
    for (kk in 1:K) {
      shift <- 0
      for (tt in 1:Nt) {
        shift <- shift + x[tt] * u[tt, kk] ^ 2
      }
      sig2k[kk] <- 1 - shift
    }
  }

  # Construct covariance matrix up to (not including) layer tt
  Sigma <- diag(sig2k)
  if (which.tt > 1) {
    for (tt in 1:(which.tt - 1)) {
      Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * x[tt]
    }
  }

  return(-1 / (t(u[which.tt,]) %*% solve(Sigma) %*% u[which.tt,])[1, 1])
}

# Sample latent responses
Sample.Z <- function(Z, Y, B.sum, Sigma, par) {
  N <- nrow(Z)
  K <- ncol(Z)
  for (kk in sample(1:K)) {
    condition.on <- (1:K)[-kk] # condition on Z[-kk]
    Sigma.inv <- solve(Sigma[condition.on, condition.on]) # faster, as it is coded in rcpp
    B11 <- Sigma[kk, kk]
    B12 <- Sigma[kk, condition.on]
    B21 <- Sigma[condition.on, kk]
    B22.inv <- Sigma.inv
    mu <- B.sum[kk] + B12 %*% B22.inv %*% (matrix(t(Z[, condition.on]), nrow = K - 1, ncol = N) - matrix(B.sum[condition.on], nrow = K - 1, ncol = N, byrow = FALSE))
    sig2 <- B11 - B12 %*% B22.inv %*% B21
    Z[Y[, kk] == 0, kk] <- extraDistr::rtnorm(n = length(mu[1, Y[, kk] == 0]), mean = mu[1, Y[, kk] == 0], sd = sqrt(sig2[1, 1]), a = -Inf, b = 0)
    Z[Y[, kk] == 1, kk] <- extraDistr::rtnorm(n = length(mu[1, Y[, kk] == 1]), mean = mu[1, Y[, kk] == 1], sd = sqrt(sig2[1, 1]), a = 0, b = Inf)
    if (par == "delta") {
      Z[, kk] <- Z[, kk] / sd(Z[, kk]) # Delta parameterization
    }
  }
  return(Z)
}

# Sample item intercepts
Sample.IInter <- function(Z, Sigma1, mu0 = 0, var0 = 10 ^ 10) {
  N <- nrow(Z)
  K <- ncol(Z)
  var1 <- 1 / (N / diag(Sigma1) + 1 / var0)
  mu1 <- var1 * (N * colMeans(Z) / diag(Sigma1) + mu0 / var0)
  draw <- rnorm(K, mu1, sqrt(var1))
  #draw <- draw - mean(draw)
  return(draw)
}

# Sample group or cluster intercepts
Sample.TInter <- function(Z, u, Sigma1, mu0 = 0, var0 = 10 ^ 10) {
  N <- nrow(Z)
  K <- ncol(u)
  Nt <- nrow(u)

  draw <- numeric(Nt)
  for (tt in 1:Nt) {
    rowmeans <- rowMeans(Z[, as.logical(u[tt,])])
    rowmeans.var <- sum(Sigma1[as.logical(u[tt,]), as.logical(u[tt,])]) / sum(as.logical(u[tt,])) ^ 2
    var1 <- 1 / (N / rowmeans.var + 1 / var0)
    mu1 <- var1 * (N * mean(rowmeans) / rowmeans.var + mu0 / var0)
    draw[tt] <- rnorm(1, mu1, sqrt(var1))
    #Z <- Z - matrix(draw[tt] * as.logical(u[tt,]), nrow = N, ncol = K, byrow = TRUE)
  }
  return(draw)
}

dinvgamma_tr <- function(x, shift = 0, scale = 1, shape, rate, a, b, log = FALSE) {
  x <- scale * x + shift

  if (length(shape) == 1) {
    shape <- rep(shape, length(x))
  }
  if (length(rate) == 1) {
    rate <- rep(rate, length(x))
  }
  if (length(a) == 1) {
    a <- rep(a, length(x))
  }
  if (length(b) == 1) {
    b <- rep(b, length(x))
  }

  if (length(shape) != length(x)) {
    stop("ERROR: 'shape' must be a number or a vector of length(x).")
  }
  if (length(rate) != length(x)) {
    stop("ERROR: 'rate' must be a number or a vector of length(x).")
  }
  if (length(a) != length(x)) {
    stop("ERROR: 'a' must be a number or a vector of length(x).")
  }
  if (length(b) != length(x)) {
    stop("ERROR: 'b' must be a number or a vector of length(x).")
  }

  #p <- ((rate^shape) / gamma(shape)) * x ^(-shape-1) * exp(-rate/x) 
  p <- shape * log(rate) - lgamma(shape) + (-shape - 1) * log(x) - rate / x
  cp.a <- pgamma(1 / a, shape, rate, lower.tail = FALSE, log.p = FALSE)
  cp.b <- pgamma(1 / b, shape, rate, lower.tail = FALSE, log.p = FALSE)

  id <- (cp.b == cp.a)
  p[id] <- -Inf
  p[!id] <- p[!id] - log(cp.b[!id] - cp.a[!id])
  p[which(!((x > a) & (x < b)))] <- -Inf

  if (log) {
    return(p)
  }
  else {
    return(exp(p))
  }
}

colSDs <- function(x) sqrt(rowMeans((t(x) - colMeans(x)) ^ 2) * ((dim(x)[1]) / (dim(x)[1] - 1)))

