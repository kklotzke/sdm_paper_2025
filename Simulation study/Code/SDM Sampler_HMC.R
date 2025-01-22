# Implementation of the Hamiltonian Metropolis Hastings algorithm for Structured Dependence Modelling
# The model is a SDM for a bi-factor dependence structure with Nt (f1 + fx) factors
#  
# Konrad Klotzke

# Note: last layer must be f1
# dtype 1: MVN-distributed data, 2: binary data 
# pre.est 'mean': posterior mean of Gibbs sampler, 'last': last sample of Gibbs sampler, "random": random sample of Gibbs sampler as set of starting values for HMC
SamplerSDM_HMC <- function(Y, u, dtype = 1, epsilon, epsilon.scale, epsilon.sd, L, L.sd, M, enforceLB.tau = TRUE, pr.tau.lower = -2, pr.tau.upper = 2, pr.sig2k.upper = 10, pr.iinter.mu0 = 0, pr.iinter.var0 = 10, a0 = 0.0001, b0 = 0.0001, XG, XG.Gibbs, Burnin.Gibbs, runHMC = TRUE, runGibbs = TRUE, Gibbs.trmethod = "ICDF", Gibbs.tau.start = NULL, Gibbs.Z.start = NULL, pre.est = "mean", pre.iinter = NULL, pre.sig2k = NULL, pre.tau = NULL, pre.Z = NULL, save.Z = FALSE, save.Z.Gibbs = FALSE, save.DistPar = FALSE, ii.print = 100, ii.print.Gibbs = 100, verbose = FALSE) {
  N <- nrow(Y)
  K <- ncol(Y)
  Nt <- nrow(u)

  if (missing(epsilon.scale)) {
    epsilon.scale <- rep(1, Nt)
  }

  if (missing(epsilon.sd)) {
    epsilon.sd <- epsilon / 10
  }

  if (missing(L.sd)) {
    L.sd <- L / 10
  }

  if (length(epsilon.scale) == 1) {
    epsilon.scale <- rep(epsilon.scale, Nt)
  }

  if (length(epsilon.scale) != Nt) {
    cat("ERROR: length(epsilon.scale) must be equal to nrow(u).", "\n")
    return(NULL)
  }

  if (is.null(Gibbs.tau.start)) {
    Gibbs.tau.start <- rep(0.01, Nt)
  }

  if (is.null(Gibbs.Z.start)) {
    Gibbs.Z.start <- matrix(0, nrow = N, ncol = K)
  }
 
  if (dtype == 1) {
    Z <- Y
    # for (kk in 1:K) {
    #   Z[, kk] <- Z[, kk] / sd(Z[, kk]) # Standardize to var(Z_k) = 1 
    # }
  }
  else if (dtype == 2) {
    pr.tau.lower <- -1
    pr.tau.upper <- 1
    # Auxiliary latent responses
    Z <- Gibbs.Z.start
  }

  if (length(pr.tau.lower) == 1) {
    pr.tau.lower <- rep(pr.tau.lower, Nt)
  }
  if (length(pr.tau.upper) == 1) {
    pr.tau.upper <- rep(pr.tau.upper, Nt)
  }

  if (save.Z) {
    Z.list <- vector(mode = "list", length = XG)
  }
  else {
    Z.list <- NULL
  }

  if (save.Z.Gibbs) {
    Z.list.Gibbs <- vector(mode = "list", length = XG.Gibbs)
  }
  else {
    Z.list.Gibbs <- NULL
  }

  if (save.DistPar) {
    if (dtype == 1) {
      DistPar.sig2k <- vector(mode = "list", length = XG)
      DistPar.sig2k.Gibbs <- vector(mode = "list", length = XG.Gibbs)
    }
    else {
      DistPar.sig2k <- NULL
      DistPar.sig2k.Gibbs <- NULL
    }
    DistPar.tau.Gibbs <- vector(mode = "list", length = XG.Gibbs)
  }
  else {
    DistPar.sig2k <- NULL
    DistPar.sig2k.Gibbs <- NULL
    DistPar.tau.Gibbs <- NULL
  }

  # Sum of squares for covariance parameters
  SSB <- numeric(Nt)

  if (dtype == 1) {
    # Sum of squares for measurement error variance parameters
    SSWk <- numeric(K)
    for (kk in 1:K) {
      SSWk[kk] <- sum((Z[, kk] - mean(Z[, kk])) ^ 2)
    }
  }

  ## Run Gibbs-sampler to obtain starting values
  if (runGibbs) {

    Samples.iinter.Gibbs <- matrix(NA, nrow = XG.Gibbs, ncol = K) # Item intercepts
    Samples.sig2k.Gibbs <- matrix(NA, nrow = XG.Gibbs, ncol = K) # Measurement error variance parameters
    Samples.tau.Gibbs <- matrix(NA, nrow = XG.Gibbs, ncol = Nt)
    Samples.iinter.Gibbs[1,] <- 0
    if (dtype == 1) {
      Samples.sig2k.Gibbs[1,] <- 1
    }
    else if (dtype == 2) {
      Samples.sig2k.Gibbs[1,] <- 0.2
    }
    Samples.tau.Gibbs[1,] <- Gibbs.tau.start
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
        Z.list.Gibbs[[ii]] <- Z
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
        # Bi-factor structure
        scale <- sum(u[tt,]) ^ 2 / sum(as.logical(u[tt,])) ^ 2 # Rescale for non-binary loadings
        if (tt != Nt) {
          vec <- which(as.logical(u[tt,]))
          tr1 <- (-sum(sig2k[as.logical(u[tt,])]) / sum(as.logical(u[tt,])) ^ 2) / scale # Blockdiagonal structure
          # if (Nt != 1) {
          #   tr2 <- -1 / (t(u[tt, vec]) %*% solve(diag(sig2k[vec]) + tau[Nt]) %*% u[tt, vec])[1, 1] # lower bound can be smaller conditional on f1
          # }
          # else {
          #   tr2 <- Inf
          # }
          tr.lb[tt] <- tr1 #min(tr1, tr2)
          if (dtype == 2) {
            tr.ub[tt] <- 1 - tau[Nt] # ensure that var(Z_k) = 1 with all positive sig2k (delta parameterization)
          }
          shift <- -tr.lb[tt] + tau[Nt]
        }
        else {
          tr.lb[tt] <- -1 / (t(u[tt,]) %*% solve(Sigma) %*% u[tt,])[1, 1]
          if (dtype == 2) {
            if (Nt == 1) {
              tr.ub[tt] <- 1
            }
            else {
              tr.ub[tt] <- 1 - max(tau[1:(Nt - 1)]) # ensure that var(Z_k) = 1 with all positive sig2k (delta parameterization)
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

      # if (dtype == 1) {
      #   Z <- Y
      #   for (kk in 1:K) {
      #     kvar <- sig2k[kk]
      #     for (tt in 1:Nt) {
      #       kvar <- kvar + tau[tt] * u[tt, kk] ^ 2
      #     }
      #     Z[, kk] <- Z[, kk] / sqrt(kvar) # Standardize to var(Z_K) = 1
      #   }
      # }

      if (ii %% ii.print.Gibbs == 0) {
        cat("-- Gibbs-Sampling Iteration ", ii, "-- ", "\n")
        print(round(colMeans(Samples.tau.Gibbs[2:ii,]), 3))
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
    L.mean <- L
    epsilon.mean <- epsilon
    #epsilon.sd <- epsilon / 10
    mh.accept <- 0
    Samples.iinter <- matrix(NA, nrow = XG, ncol = K) # Item intercepts
    Samples.sig2k <- matrix(NA, nrow = XG, ncol = K) # Measurement error variance parameters
    Samples.tau <- matrix(NA, nrow = XG, ncol = Nt)
    Samples.p <- matrix(NA, nrow = XG, ncol = Nt)
    Samples.iinter[1,] <- iinter.est
    Samples.sig2k[1,] <- sig2k.est
    Samples.tau[1,] <- tau.est
    sig2k.shift <- numeric(K)
    if (missing(M)) {
      #M <- diag(Nt)
      if (Nt == 1) {
        x <- var(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs, 1])
        M <- matrix(1 / x, nrow = 1, ncol = 1)
      }
      else {
        x <- diag((cov(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs,])))
        M <- diag(1 / x)
      }
      print(1 / x)
    }
    #M <- diag( 1/( diag( (cov(Samples.tau.Gibbs[(Burnin.Gibbs + 1):XG.Gibbs,])) ) ) )
    Samples.p[1,] <- MASS::mvrnorm(n = 1, mu = rep(0, Nt), Sigma = M) # Initial momentum

    for (ii in 2:XG) {
      iinter <- Samples.iinter[ii - 1,] #seq(-1,1,length.out = K)
      sig2k <- Samples.sig2k[ii - 1,]

      # Current position parameter(s)
      tau <- Samples.tau[ii - 1,]

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

      if (dtype == 2) {
        sig2k <- NULL # Computed during evaluation 
      }

      if (save.Z) {
        Z.list[[ii]] <- Z
      }

      if (save.DistPar && (dtype == 1)) {
        DistPar.sig2k[[ii]] <- matrix(NA, nrow = K, ncol = 6) # shift, scale, shape, rate, a, b
      }

      # Sample momentum parameter(s)
      p <- MASS::mvrnorm(n = 1, mu = rep(0, Nt), Sigma = M)

      # Vary number of steps for better parameter space exploration / to prevent getting stuck in a loop
      L <- round(runif(1, L.mean - L.sd, L.mean + L.sd))

      # Vary step size for better parameter space exploration / to prevent getting stuck in a loop
      epsilon <- runif(1, epsilon.mean - epsilon.sd, epsilon.mean + epsilon.sd)

      # if (slimit < epsilon) {
      #   cat("WARNING: epsilon exceeds estimated stability limit (epsilon=", epsilon, ", slimit=", slimit, ").", "\n", sep = "")
      # }

      # Position parameter(s) to update
      tau.star <- tau

      # Momentum parameter(s) to update
      p.star <- p

      # Half step for momentum parameter(s)
      G <- getGrad(x = tau.star, Z = Z, sig2k = sig2k, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = sig2k.shift)
      p.star <- p.star - epsilon.scale * epsilon * G / 2
      
      ll <- 1
      while (ll <= L) {
        if ((ll == 1) | ((ll > 1) & !all(G == 0))) {
          if (any(!is.finite(tau.star)))
            browser()
          
          # Full step for position parameter(s)
          if (!enforceLB.tau) {
            tau.star <- tau.star + epsilon.scale * epsilon * p.star / diag(M)
          }
          else {
            for (tt in 1:Nt) {
              # Full step for position parameter(s) 
              tau.star[tt] <- tau.star[tt] + epsilon.scale[tt] * epsilon * p.star[tt] / (diag(M)[tt])
              lb <- getLB(x = tau.star, sig2k = sig2k, u = u, which.tt = tt)
              while (tau.star[tt] < lb) { # Enforce lower bound
                #browser()
                tau.star[tt] <- lb + (lb - tau.star[tt])
                p.star[tt] <- -p.star[tt]
              }
            }
          }
          
          # Full step for momentum parameter(s), unless last iteration
          if (ll != L) {
            G <- getGrad(x = tau.star, Z = Z, sig2k = sig2k, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = sig2k.shift)
            p.star <- p.star - epsilon.scale * epsilon * G
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
      G <- getGrad(x = tau.star, Z = Z, sig2k = sig2k, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = sig2k.shift)
      p.star <- p.star - epsilon.scale * epsilon * G / 2


      # Make proposal distribution symmetric
      p.star <- -p.star

      # cat("tau.star:", round(tau.star, 3), "|", "p.star:", round(p.star, 3), "\n")

      # Metropolis step
      U.min1 <- -lpost(x = tau, Z = Z, sig2k = sig2k, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = sig2k.shift)
      K.min1 <- sum((p ^ 2) / (2 * diag(M))) # (t(p) %*% solve(M) %*% p)[1, 1] / 2
      H.min1 <- U.min1 + K.min1
      U.cand <- -lpost(x = tau.star, Z = Z, sig2k = sig2k, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = sig2k.shift)
      K.cand <- sum((p.star ^ 2) / (2 * diag(M))) #  (t(p.star) %*% solve(M) %*% p.star)[1, 1] / 2
      H.cand <- U.cand + K.cand

      ar <- exp(H.min1 - H.cand)
      if (is.na(ar)) {
        ar <- 0
      }
      ar <- min(1, ar)
      un <- runif(1, 0, 1)
      if (un > ar) {
        # Reject
        Samples.tau[ii,] <- tau
        Samples.p[ii,] <- p
        if (verbose) {
          cat("ii:", ii, "|", "Udiff:", round(U.cand - U.min1, 2), "|", "Kdiff:"
                , round(K.cand - K.min1, 2)
                , "|", "REJECT", "\n")
          print(round(tau.star, 3))
        }
      }
      else {
        # Accept
        Samples.tau[ii,] <- tau <- tau.star
        Samples.p[ii,] <- p.star
        mh.accept <- mh.accept + 1
        if (verbose) {
          cat("ii:", ii, "|", "Udiff:", round(U.cand - U.min1, 2), "|", "Kdiff:"
                  , round(K.cand - K.min1, 2)
                  , "|", "ACCEPT", "\n")
          print(round(tau.star, 3))
        }
      }

      # Measurement error variance parameters (Gibbs-sampling step)
      for (kk in 1:K) {
        shift <- 0
        for (tt in 1:Nt) {
          shift <- shift + tau[tt] * u[tt, kk] ^ 2
        }
        sig2k.shift[kk] <- shift
        if (dtype == 1) {
          # Sample measurement error variance parameter from proper TS-IG posterior
          lower <- shift + max(-shift, 0) # sig2k > 0
          upper <- shift + pr.sig2k.upper
          if (save.DistPar) {
            DistPar.sig2k[[ii]][kk, 1] <- shift
            DistPar.sig2k[[ii]][kk, 2] <- 1
            DistPar.sig2k[[ii]][kk, 3] <- a0
            DistPar.sig2k[[ii]][kk, 4] <- b0
            DistPar.sig2k[[ii]][kk, 5] <- lower
            DistPar.sig2k[[ii]][kk, 6] <- upper
          }
          CD.L <- invgamma::pinvgamma(lower, shape = a0 + N / 2, rate = b0 + SSWk[kk] / 2, lower.tail = TRUE) # Lower bound sig2k = 0
          CD.U <- invgamma::pinvgamma(upper, shape = a0 + N / 2, rate = b0 + SSWk[kk] / 2, lower.tail = TRUE) # Upper bound sig2k = pr.sig2k.upper
          sig2k[kk] <- invgamma::qinvgamma(runif(1) * (CD.U - CD.L) + CD.L, shape = a0 + N / 2, rate = b0 + SSWk[kk] / 2) - shift
        }
        else if (dtype == 2) {
          # Compute measurement error variance parameter
          sig2k[kk] <- 1 - shift
        }
        Samples.sig2k[ii, kk] <- sig2k[kk]
      }

      # Update covariance structure
      Sigma <- diag(sig2k)
      for (tt in 1:Nt) {
        Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * tau[tt]
      }

      # Item intercepts (Gibbs-sampling step)
      Samples.iinter[ii,] <- Sample.IInter(Z = Z, Sigma1 = Sigma, mu0 = pr.iinter.mu0, var0 = pr.iinter.var0)
      if (dtype == 2) {
        Samples.iinter[ii,] <- Samples.iinter[ii,] - mean(Samples.iinter[ii,])
      }

      if (ii %% ii.print == 0) {
        cat("-- HMC Iteration ", ii, "-- ", "\n")
        cat("ESS:", round(coda::effectiveSize(coda::as.mcmc(Samples.tau[2:ii,])), 0), "\n")
        print(mh.accept)
        if (Nt == 1) {
          print(round(mean(Samples.tau[2:ii]), 3))
        }
        else {
          print(round(colMeans(Samples.tau[2:ii,]), 3))
        }
        #print(round(mean(Samples.tau[1:ii]), 4))
        #print(round(Samples.tau[ii,], 3))
        #print(sig2k)
      }
    }
  }
  else {
    mh.accept <- NULL
    Samples.iinter <- NULL
    Samples.sig2k <- NULL
    Samples.tau <- NULL
    Samples.p <- NULL
    M <- NULL
  }

  return(list("Samples.tau" = Samples.tau, "Samples.p" = Samples.p, "Samples.iinter" = Samples.iinter, "Samples.sig2k" = Samples.sig2k, "mh.accept" = mh.accept
  , "Samples.tau.Gibbs" = Samples.tau.Gibbs, "Samples.iinter.Gibbs" = Samples.iinter.Gibbs, "Samples.sig2k.Gibbs" = Samples.sig2k.Gibbs, "runGibbs" = runGibbs, "pre.est" = pre.est
  , "u" = u, "epsilon" = epsilon, "epsilon.scale" = epsilon.scale, "epsilon.sd" = "epsilon.sd", "L" = L, "L.sd" = L.sd, "M" = M, "enforceLB.tau" = enforceLB.tau, "dtype" = dtype, "slimit" = slimit
  , "pr.tau.lower" = pr.tau.lower, "pr.tau.upper" = pr.tau.upper, "pr.sig2k.upper" = pr.sig2k.upper, "pr.iinter.mu0" = pr.iinter.mu0, "pr.iinter.var0" = pr.iinter.var0, "a0" = a0, "b0" = b0
  , "Z.list" = Z.list, "Z.list.Gibbs" = Z.list.Gibbs, "DistPar.sig2k" = DistPar.sig2k, "DistPar.sig2k.Gibbs" = DistPar.sig2k.Gibbs, "DistPar.tau.Gibbs" = DistPar.tau.Gibbs))
}

# Log posterior 
lpost <- function(x, Z, sig2k = NULL, iinter, u, pr.tau.lower = -1, pr.tau.upper = 1, pr.sig2k.upper = NULL, pr.iinter.mu0, pr.iinter.var0, a0 = NULL, b0 = NULL, sig2k.shift = NULL, neg = FALSE
, pr.tau = "Uniform", tau.shift = NULL, tau.scale = NULL) {
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

  mu <- iinter
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
getGrad <- function(x, Z, sig2k, iinter, u, pr.tau.lower, pr.tau.upper, pr.sig2k.upper = NULL, pr.iinter.mu0, pr.iinter.var0, a0 = NULL, b0 = NULL, sig2k.shift = NULL) {
  G <- as.vector(rootSolve::gradient(f = lpost, x = x, Z = Z, sig2k = sig2k, iinter = iinter, u = u, pr.tau.lower = pr.tau.lower, pr.tau.upper = pr.tau.upper, pr.sig2k.upper = pr.sig2k.upper, pr.iinter.mu0 = pr.iinter.mu0, pr.iinter.var0 = pr.iinter.var0, a0 = a0, b0 = b0, sig2k.shift = sig2k.shift, neg = TRUE))
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


#####################################################################################################################
# Functions used for BF testing

# # pr.tau (prior for tau): "TS-IG", "Uniform"
# lpost.BCSM <- function(Sample, data, u, dtype, Z.list = NULL, DistPar.sig2k = NULL, DistPar.tau = NULL, pr.tau, pr.tau.lower = NULL, pr.tau.upper = NULL, pr.sig2k.upper = NULL, pr.iinter.mu0, pr.iinter.var0) {
#   Y <- data$Y
#   u <- data$u
#   N <- nrow(Y)
#   K <- ncol(Y)
#   Nt <- nrow(u)
#   ii <- Sample[1]
#   print(ii)
#   if (dtype == 1) {
#     if ((ii == as.integer(ii)) && (ii >= 1) && (ii <= length(DistPar.sig2k))) {
#       sfrom <- "post"
#     } else {
#       sfrom <- "prop"
#       ii <- sample(1:length(DistPar.sig2k), 1) # If sample is from proposal distribution, select random sample of Z
#     }
#   } else if (dtype == 2) {
#     if ((ii == as.integer(ii)) && (ii >= 1) && (ii <= length(Z.list))) {
#       sfrom <- "post"
#     } else {
#       sfrom <- "prop"
#       ii <- sample(1:length(Z.list), 1) # If sample is from proposal distribution, select random sample of Z
#     }
#   }
#   tau <- Sample[2:(Nt + 1)]
#   sig2k <- Sample[(Nt + 2):(Nt + K + 1)]
#   iinter <- Sample[(Nt + K + 2):length(Sample)]
#   if (dtype == 1) {
#     Z <- Y
#   } else if (dtype == 2) {
#     Z <- Z.list[[ii]] # Integrate over latent responses 
#   }
#   mu <- iinter
#   Sigma <- diag(sig2k)
#   for (tt in 1:Nt) {
#     Sigma <- Sigma + u[tt,] %*% t(u[tt,]) * tau[tt]
#   }
#   # Enforce PD
#   if (all(is.finite(Sigma)) && all(eigen(Sigma)$values > 0)) {
#     llik <- sum(Rfast::dmvnorm(Z, mu = mu, sigma = Sigma, logged = TRUE))
#     if (pr.tau == "TS-IG") {
#       if (sfrom == "post") {
#         llik <- llik + sum(dinvgamma_tr(tau, shift = DistPar.tau[[ii]][, 1], scale = DistPar.tau[[ii]][, 2], shape = DistPar.tau[[ii]][, 3], rate = DistPar.tau[[ii]][, 4], a = DistPar.tau[[ii]][, 5], b = DistPar.tau[[ii]][, 6], log = TRUE))
#       } else if (sfrom == "prop") {
#         llik <- llik + sum(dinvgamma_tr(tau, shift = 0, scale = 1, shape = DistPar.tau[[1]][, 3], rate = DistPar.tau[[1]][, 4], a = pr.tau.lower, b = pr.tau.upper, log = TRUE)) # Not nice
#       }
#     } else if (pr.tau == "Uniform") {
#       llik <- llik + sum(dunif(tau, min = pr.tau.lower, max = pr.tau.upper, log = TRUE))
#     }
#     llik <- llik + sum(dnorm(iinter, mean = pr.iinter.mu0, sd = sqrt(pr.iinter.var0), log = TRUE))
#     if (dtype == 1) {
#       if (sfrom == "post") {
#         llik <- llik + sum(dinvgamma_tr(sig2k, shift = DistPar.sig2k[[ii]][, 1], scale = DistPar.sig2k[[ii]][, 2], shape = DistPar.sig2k[[ii]][, 3], rate = DistPar.sig2k[[ii]][, 4], a = DistPar.sig2k[[ii]][, 5], b = DistPar.sig2k[[ii]][, 6], log = TRUE))
#       } else if (sfrom == "prop") {
#         llik <- llik + sum(dinvgamma_tr(sig2k, shift = 0, scale = 1, shape = DistPar.sig2k[[1]][, 3], rate = DistPar.sig2k[[1]][, 4], a = 0, b = pr.sig2k.upper, log = TRUE))
#       }
#     }
#   }
#   else {
#     llik <- -Inf # Not PD
#   }
#   return(llik)
# }

