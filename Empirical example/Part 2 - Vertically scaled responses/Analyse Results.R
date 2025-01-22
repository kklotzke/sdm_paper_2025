source("../Tools.R")

library(coda)

results.UD <- readRDS("Results/Run_UD_VS.RDS")
results.MD6 <- readRDS("Results/Run_MD6_VS.RDS")
results.MD2 <- readRDS("Results/Run_MD2_VS.RDS")
results.BF6 <- readRDS("Results/Run_BF6_VS.RDS")
results.BF2 <- readRDS("Results/Run_BF2_VS.RDS")

## Correlation parameters
# UD
Nt <- ncol(results.UD$gdata[[1]]$est.tau) / 2
table_tau.UD <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_SD.UD <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_HPD0.UD <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tau.UD[gg, 1:Nt] <- results.UD$gdata[[gg]]$est.tau[,1:Nt]
  table_tau_SD.UD[gg, 1:Nt] <- results.UD$gdata[[gg]]$est.tau[,(Nt+1):(2*Nt)]
  table_tau_HPD0.UD[gg, Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.UD$gdata[[gg]]$Samples.tau[1, ], results.UD$gdata[[gg]]$Samples.tau[2, ])), prob = 0.95), 1, hasNotZero))
}
# MD6
Nt <- ncol(results.MD6$gdata[[1]]$est.tau) / 2
table_tau.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_SD.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_HPD0.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tau.MD6[gg, 1:Nt] <- results.MD6$gdata[[gg]]$est.tau[,1:Nt]
  table_tau_SD.MD6[gg, 1:Nt] <- results.MD6$gdata[[gg]]$est.tau[,(Nt+1):(2*Nt)]
  table_tau_HPD0.MD6[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.MD6$gdata[[gg]]$Samples.tau), prob = 0.95), 1, hasNotZero))
}
# MD2
Nt <- ncol(results.MD2$gdata[[1]]$est.tau) / 2
table_tau.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_SD.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_HPD0.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tau.MD2[gg, 1:Nt] <- results.MD2$gdata[[gg]]$est.tau[,1:Nt]
  table_tau_SD.MD2[gg, 1:Nt] <- results.MD2$gdata[[gg]]$est.tau[,(Nt+1):(2*Nt)]
  table_tau_HPD0.MD2[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.MD2$gdata[[gg]]$Samples.tau), prob = 0.95), 1, hasNotZero))
}
# BF6
Nt <- ncol(results.BF6$gdata[[1]]$est.tau) / 2
table_tau.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_SD.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_HPD0.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tau.BF6[gg, 1:Nt] <- results.BF6$gdata[[gg]]$est.tau[,1:Nt]
  table_tau_SD.BF6[gg, 1:Nt] <- results.BF6$gdata[[gg]]$est.tau[,(Nt+1):(2*Nt)]
  table_tau_HPD0.BF6[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.BF6$gdata[[gg]]$Samples.tau), prob = 0.95), 1, hasNotZero))
}
# BF2
Nt <- ncol(results.BF2$gdata[[1]]$est.tau) / 2
table_tau.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_SD.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tau_HPD0.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tau.BF2[gg, 1:Nt] <- results.BF2$gdata[[gg]]$est.tau[,1:Nt]
  table_tau_SD.BF2[gg, 1:Nt] <- results.BF2$gdata[[gg]]$est.tau[,(Nt+1):(2*Nt)]
  table_tau_HPD0.BF2[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.BF2$gdata[[gg]]$Samples.tau), prob = 0.95), 1, hasNotZero))
}

# Print results
round(rbind(table_tau.UD, rep(NA, ncol(table_tau.UD)), table_tau_SD.UD, rep(NA, ncol(table_tau.UD)), table_tau_HPD0.UD), 2)
round(rbind(table_tau.MD6, rep(NA, ncol(table_tau.MD6)), table_tau_SD.MD6, rep(NA, ncol(table_tau.MD6)), table_tau_HPD0.MD6), 2)
round(rbind(table_tau.MD2, rep(NA, ncol(table_tau.MD2)), table_tau_SD.MD2, rep(NA, ncol(table_tau.MD2)), table_tau_HPD0.MD2), 2)
round(rbind(table_tau.BF6, rep(NA, ncol(table_tau.BF6)), table_tau_SD.BF6, rep(NA, ncol(table_tau.BF6)), table_tau_HPD0.BF6), 2)
round(rbind(table_tau.BF2, rep(NA, ncol(table_tau.BF2)), table_tau_SD.BF2, rep(NA, ncol(table_tau.BF2)), table_tau_HPD0.BF2), 2)



## Group means (95%-CI)
# UD
table_ginter.UD <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.UD <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.UD <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.UD[gg, 1] <- results.UD$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.UD[gg, 1] <- results.UD$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.UD[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.UD$gdata[[gg]]$Samples.ginter[1, ], results.UD$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.95), 1, hasNotZero))
}
# MD6
table_ginter.MD6 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.MD6 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.MD6 <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.MD6[gg, 1] <- results.MD6$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.MD6[gg, 1] <- results.MD6$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.MD6[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.MD6$gdata[[gg]]$Samples.ginter[1, ], results.MD6$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.95), 1, hasNotZero))
}
# MD2
table_ginter.MD2 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.MD2 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.MD2 <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.MD2[gg, 1] <- results.MD2$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.MD2[gg, 1] <- results.MD2$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.MD2[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.MD2$gdata[[gg]]$Samples.ginter[1, ], results.MD2$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.95), 1, hasNotZero))
}
# BF6
table_ginter.BF6 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.BF6 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.BF6 <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.BF6[gg, 1] <- results.BF6$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.BF6[gg, 1] <- results.BF6$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.BF6[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.BF6$gdata[[gg]]$Samples.ginter[1, ], results.BF6$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.95), 1, hasNotZero))
}
# BF2
table_ginter.BF2 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.BF2 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.BF2 <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.BF2[gg, 1] <- results.BF2$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.BF2[gg, 1] <- results.BF2$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.BF2[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.BF2$gdata[[gg]]$Samples.ginter[1, ], results.BF2$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.95), 1, hasNotZero))
}

# Print results
round(rbind(table_ginter.UD, rep(NA, ncol(table_ginter.UD)), table_ginter_SD.UD, rep(NA, ncol(table_ginter.UD)), table_ginter_HPD0.UD), 2)
round(rbind(table_ginter.MD6, rep(NA, ncol(table_ginter.MD6)), table_ginter_SD.MD6, rep(NA, ncol(table_ginter.MD6)), table_ginter_HPD0.MD6), 2)
round(rbind(table_ginter.MD2, rep(NA, ncol(table_ginter.MD2)), table_ginter_SD.MD2, rep(NA, ncol(table_ginter.MD2)), table_ginter_HPD0.MD2), 2)
round(rbind(table_ginter.BF6, rep(NA, ncol(table_ginter.BF6)), table_ginter_SD.BF6, rep(NA, ncol(table_ginter.BF6)), table_ginter_HPD0.BF6), 2)
round(rbind(table_ginter.BF2, rep(NA, ncol(table_ginter.BF2)), table_ginter_SD.BF2, rep(NA, ncol(table_ginter.BF2)), table_ginter_HPD0.BF2), 2)



## Subtest means (95%-CI)
# UD
Nt <- ncol(results.UD$gdata[[1]]$est.tinter) / 2
table_tinter.UD <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.UD <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.UD <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.UD[gg, 1:Nt] <- results.UD$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.UD[gg, 1:Nt] <- results.UD$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.UD[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.UD$gdata[[gg]]$Samples.tinter[1, ], results.UD$gdata[[gg]]$Samples.tinter[2, ])), prob = 0.95), 1, hasNotZero))
}
# MD6
Nt <- ncol(results.MD6$gdata[[1]]$est.tinter) / 2
table_tinter.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.MD6[gg, 1:Nt] <- results.MD6$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.MD6[gg, 1:Nt] <- results.MD6$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.MD6[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.MD6$gdata[[gg]]$Samples.tinter), prob = 0.95), 1, hasNotZero))
}
# MD2
Nt <- ncol(results.MD2$gdata[[1]]$est.tinter) / 2
table_tinter.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.MD2[gg, 1:Nt] <- results.MD2$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.MD2[gg, 1:Nt] <- results.MD2$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.MD2[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.MD2$gdata[[gg]]$Samples.tinter), prob = 0.95), 1, hasNotZero))
}
# BF6
Nt <- ncol(results.BF6$gdata[[1]]$est.tinter) / 2
table_tinter.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.BF6[gg, 1:Nt] <- results.BF6$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.BF6[gg, 1:Nt] <- results.BF6$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.BF6[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.BF6$gdata[[gg]]$Samples.tinter), prob = 0.95), 1, hasNotZero))
}
# BF2
Nt <- ncol(results.BF2$gdata[[1]]$est.tinter) / 2
table_tinter.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.BF2[gg, 1:Nt] <- results.BF2$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.BF2[gg, 1:Nt] <- results.BF2$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.BF2[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.BF2$gdata[[gg]]$Samples.tinter), prob = 0.95), 1, hasNotZero))
}

# Print results
round(rbind(table_tinter.UD, rep(NA, ncol(table_tinter.UD)), table_tinter_SD.UD, rep(NA, ncol(table_tinter.UD)), table_tinter_HPD0.UD), 2)
round(rbind(table_tinter.MD6, rep(NA, ncol(table_tinter.MD6)), table_tinter_SD.MD6, rep(NA, ncol(table_tinter.MD6)), table_tinter_HPD0.MD6), 2)
round(rbind(table_tinter.MD2, rep(NA, ncol(table_tinter.MD2)), table_tinter_SD.MD2, rep(NA, ncol(table_tinter.MD2)), table_tinter_HPD0.MD2), 2)
round(rbind(table_tinter.BF6, rep(NA, ncol(table_tinter.BF6)), table_tinter_SD.BF6, rep(NA, ncol(table_tinter.BF6)), table_tinter_HPD0.BF6), 2)
round(rbind(table_tinter.BF2, rep(NA, ncol(table_tinter.BF2)), table_tinter_SD.BF2, rep(NA, ncol(table_tinter.BF2)), table_tinter_HPD0.BF2), 2)



## Group means (90%-CI)
# UD
table_ginter.UD <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.UD <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.UD <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.UD[gg, 1] <- results.UD$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.UD[gg, 1] <- results.UD$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.UD[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.UD$gdata[[gg]]$Samples.ginter[1, ], results.UD$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.90), 1, hasNotZero))
}
# MD6
table_ginter.MD6 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.MD6 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.MD6 <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.MD6[gg, 1] <- results.MD6$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.MD6[gg, 1] <- results.MD6$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.MD6[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.MD6$gdata[[gg]]$Samples.ginter[1, ], results.MD6$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.90), 1, hasNotZero))
}
# MD2
table_ginter.MD2 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.MD2 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.MD2 <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.MD2[gg, 1] <- results.MD2$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.MD2[gg, 1] <- results.MD2$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.MD2[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.MD2$gdata[[gg]]$Samples.ginter[1, ], results.MD2$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.90), 1, hasNotZero))
}
# BF6
table_ginter.BF6 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.BF6 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.BF6 <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.BF6[gg, 1] <- results.BF6$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.BF6[gg, 1] <- results.BF6$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.BF6[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.BF6$gdata[[gg]]$Samples.ginter[1, ], results.BF6$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.90), 1, hasNotZero))
}
# BF2
table_ginter.BF2 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_SD.BF2 <- matrix(NA, nrow = 4, ncol = 1)
table_ginter_HPD0.BF2 <- matrix(NA, nrow = 4, ncol = 1)
for (gg in 1:4) {
  table_ginter.BF2[gg, 1] <- results.BF2$gdata[[gg]]$est.ginter[,1]
  table_ginter_SD.BF2[gg, 1] <- results.BF2$gdata[[gg]]$est.ginter[,2]
  table_ginter_HPD0.BF2[gg, 1] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.BF2$gdata[[gg]]$Samples.ginter[1, ], results.BF2$gdata[[gg]]$Samples.ginter[2, ])), prob = 0.90), 1, hasNotZero))
}

# Print results
round(rbind(table_ginter.UD, rep(NA, ncol(table_ginter.UD)), table_ginter_SD.UD, rep(NA, ncol(table_ginter.UD)), table_ginter_HPD0.UD), 2)
round(rbind(table_ginter.MD6, rep(NA, ncol(table_ginter.MD6)), table_ginter_SD.MD6, rep(NA, ncol(table_ginter.MD6)), table_ginter_HPD0.MD6), 2)
round(rbind(table_ginter.MD2, rep(NA, ncol(table_ginter.MD2)), table_ginter_SD.MD2, rep(NA, ncol(table_ginter.MD2)), table_ginter_HPD0.MD2), 2)
round(rbind(table_ginter.BF6, rep(NA, ncol(table_ginter.BF6)), table_ginter_SD.BF6, rep(NA, ncol(table_ginter.BF6)), table_ginter_HPD0.BF6), 2)
round(rbind(table_ginter.BF2, rep(NA, ncol(table_ginter.BF2)), table_ginter_SD.BF2, rep(NA, ncol(table_ginter.BF2)), table_ginter_HPD0.BF2), 2)



## Subtest means (90%-CI)
# UD
Nt <- ncol(results.UD$gdata[[1]]$est.tinter) / 2
table_tinter.UD <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.UD <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.UD <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.UD[gg, 1:Nt] <- results.UD$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.UD[gg, 1:Nt] <- results.UD$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.UD[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(c(results.UD$gdata[[gg]]$Samples.tinter[1, ], results.UD$gdata[[gg]]$Samples.tinter[2, ])), prob = 0.90), 1, hasNotZero))
}
# MD6
Nt <- ncol(results.MD6$gdata[[1]]$est.tinter) / 2
table_tinter.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.MD6 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.MD6[gg, 1:Nt] <- results.MD6$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.MD6[gg, 1:Nt] <- results.MD6$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.MD6[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.MD6$gdata[[gg]]$Samples.tinter), prob = 0.90), 1, hasNotZero))
}
# MD2
Nt <- ncol(results.MD2$gdata[[1]]$est.tinter) / 2
table_tinter.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.MD2 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.MD2[gg, 1:Nt] <- results.MD2$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.MD2[gg, 1:Nt] <- results.MD2$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.MD2[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.MD2$gdata[[gg]]$Samples.tinter), prob = 0.90), 1, hasNotZero))
}
# BF6
Nt <- ncol(results.BF6$gdata[[1]]$est.tinter) / 2
table_tinter.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.BF6 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.BF6[gg, 1:Nt] <- results.BF6$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.BF6[gg, 1:Nt] <- results.BF6$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.BF6[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.BF6$gdata[[gg]]$Samples.tinter), prob = 0.90), 1, hasNotZero))
}
# BF2
Nt <- ncol(results.BF2$gdata[[1]]$est.tinter) / 2
table_tinter.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_SD.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
table_tinter_HPD0.BF2 <- matrix(NA, nrow = 4, ncol = Nt)
for (gg in 1:4) {
  table_tinter.BF2[gg, 1:Nt] <- results.BF2$gdata[[gg]]$est.tinter[,1:Nt]
  table_tinter_SD.BF2[gg, 1:Nt] <- results.BF2$gdata[[gg]]$est.tinter[,(Nt+1):(2*Nt)]
  table_tinter_HPD0.BF2[gg, 1:Nt] <- as.numeric(apply(coda::HPDinterval(coda::mcmc(results.BF2$gdata[[gg]]$Samples.tinter), prob = 0.90), 1, hasNotZero))
}

# Print results
round(rbind(table_tinter.UD, rep(NA, ncol(table_tinter.UD)), table_tinter_SD.UD, rep(NA, ncol(table_tinter.UD)), table_tinter_HPD0.UD), 2)
round(rbind(table_tinter.MD6, rep(NA, ncol(table_tinter.MD6)), table_tinter_SD.MD6, rep(NA, ncol(table_tinter.MD6)), table_tinter_HPD0.MD6), 2)
round(rbind(table_tinter.MD2, rep(NA, ncol(table_tinter.MD2)), table_tinter_SD.MD2, rep(NA, ncol(table_tinter.MD2)), table_tinter_HPD0.MD2), 2)
round(rbind(table_tinter.BF6, rep(NA, ncol(table_tinter.BF6)), table_tinter_SD.BF6, rep(NA, ncol(table_tinter.BF6)), table_tinter_HPD0.BF6), 2)
round(rbind(table_tinter.BF2, rep(NA, ncol(table_tinter.BF2)), table_tinter_SD.BF2, rep(NA, ncol(table_tinter.BF2)), table_tinter_HPD0.BF2), 2)


## Effect sizes
# UD
for (gg in 1:3) {
  print(round(table_ginter.UD[gg + 1] - table_ginter.UD[gg], 2))  
}
# MD6
for (gg in 1:3) {
  print(round(table_ginter.MD6[gg + 1] - table_ginter.MD6[gg], 2))  
}
for (gg in 1:3) {
  print(round((table_tinter.MD6[gg + 1, ] + table_ginter.MD6[gg + 1]) - (table_tinter.MD6[gg, ] + table_ginter.MD6[gg]), 2))  
}
# MD2
for (gg in 1:3) {
  print(round(table_ginter.MD2[gg + 1] - table_ginter.MD2[gg], 2))  
}
for (gg in 1:3) {
  print(round((table_tinter.MD2[gg + 1, ] + table_ginter.MD2[gg + 1]) - (table_tinter.MD2[gg, ] + table_ginter.MD2[gg]), 2))  
}
# BF6
for (gg in 1:3) {
  print(round(table_ginter.BF6[gg + 1] - table_ginter.BF6[gg], 2))  
}
for (gg in 1:3) {
  print(round((table_tinter.BF6[gg + 1, ] + table_ginter.BF6[gg + 1]) - (table_tinter.BF6[gg, ] + table_ginter.BF6[gg]), 2))  
}
# BF2
for (gg in 1:3) {
  print(round(table_ginter.BF2[gg + 1] - table_ginter.BF2[gg], 2))  
}
for (gg in 1:3) {
  print(round((table_tinter.BF2[gg + 1, ] + table_ginter.BF2[gg + 1]) - (table_tinter.BF2[gg, ] + table_ginter.BF2[gg]), 2))  
}

## Marginal likelihood 
table_results <- matrix(NA, nrow = 5, ncol = 1)
rownames(table_results) <- c("UD", "MD6", "MD2", "BF6", "BF2")

table_results[1, ] <- sum(results.UD$marllik)
table_results[2, ] <- sum(results.MD6$marllik)
table_results[3, ] <- sum(results.MD2$marllik)
table_results[4, ] <- sum(results.BF6$marllik)
table_results[5, ] <- sum(results.BF2$marllik)

table_names <- matrix(rownames(table_results), nrow = 5, ncol = 1)
table_sorted <- matrix(table_names[apply(table_results, 2, order, decreasing = TRUE)], nrow = 5, ncol = 1)
table_results
table_sorted # Sorted by log marginal likelihood per group (decreasing order)
