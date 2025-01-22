library(Cairo)
library(lattice)

source("Tools.R")

############ Table 'Estimating dependence structure' ############
sim.values <- seq(from = -0.05, to = 0.10, by = 0.01)

## N = 500
filenames <- paste0("../Results/Run_", c("m005", "m004", "m003", "m002", "m001", "000", "001", "002", "003", "004", "005", "006", "007", "008", "009", "010"), "_set1.RDS")
Est.500 <- matrix(NA, nrow = 16, ncol = 25)
colnames(Est.500) <- c("Sim", "SDM.Est.f1", "SDM.Est.f3_1", "SDM.Est.f3_2", "SDM.Est.f3_3", "SDM.SD.f1", "SDM.SD.f3_1", "SDM.SD.f3_2", "SDM.SD.f3_3", "SDM.RMSE.f1", "SDM.RMSE.f3_1", "SDM.RMSE.f3_2", "SDM.RMSE.f3_3", 
                       "CFA.Est.f1", "CFA.Est.f3_1", "CFA.Est.f3_2", "CFA.Est.f3_3", "CFA.SD.f1", "CFA.SD.f3_1", "CFA.SD.f3_2", "CFA.SD.f3_3", "CFA.RMSE.f1", "CFA.RMSE.f3_1", "CFA.RMSE.f3_2", "CFA.RMSE.f3_3")
for (ii in 1:length(sim.values)) {
  results.BCSM <- readRDS(filenames[ii])$results.BCSM$est.M1
  results.lavaan <- readRDS(filenames[ii])$results.lavaan$lavaan.est.M1
  Est.500[ii, 1] <- sim.values[ii]
  Est.500[ii, 2] <- mean(results.BCSM[, 4]) - 0.3
  Est.500[ii, 3:5] <- colMeans(results.BCSM[, 1:3]) - sim.values[ii]
  Est.500[ii, 6] <- sd(results.BCSM[, 4])
  Est.500[ii, 7:9] <- colSDs(results.BCSM[, 1:3])
  Est.500[ii, 10] <- sqrt(mean((results.BCSM[, 4] - 0.3)^2)) 
  Est.500[ii, 11:13] <- sqrt(colMeans((results.BCSM[, 1:3] - sim.values[ii])^2))
  Est.500[ii, 14:17] <- colMeans(results.lavaan[, 1:4]) - c(0.3, rep(sim.values[ii], 3))
  Est.500[ii, 18:21] <- colSDs(results.lavaan[, 1:4])
  Est.500[ii, 22] <- sqrt(mean((results.lavaan[, 1] - 0.3)^2))
  Est.500[ii, 23:25] <- sqrt(colMeans((results.lavaan[, 2:4] - sim.values[ii])^2))
}

## N = 300
filenames <- paste0("../Results/Run_", c("m005", "m004", "m003", "m002", "m001", "000", "001", "002", "003", "004", "005", "006", "007", "008", "009", "010"), "_set1_N300.RDS")
Est.300 <- matrix(NA, nrow = 16, ncol = 25)
colnames(Est.300) <- c("Sim", "SDM.Est.f1", "SDM.Est.f3_1", "SDM.Est.f3_2", "SDM.Est.f3_3", "SDM.SD.f1", "SDM.SD.f3_1", "SDM.SD.f3_2", "SDM.SD.f3_3", "SDM.RMSE.f1", "SDM.RMSE.f3_1", "SDM.RMSE.f3_2", "SDM.RMSE.f3_3", 
                       "CFA.Est.f1", "CFA.Est.f3_1", "CFA.Est.f3_2", "CFA.Est.f3_3", "CFA.SD.f1", "CFA.SD.f3_1", "CFA.SD.f3_2", "CFA.SD.f3_3", "CFA.RMSE.f1", "CFA.RMSE.f3_1", "CFA.RMSE.f3_2", "CFA.RMSE.f3_3")
for (ii in 1:length(sim.values)) {
  results.BCSM <- readRDS(filenames[ii])$results.BCSM$est.M1
  results.lavaan <- readRDS(filenames[ii])$results.lavaan$lavaan.est.M1
  Est.300[ii, 1] <- sim.values[ii]
  Est.300[ii, 2] <- mean(results.BCSM[, 4]) - 0.3
  Est.300[ii, 3:5] <- colMeans(results.BCSM[, 1:3]) - sim.values[ii]
  Est.300[ii, 6] <- sd(results.BCSM[, 4])
  Est.300[ii, 7:9] <- colSDs(results.BCSM[, 1:3])
  Est.300[ii, 10] <- sqrt(mean((results.BCSM[, 4] - 0.3)^2))
  Est.300[ii, 11:13] <- sqrt(colMeans((results.BCSM[, 1:3] - sim.values[ii])^2))
  Est.300[ii, 14:17] <- colMeans(results.lavaan[, 1:4]) - c(0.3, rep(sim.values[ii], 3))
  Est.300[ii, 18:21] <- colSDs(results.lavaan[, 1:4])
  Est.300[ii, 22] <- sqrt(mean((results.lavaan[, 1] - 0.3)^2))
  Est.300[ii, 23:25] <- sqrt(colMeans((results.lavaan[, 2:4] - sim.values[ii])^2))
}


## N = 100
filenames <- paste0("../Results/Run_", c("m005", "m004", "m003", "m002", "m001", "000", "001", "002", "003", "004", "005", "006", "007", "008", "009", "010"), "_set1_N100.RDS")
Est.100 <- matrix(NA, nrow = 16, ncol = 25)
colnames(Est.100) <- c("Sim", "SDM.Est.f1", "SDM.Est.f3_1", "SDM.Est.f3_2", "SDM.Est.f3_3", "SDM.SD.f1", "SDM.SD.f3_1", "SDM.SD.f3_2", "SDM.SD.f3_3", "SDM.RMSE.f1", "SDM.RMSE.f3_1", "SDM.RMSE.f3_2", "SDM.RMSE.f3_3", 
                       "CFA.Est.f1", "CFA.Est.f3_1", "CFA.Est.f3_2", "CFA.Est.f3_3", "CFA.SD.f1", "CFA.SD.f3_1", "CFA.SD.f3_2", "CFA.SD.f3_3", "CFA.RMSE.f1", "CFA.RMSE.f3_1", "CFA.RMSE.f3_2", "CFA.RMSE.f3_3")
for (ii in 1:length(sim.values)) {
  results.BCSM <- readRDS(filenames[ii])$results.BCSM$est.M1
  results.lavaan <- readRDS(filenames[ii])$results.lavaan$lavaan.est.M1
  Est.100[ii, 1] <- sim.values[ii]
  Est.100[ii, 2] <- mean(results.BCSM[, 4]) - 0.3
  Est.100[ii, 3:5] <- colMeans(results.BCSM[, 1:3]) - sim.values[ii]
  Est.100[ii, 6] <- sd(results.BCSM[, 4])
  Est.100[ii, 7:9] <- colSDs(results.BCSM[, 1:3])
  Est.100[ii, 10] <- sqrt(mean((results.BCSM[, 4] - 0.3)^2))
  Est.100[ii, 11:13] <- sqrt(colMeans((results.BCSM[, 1:3] - sim.values[ii])^2))
  Est.100[ii, 14:17] <- colMeans(results.lavaan[, 1:4])- c(0.3, rep(sim.values[ii], 3))
  Est.100[ii, 18:21] <- colSDs(results.lavaan[, 1:4])
  Est.100[ii, 22] <- sqrt(mean((results.lavaan[, 1] - 0.3)^2))
  Est.100[ii, 23:25] <- sqrt(colMeans((results.lavaan[, 2:4] - sim.values[ii])^2))
}

round(Est.500, 2)
round(Est.300, 2)
round(Est.100, 2)


# Sim, SDM.Bias, SDM.RMSE, Lavaan.Bias, Lavaan.RMSE



############ Plot 'Model selection' ############
## N = 500
# BCSM
summary.BCSM <- numeric(0)
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m005_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m004_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m003_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m002_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m001_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_000_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_001_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_002_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_003_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_004_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_005_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_006_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_007_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_008_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_009_set1.RDS")$summary.BCSM[-1])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_010_set1.RDS")$summary.BCSM[-1])
# Lavaan
summary.lavaan <- numeric(0)
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m005_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m004_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m003_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m002_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m001_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_000_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_001_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_002_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_003_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_004_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_005_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_006_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_007_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_008_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_009_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_010_set1.RDS")$summary.lavaan[-1])

## N = 300
# BCSM
summary.BCSM.N300 <- numeric(0)
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_m005_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_m004_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_m003_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_m002_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_m001_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_000_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_001_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_002_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_003_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_004_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_005_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_006_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_007_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_008_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_009_set1_N300.RDS")$summary.BCSM[-1])
summary.BCSM.N300 <- rbind(summary.BCSM.N300,readRDS("../Results/Run_010_set1_N300.RDS")$summary.BCSM[-1])
# Lavaan
summary.lavaan.N300 <- numeric(0)
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_m005_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_m004_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_m003_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_m002_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_m001_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_000_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_001_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_002_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_003_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_004_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_005_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_006_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_007_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_008_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_009_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan.N300 <- rbind(summary.lavaan.N300,readRDS("../Results/Run_010_set1_N300.RDS")$summary.lavaan[-1])

## N = 100
# BCSM
summary.BCSM.N100 <- numeric(0)
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_m005_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_m004_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_m003_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_m002_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_m001_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_000_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_001_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_002_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_003_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_004_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_005_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_006_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_007_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_008_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_009_set1_N100.RDS")$summary.BCSM[-1])
summary.BCSM.N100 <- rbind(summary.BCSM.N100,readRDS("../Results/Run_010_set1_N100.RDS")$summary.BCSM[-1])
# Lavaan
summary.lavaan.N100 <- numeric(0)
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_m005_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_m004_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_m003_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_m002_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_m001_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_000_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_001_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_002_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_003_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_004_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_005_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_006_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_007_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_008_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_009_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan.N100 <- rbind(summary.lavaan.N100,readRDS("../Results/Run_010_set1_N100.RDS")$summary.lavaan[-1])


opar <- par() 
#CairoWin()

cairo_pdf("Sim_study_model_selection.pdf", width = 14, height = 20, pointsize = 14)
par(opar)
par(mfrow=c(3,2))
par(mar=c(10,4.3,2,2.1))
par(oma=c(0,0,2,0))
## Plot BF vs frequentist fit indices
lnames <- c(expression("BF"["tot"]), expression("BF"["ESS"]), expression("RMSEA"[Delta]), expression("RMSEA"[".05"]), expression("RMSEA"[".10"]), expression("TLI"[".95"]))

## N = 500
# Select Model M0
col <- c("grey20", "grey20", "black", "grey20", "black", "grey20", "black", "black", "black")
col.lines <- c("grey60", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey20", "grey20")
pch <- c(23, 23, 23, 3, 3, 24, 24, 15, 16) # 1: RMSEA (0.05), 2: RMSEA (0.10), 3: RMSEA (delta), 4: SRMR (0.08), 6: TLI (0.95), 8: BF (total), 9: BF (ESS)
lty <- c(3, 2, 1, 2, 1, 2, 1, 1, 1)
cex <- 1.6
cex.main <- 2
cex.lab = 2
cex.axis = 1.3
lwd <- 2
colnames(summary.BCSM)[c(9,7)]
colnames(summary.lavaan)[c(2,4,6,8,12)]
plot(summary.lavaan[, 1], summary.lavaan[, 2], type = "l", pch = pch[1], lty = lty[1], col = col.lines[1], ylim = c(0, 100), cex = cex, lwd = lwd, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, xaxt = "n", main = "Select Model M0", xlab = "Conditional within-cluster correlation", ylab = "% M0 Selected")
axis(1, at = summary.lavaan[, 1], labels = c("-.05", "-.04", "-.03", "-.02", "-.01", ".0", ".01", ".02", ".03", ".04", ".05", ".06", ".07", ".08", ".09", ".10"), cex.axis = 1.3)
lines(summary.lavaan[, 1], summary.lavaan[, 4], lty = lty[2], col = col.lines[2], lwd = lwd)
lines(summary.lavaan[, 1], summary.lavaan[, 6], lty = lty[3], col = col.lines[3], lwd = lwd)
#lines(summary.lavaan[, 1], summary.lavaan[, 8], lty = lty[4], col = col.lines[4], lwd = lwd)
#lines(summary.lavaan[, 1], summary.lavaan[, 10], lty = lty[5], col = col.lines[5], lwd = lwd)
lines(summary.lavaan[, 1], summary.lavaan[, 12], lty = lty[6], col = col.lines[6], lwd = lwd)
#lines(summary.lavaan[, 1], summary.lavaan[, 14], lty = lty[7], col = col.lines[7], lwd = lwd)
lines(summary.lavaan[, 1], summary.BCSM[, 9], lty = lty[8], col = col.lines[8], lwd = lwd)
lines(summary.lavaan[, 1], summary.BCSM[, 7], lty = lty[9], col = col.lines[9], lwd = lwd)
points(summary.lavaan[, 1], summary.lavaan[, 2], pch = pch[1], col = col[1], cex = cex, bg = "white")
points(summary.lavaan[, 1], summary.lavaan[, 4], pch = pch[2], col = col[2], cex = cex, bg = "white")
points(summary.lavaan[, 1], summary.lavaan[, 6], pch = pch[3], col = col[3], cex = cex, bg = "white")
#points(summary.lavaan[, 1], summary.lavaan[, 8], pch = pch[4], col = col[4], cex = cex)
#points(summary.lavaan[, 1], summary.lavaan[, 10], pch = pch[5], col = col[5], cex = cex)
points(summary.lavaan[, 1], summary.lavaan[, 12], pch = pch[6], col = col[6], cex = cex, bg = "white")
#points(summary.lavaan[, 1], summary.lavaan[, 14], pch = pch[7], col = col[7], cex = cex)
points(summary.lavaan[, 1], summary.BCSM[, 9], pch = pch[8], col = col[8], cex = cex)
points(summary.lavaan[, 1], summary.BCSM[, 7], pch = pch[9], col = col[9], cex = cex)

# Select Model M1
col <- c("grey20", "grey20", "black", "grey20", "black", "grey20", "black", "black", "black")
col.lines <- c("grey70", "grey70", "grey40", "grey70", "grey40", "grey70", "grey40", "grey20", "grey20")
pch <- c(23, 23, 23, 3, 3, 24, 24, 15, 16) # 1: RMSEA (0.05), 2: RMSEA (0.10), 3: RMSEA (delta), 4: SRMR (0.08), 6: TLI (0.95), 8: BF (total), 9: BF (ESS)
lty <- c(3, 2, 1, 2, 1, 2, 1, 1, 1)
cex <- 1.5
lwd <- 2
colnames(summary.BCSM)[c(8,6)]
colnames(summary.lavaan)[c(3,5,7,13)]
plot(summary.lavaan[, 1], summary.lavaan[, 3], type = "l", pch = pch[1], lty = lty[1], col = col.lines[1], ylim = c(0, 100), cex = cex, lwd = lwd, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, xaxt = "n", main = "Select Model M1", xlab = "Conditional within-cluster correlation", ylab = "% M1 Selected")
axis(1, at = summary.lavaan[, 1], labels = c("-.05", "-.04", "-.03", "-.02", "-.01", ".0", ".01", ".02", ".03", ".04", ".05", ".06", ".07", ".08", ".09", ".10"), cex.axis = 1.3)
lines(summary.lavaan[, 1], summary.lavaan[, 5], lty = lty[2], col = col.lines[2], lwd = lwd)
lines(summary.lavaan[, 1], summary.lavaan[, 7], lty = lty[3], col = col.lines[3], lwd = lwd)
#lines(summary.lavaan[, 1], summary.lavaan[, 9], lty = lty[4], col = col.lines[4], lwd = lwd)
#lines(summary.lavaan[, 1], summary.lavaan[, 10], lty = lty[5], col = col.lines[5], lwd = lwd)
lines(summary.lavaan[, 1], summary.lavaan[, 13], lty = lty[6], col = col.lines[6], lwd = lwd)
#lines(summary.lavaan[, 1], summary.lavaan[, 14], lty = lty[7], col = col.lines[7], lwd = lwd)
lines(summary.lavaan[, 1], summary.BCSM[, 8], lty = lty[8], col = col.lines[8], lwd = lwd)
lines(summary.lavaan[, 1], summary.BCSM[, 6], lty = lty[9], col = col.lines[9], lwd = lwd)
points(summary.lavaan[, 1], summary.lavaan[, 3], pch = pch[1], col = col[1], cex = cex, bg = "white")
points(summary.lavaan[, 1], summary.lavaan[, 5], pch = pch[2], col = col[2], cex = cex, bg = "white")
points(summary.lavaan[, 1], summary.lavaan[, 7], pch = pch[3], col = col[3], cex = cex, bg = "white")
#points(summary.lavaan[, 1], summary.lavaan[, 9], pch = pch[4], col = col[4], cex = cex)
#points(summary.lavaan[, 1], summary.lavaan[, 10], pch = pch[5], col = col[5], cex = cex)
points(summary.lavaan[, 1], summary.lavaan[, 13], pch = pch[6], col = col[6], cex = cex, bg = "white")
#points(summary.lavaan[, 1], summary.lavaan[, 14], pch = pch[7], col = col[7], cex = cex)
points(summary.lavaan[, 1], summary.BCSM[, 8], pch = pch[8], col = col[8], cex = cex)
points(summary.lavaan[, 1], summary.BCSM[, 6], pch = pch[9], col = col[9], cex = cex)

## N = 300
# Select Model M0
plot(summary.lavaan.N300[, 1], summary.lavaan.N300[, 2], type = "l", pch = pch[1], lty = lty[1], col = col.lines[1], ylim = c(0, 100), cex = cex, lwd = lwd, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, xaxt = "n", main = "Select Model M0", xlab = "Conditional within-cluster correlation", ylab = "% M0 Selected")
axis(1, at = summary.lavaan.N300[, 1], labels = c("-.05", "-.04", "-.03", "-.02", "-.01", ".0", ".01", ".02", ".03", ".04", ".05", ".06", ".07", ".08", ".09", ".10"), cex.axis = 1.3)
lines(summary.lavaan.N300[, 1], summary.lavaan.N300[, 4], lty = lty[2], col = col.lines[2], lwd = lwd)
lines(summary.lavaan.N300[, 1], summary.lavaan.N300[, 6], lty = lty[3], col = col.lines[3], lwd = lwd)
lines(summary.lavaan.N300[, 1], summary.lavaan.N300[, 12], lty = lty[6], col = col.lines[6], lwd = lwd)
lines(summary.lavaan.N300[, 1], summary.BCSM.N300[, 9], lty = lty[8], col = col.lines[8], lwd = lwd)
lines(summary.lavaan.N300[, 1], summary.BCSM.N300[, 7], lty = lty[9], col = col.lines[9], lwd = lwd)
points(summary.lavaan.N300[, 1], summary.lavaan.N300[, 2], pch = pch[1], col = col[1], cex = cex, bg = "white")
points(summary.lavaan.N300[, 1], summary.lavaan.N300[, 4], pch = pch[2], col = col[2], cex = cex, bg = "white")
points(summary.lavaan.N300[, 1], summary.lavaan.N300[, 6], pch = pch[3], col = col[3], cex = cex, bg = "white")
points(summary.lavaan.N300[, 1], summary.lavaan.N300[, 12], pch = pch[6], col = col[6], cex = cex, bg = "white")
points(summary.lavaan.N300[, 1], summary.BCSM.N300[, 9], pch = pch[8], col = col[8], cex = cex)
points(summary.lavaan.N300[, 1], summary.BCSM.N300[, 7], pch = pch[9], col = col[9], cex = cex)

# Select Model M1
plot(summary.lavaan.N300[, 1], summary.lavaan.N300[, 3], type = "l", pch = pch[1], lty = lty[1], col = col.lines[1], ylim = c(0, 100), cex = cex, lwd = lwd, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, xaxt = "n", main = "Select Model M1", xlab = "Conditional within-cluster correlation", ylab = "% M1 Selected")
axis(1, at = summary.lavaan.N300[, 1], labels = c("-.05", "-.04", "-.03", "-.02", "-.01", ".0", ".01", ".02", ".03", ".04", ".05", ".06", ".07", ".08", ".09", ".10"), cex.axis = 1.3)
lines(summary.lavaan.N300[, 1], summary.lavaan.N300[, 5], lty = lty[2], col = col.lines[2], lwd = lwd)
lines(summary.lavaan.N300[, 1], summary.lavaan.N300[, 7], lty = lty[3], col = col.lines[3], lwd = lwd)
lines(summary.lavaan.N300[, 1], summary.lavaan.N300[, 13], lty = lty[6], col = col.lines[6], lwd = lwd)
lines(summary.lavaan.N300[, 1], summary.BCSM.N300[, 8], lty = lty[8], col = col.lines[8], lwd = lwd)
lines(summary.lavaan.N300[, 1], summary.BCSM.N300[, 6], lty = lty[9], col = col.lines[9], lwd = lwd)
points(summary.lavaan.N300[, 1], summary.lavaan.N300[, 3], pch = pch[1], col = col[1], cex = cex, bg = "white")
points(summary.lavaan.N300[, 1], summary.lavaan.N300[, 5], pch = pch[2], col = col[2], cex = cex, bg = "white")
points(summary.lavaan.N300[, 1], summary.lavaan.N300[, 7], pch = pch[3], col = col[3], cex = cex, bg = "white")
points(summary.lavaan.N300[, 1], summary.lavaan.N300[, 13], pch = pch[6], col = col[6], cex = cex, bg = "white")
points(summary.lavaan.N300[, 1], summary.BCSM.N300[, 8], pch = pch[8], col = col[8], cex = cex)
points(summary.lavaan.N300[, 1], summary.BCSM.N300[, 6], pch = pch[9], col = col[9], cex = cex)

## N = 100
# Select Model M0
plot(summary.lavaan.N100[, 1], summary.lavaan.N100[, 2], type = "l", pch = pch[1], lty = lty[1], col = col.lines[1], ylim = c(0, 100), cex = cex, lwd = lwd, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, xaxt = "n", main = "Select Model M0", xlab = "Conditional within-cluster correlation", ylab = "% M0 Selected")
axis(1, at = summary.lavaan.N100[, 1], labels = c("-.05", "-.04", "-.03", "-.02", "-.01", ".0", ".01", ".02", ".03", ".04", ".05", ".06", ".07", ".08", ".09", ".10"), cex.axis = 1.3)
lines(summary.lavaan.N100[, 1], summary.lavaan.N100[, 4], lty = lty[2], col = col.lines[2], lwd = lwd)
lines(summary.lavaan.N100[, 1], summary.lavaan.N100[, 6], lty = lty[3], col = col.lines[3], lwd = lwd)
lines(summary.lavaan.N100[, 1], summary.lavaan.N100[, 12], lty = lty[6], col = col.lines[6], lwd = lwd)
lines(summary.lavaan.N100[, 1], summary.BCSM.N100[, 9], lty = lty[8], col = col.lines[8], lwd = lwd)
lines(summary.lavaan.N100[, 1], summary.BCSM.N100[, 7], lty = lty[9], col = col.lines[9], lwd = lwd)
points(summary.lavaan.N100[, 1], summary.lavaan.N100[, 2], pch = pch[1], col = col[1], cex = cex, bg = "white")
points(summary.lavaan.N100[, 1], summary.lavaan.N100[, 4], pch = pch[2], col = col[2], cex = cex, bg = "white")
points(summary.lavaan.N100[, 1], summary.lavaan.N100[, 6], pch = pch[3], col = col[3], cex = cex, bg = "white")
points(summary.lavaan.N100[, 1], summary.lavaan.N100[, 12], pch = pch[6], col = col[6], cex = cex, bg = "white")
points(summary.lavaan.N100[, 1], summary.BCSM.N100[, 9], pch = pch[8], col = col[8], cex = cex)
points(summary.lavaan.N100[, 1], summary.BCSM.N100[, 7], pch = pch[9], col = col[9], cex = cex)

# Select Model M1
plot(summary.lavaan.N100[, 1], summary.lavaan.N100[, 3], type = "l", pch = pch[1], lty = lty[1], col = col.lines[1], ylim = c(0, 100), cex = cex, lwd = lwd, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, xaxt = "n", main = "Select Model M1", xlab = "Conditional within-cluster correlation", ylab = "% M1 Selected")
axis(1, at = summary.lavaan.N100[, 1], labels = c("-.05", "-.04", "-.03", "-.02", "-.01", ".0", ".01", ".02", ".03", ".04", ".05", ".06", ".07", ".08", ".09", ".10"), cex.axis = 1.3)
lines(summary.lavaan.N100[, 1], summary.lavaan.N100[, 5], lty = lty[2], col = col.lines[2], lwd = lwd)
lines(summary.lavaan.N100[, 1], summary.lavaan.N100[, 7], lty = lty[3], col = col.lines[3], lwd = lwd)
lines(summary.lavaan.N100[, 1], summary.lavaan.N100[, 13], lty = lty[6], col = col.lines[6], lwd = lwd)
lines(summary.lavaan.N100[, 1], summary.BCSM.N100[, 8], lty = lty[8], col = col.lines[8], lwd = lwd)
lines(summary.lavaan.N100[, 1], summary.BCSM.N100[, 6], lty = lty[9], col = col.lines[9], lwd = lwd)
points(summary.lavaan.N100[, 1], summary.lavaan.N100[, 3], pch = pch[1], col = col[1], cex = cex, bg = "white")
points(summary.lavaan.N100[, 1], summary.lavaan.N100[, 5], pch = pch[2], col = col[2], cex = cex, bg = "white")
points(summary.lavaan.N100[, 1], summary.lavaan.N100[, 7], pch = pch[3], col = col[3], cex = cex, bg = "white")
points(summary.lavaan.N100[, 1], summary.lavaan.N100[, 13], pch = pch[6], col = col[6], cex = cex, bg = "white")
points(summary.lavaan.N100[, 1], summary.BCSM.N100[, 8], pch = pch[8], col = col[8], cex = cex)
points(summary.lavaan.N100[, 1], summary.BCSM.N100[, 6], pch = pch[9], col = col[9], cex = cex)

# Title
mtext("N = 500", outer = TRUE, cex = 1.45, font = 2, line = 0)#, at = 0.05)
mtext("N = 300", outer = TRUE, cex = 1.45, font = 2, line = -42.5)#, at = 0.05)
mtext("N = 100", outer = TRUE, cex = 1.45, font = 2, line = -85)#, at = 0.05)

# Legend
legend(x = -0.108, y = -14.5, horiz = FALSE, ncol = 2, legend = lnames, pch = pch[c(8, 9, 3, 1, 2, 6)], lty = lty[c(8, 9, 3, 1, 2, 6)], pt.bg = "white", bty = "y", cex = 1.6, lwd = 2, xpd = NA)
dev.off()



############ Table 'Model selection' ############
which.cat <- c("000", "001", "002", "003", "004", "005", "006", "007")

## N = 500
filenames <- paste0("../Results/Run_", which.cat, "_set1.RDS")
MS.500 <- matrix(NA, ncol = length(filenames), nrow = 4)
colnames(MS.500) <- which.cat
for (ii in 1:length(filenames)) {
  tmp <- readRDS(filenames[ii])
  cat0 <- tmp$results.BCSM$logBF
  MS.500[1, ii] <- 1 - length(which(cat0 < 0)) / length(cat0)# * 100
  MS.500[2, ii] <- 1 - length(which(cat0 < 1)) / length(cat0)# * 100
  MS.500[3, ii] <- 1 - length(which(cat0 < 3)) / length(cat0)# * 100
  MS.500[4, ii] <- 1 - length(which(cat0 < 5)) / length(cat0)# * 100
}

## N = 300
filenames <- paste0("../Results/Run_", which.cat, "_set1_N300.RDS")
MS.300 <- matrix(NA, ncol = length(filenames), nrow = 4)
colnames(MS.300) <- which.cat
for (ii in 1:length(filenames)) {
  tmp <- readRDS(filenames[ii])
  cat0 <- tmp$results.BCSM$logBF
  MS.300[1, ii] <- 1 - length(which(cat0 < 0)) / length(cat0)# * 100
  MS.300[2, ii] <- 1 - length(which(cat0 < 1)) / length(cat0)# * 100
  MS.300[3, ii] <- 1 - length(which(cat0 < 3)) / length(cat0)# * 100
  MS.300[4, ii] <- 1 - length(which(cat0 < 5)) / length(cat0)# * 100
}

## N = 100
filenames <- paste0("../Results/Run_", which.cat, "_set1_N100.RDS")
MS.100 <- matrix(NA, ncol = length(filenames), nrow = 4)
colnames(MS.100) <- which.cat
for (ii in 1:length(filenames)) {
  tmp <- readRDS(filenames[ii])
  cat0 <- tmp$results.BCSM$logBF
  MS.100[1, ii] <- 1 - length(which(cat0 < 0)) / length(cat0)# * 100
  MS.100[2, ii] <- 1 - length(which(cat0 < 1)) / length(cat0)# * 100
  MS.100[3, ii] <- 1 - length(which(cat0 < 3)) / length(cat0)# * 100
  MS.100[4, ii] <- 1 - length(which(cat0 < 5)) / length(cat0)# * 100
}


MS.500
MS.300
MS.100


############ Table 'Model selection Appendix' ############
## N = 500
summary.BCSM <- numeric(0)
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m005_set1.RDS")$summary.BCSM[c(2,7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m004_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m003_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m002_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m001_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_000_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_001_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_002_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_003_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_004_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_005_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_006_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_007_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_008_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_009_set1.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_010_set1.RDS")$summary.BCSM[c(2, 7:10)])
# Lavaan
summary.lavaan <- numeric(0)
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m005_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m004_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m003_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m002_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m001_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_000_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_001_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_002_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_003_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_004_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_005_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_006_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_007_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_008_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_009_set1.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_010_set1.RDS")$summary.lavaan[-1])

summary.BCSM
summary.lavaan

## N = 300
summary.BCSM <- numeric(0)
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m005_set1_N300.RDS")$summary.BCSM[c(2,7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m004_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m003_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m002_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m001_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_000_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_001_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_002_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_003_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_004_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_005_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_006_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_007_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_008_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_009_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_010_set1_N300.RDS")$summary.BCSM[c(2, 7:10)])
# Lavaan
summary.lavaan <- numeric(0)
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m005_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m004_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m003_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m002_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m001_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_000_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_001_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_002_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_003_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_004_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_005_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_006_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_007_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_008_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_009_set1_N300.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_010_set1_N300.RDS")$summary.lavaan[-1])

summary.BCSM
summary.lavaan

## N = 100
summary.BCSM <- numeric(0)
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m005_set1_N100.RDS")$summary.BCSM[c(2,7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m004_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m003_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m002_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_m001_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_000_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_001_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_002_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_003_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_004_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_005_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_006_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_007_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_008_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_009_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
summary.BCSM <- rbind(summary.BCSM,readRDS("../Results/Run_010_set1_N100.RDS")$summary.BCSM[c(2, 7:10)])
# Lavaan
summary.lavaan <- numeric(0)
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m005_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m004_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m003_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m002_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_m001_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_000_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_001_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_002_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_003_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_004_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_005_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_006_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_007_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_008_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_009_set1_N100.RDS")$summary.lavaan[-1])
summary.lavaan <- rbind(summary.lavaan,readRDS("../Results/Run_010_set1_N100.RDS")$summary.lavaan[-1])

summary.BCSM
summary.lavaan




