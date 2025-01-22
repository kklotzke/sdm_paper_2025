rm(list = ls())

tau.f3 <- 0.04
tau.sim <- c(tau.f3 , tau.f3, tau.f3, 0.3)
rep <- 250

set.seed(2)
filename = "../Results/Run_004_set1.RDS"

source("SDM BF Testing - Core.R")
