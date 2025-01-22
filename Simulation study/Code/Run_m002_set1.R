rm(list = ls())

tau.f3 <- -0.02
tau.sim <- c(tau.f3 , tau.f3, tau.f3, 0.3)
rep <- 250

set.seed(2)
filename = "../Results/Run_m002_set1.RDS"

source("SDM BF Testing - Core.R")
