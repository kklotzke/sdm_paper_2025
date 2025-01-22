rm(list = ls())

N <- 100

tau.f3 <- 0.08
tau.sim <- c(tau.f3 , tau.f3, tau.f3, 0.3)
rep <- 250

set.seed(2)
filename = "../Results/Run_008_set1_N100.RDS"

source("SDM BF Testing - Core.R")
