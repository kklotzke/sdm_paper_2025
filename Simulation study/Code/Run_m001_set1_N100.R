rm(list = ls())

N <- 100

tau.f3 <- -0.01
tau.sim <- c(tau.f3 , tau.f3, tau.f3, 0.3)
rep <- 250

set.seed(2)
filename = "../Results/Run_m001_set1_N100.RDS"

#state <- readRDS(filename)
#which.reps <- 35:rep#(state$summary.lavaan[1]+1):rep # Resume calculations
source("SDM BF Testing - Core.R")
