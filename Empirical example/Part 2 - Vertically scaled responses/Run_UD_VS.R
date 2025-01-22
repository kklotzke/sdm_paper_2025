rm(list = ls())

# Load data
load(file = "../RISE.Rdata")

# Unidimensional factor structure
u.fs <- matrix(1, nrow = 1, ncol = 204)

# HMC tuning parameters 
epsilon <- list() 
epsilon[[1]] <- 0.2
epsilon[[2]] <- 0.2
epsilon[[3]] <- 0.2
epsilon[[4]] <- 0.2

epsilon.scale <- list() 
epsilon.scale[[1]] <- 1
epsilon.scale[[2]] <- 1
epsilon.scale[[3]] <- 1
epsilon.scale[[4]] <- 1

epsilon.sd <- list() 
epsilon.sd[[1]] <- 0
epsilon.sd[[2]] <- 0
epsilon.sd[[3]] <- 0
epsilon.sd[[4]] <- 0

L <- list() 
L[[1]] <- 6
L[[2]] <- 6
L[[3]] <- 6
L[[4]] <- 6

L.sd <- list() 
L.sd[[1]] <- 1
L.sd[[2]] <- 1
L.sd[[3]] <- 1
L.sd[[4]] <- 1

ftype = "UD"
filename = "Results/Run_UD_VS.RDS"

# Fit model 
set.seed(1)
source("RISE Part2 - Core.R")
