rm(list = ls())

# Load data
load(file = "../RISE.Rdata")

# MD2 factor structure
u.fs <- u[7:8, ]

# HMC tuning parameters 
epsilon <- list() 
epsilon[[1]] <- 0.3
epsilon[[2]] <- 0.3
epsilon[[3]] <- 0.3
epsilon[[4]] <- 0.3

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

ftype = "MD"
filename = "Results/Run_MD2_VS.RDS"

# Fit model 
set.seed(1)
source("RISE Part2 - Core.R")
