rm(list = ls())

# Load data
load(file = "../RISE.Rdata")

# Unidimensional factor structure
u.fs <- matrix(1, nrow = 1, ncol = 204)

# HMC tuning parameters 
epsilon <- 0.4
epsilon.scale <- 1
epsilon.sd <- 0#epsilon * 0.02
L <- 6
L.sd <- 1 

ftype = "UD"
filename = "Results/Run_UD.RDS"

state <- readRDS(filename)
which.reps <- 1:4

# Fit model 
set.seed(1)
source("RISE Part1 - Core.R")
