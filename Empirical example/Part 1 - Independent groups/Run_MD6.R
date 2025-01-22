rm(list = ls())

# Load data
load(file = "../RISE.Rdata")

# MD6 factor structure
u.fs <- u[1:6, ]

# HMC tuning parameters 
epsilon <- 0.2
epsilon.scale <- 1
epsilon.sd <- 0#epsilon * 0.02
L <- 6
L.sd <- 1 

ftype = "MD"
filename = "Results/Run_MD6.RDS"

state <- readRDS(filename)
which.reps <- 1:4

# Fit model 
set.seed(1)
source("RISE Part1 - Core.R")
