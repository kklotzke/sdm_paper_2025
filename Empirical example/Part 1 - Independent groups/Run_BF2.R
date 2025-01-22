rm(list = ls())

# Load data
load(file = "../RISE.Rdata")

# BF2 factor structure
u.fs <- u[7:9, ]

# HMC tuning parameters 
epsilon <- 0.1
epsilon.scale <- 1
epsilon.sd <- 0#epsilon * 0.02
L <- 6
L.sd <- 1 

ftype = "BF"
filename = "Results/Run_BF2.RDS"

state <- readRDS(filename)
which.reps <- 1:4

# Fit model 
set.seed(1)
source("RISE Part1 - Core.R")
