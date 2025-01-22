rm(list = ls())

# Load data
load(file = "../RISE.Rdata")

# BF6 factor structure
u.fs <- u[c(1:6,9), ]

fix.ginter <- TRUE
fix.tinter <- FALSE
fix.cor <- FALSE

K <- 60
Nt <- nrow(u.fs)
Ng <- 4 

# Number of parameters
Np <- Ng*Nt + Ng*K + 0 + (Ng-1)*(Nt-1) # Ng*Nt covariance parameter(s), Ng*K item intercepts, 0 group intercepts, (Ng-1)*(Nt-1) subtest intercepts 

# HMC tuning parameters 
epsilon <- list() 
epsilon[[1]] <- 0.065
epsilon[[2]] <- 0.065
epsilon[[3]] <- 0.065
epsilon[[4]] <- 0.065

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

ftype = "BF"
filename = "Results/Run_BF6_VS_Fix_ginter.RDS"

# Fit model 
set.seed(1)
source("RISE Part3 - Core.R")
