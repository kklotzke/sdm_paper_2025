# Konrad Klotzke

# Randomly select number of items from (non-overlapping) dimensions under multidimensional or bi-factor model
# u: Classification matrix
# nItems: Number of items to select in each cluster
selectRandomItemsBF<- function(u, nItems) {
  Nt <- nrow(u)
  p <- ncol(u)
  u.new <- matrix(0, nrow = Nt, ncol = (Nt - 1) * nItems)
  p.new <- ncol(u.new)
  selected.items <- vector("list", length = Nt)
  which.f1 <- NA
  items <- numeric()
  items.new <- 1:p.new 
  for (tt in 1:Nt) {
    vec <- which(u[tt, ] == 1)
    if (length(vec) == p) {
      which.f1 <- tt
      u.new[which.f1, ] <- 1
    }
    else {
      selected.items[[tt]] <- vec <- sample(vec, nItems)
      items <- c(items, vec)

    }
  }
    
  items <- sort(items)
  print (items)
  print(which.f1)
    
  for (tt in 1:Nt) {
    vec <- selected.items[[tt]]
    if (length(vec) != p) {
      for (kk in 1:p.new) {
        if (items[kk] %in% vec) {
          u.new[tt, kk] <- 1 
        }
      }
    }
  }
  
  rownames(u.new) <- rownames(u)

  return(list("u" = u.new, "selected.items.cluster" = selected.items, "selected.items.all" = items, "nItems" = nItems))
}


# ## Test
# u <- matrix(0, 3, 7)
# u[1, ] <- c(1, 0, 0, 1, 0, 0, 1)
# u[2, ] <- c(0, 0, 1, 0, 1, 1, 0)
# u[3, ] <- 1
# 
# sri <- selectRandomItemsBF(u, 2)
# sri$selected.items.cluster
# sri$u

FM.threshold <- function(M0, M1, thr, smaller = TRUE) {
  rep <- length(M0)
  select <- matrix(0, nrow = rep, ncol = 2) # Select M0, Select M1
  for (rr in 1:rep) {
    if (all(c(M0[rr], M1[rr]) < thr)) { # If both estimates are below threshold, select simpler model
      select[rr, 1] <- 1 # Select M0
    }
    else if ((M0[rr] >= thr) && (M1[rr] < thr)) { 
      select[rr, 2] <- 1 # Select M1
    }
  }
  if(rep > 1) {
    select.noties <- select[rowSums(select[, ]) > 0, ]
    if(is.null(nrow(select.noties))) {
      select.noties <- matrix(select.noties, nrow = 1, ncol = 2)
    }
  } 
  else {
    if (all(select == 0)) {
      select.noties <- matrix(c(NA, NA), nrow = 1, ncol = 2)
    }
    else {
      select.noties <- select
    }
  }
  return(list("select" = select, "select.noties" = select.noties))
} 

FM.delta <- function(M0, M1) {
  rep <- length(M0)
  select <- matrix(0, nrow = rep, ncol = 2) # Select M0, Select M1
  for (rr in 1:rep) {
    M.delta <- M1[rr] - M0[rr]
    if (M.delta >= 0) { # If tied, select simpler model
      select[rr, 1] <- 1 # Select M0
    }
    else {
      select[rr, 2] <- 1 # Select M1
    }
  }
  return(list("select" = select))
} 

colSDs <- function(x) sqrt(rowMeans((t(x) - colMeans(x)) ^ 2) * ((dim(x)[1]) / (dim(x)[1] - 1)))
