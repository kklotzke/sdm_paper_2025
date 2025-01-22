# Draw random sample from RISE data
# 500 students per group
# 10 items per subtest

source("../Tools.R")

# Load data
# load(file = "../RISE.Rdata")

# Group 6-9
Y6 <- data6[,-c(1,2)]
Y7 <- data6[,-c(1,2)]
Y8 <- data6[,-c(1,2)]
Y9 <- data6[,-c(1,2)]

# Randomly select 10 items per subtest
set.seed(12)
nItems <- 10
u.MD6 <- u[1:6, ] # MD6 structure specifies the six subtests 
sri <- selectRandomItemsBF(u.MD6, nItems)
set <- sri$selected.items.all

Y6 <- Y6[, set] 
Y7 <- Y7[, set] 
Y8 <- Y8[, set] 
Y9 <- Y9[, set] 
u <- u.fs[, set] 
if (nrow(u.fs) == 1) {
  u <- matrix(u, nrow = 1, ncol = length(set))
}

# Randomly select 500 students per grade
nStudents <- 500
Y6 <- Y6[sample(1:nrow(Y6), nStudents, replace = FALSE), ]
Y7 <- Y7[sample(1:nrow(Y7), nStudents, replace = FALSE), ]
Y8 <- Y8[sample(1:nrow(Y8), nStudents, replace = FALSE), ]
Y9 <- Y9[sample(1:nrow(Y9), nStudents, replace = FALSE), ]