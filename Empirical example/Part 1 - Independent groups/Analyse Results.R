source("../Tools.R")

results.UD <- readRDS("Results/Run_UD.RDS")
results.MD6 <- readRDS("Results/Run_MD6.RDS")
results.MD2 <- readRDS("Results/Run_MD2.RDS")
results.BF6 <- readRDS("Results/Run_BF6.RDS")
results.BF2 <- readRDS("Results/Run_BF2.RDS")

table_results <- matrix(NA, nrow = 5, ncol = 4)
rownames(table_results) <- c("UD", "MD6", "MD2", "BF6", "BF2")
colnames(table_results) <- c("Grade 6", "Grade 7", "Grade 8", "Grade 9")
table_results.ESS <- table_results

table_results[1, ] <- results.UD$marllik
table_results[2, ] <- results.MD6$marllik
table_results[3, ] <- results.MD2$marllik
table_results[4, ] <- results.BF6$marllik
table_results[5, ] <- results.BF2$marllik

table_results.ESS[1, ] <- results.UD$marllik.ESS
table_results.ESS[2, ] <- results.MD6$marllik.ESS
table_results.ESS[3, ] <- results.MD2$marllik.ESS
table_results.ESS[4, ] <- results.BF6$marllik.ESS
table_results.ESS[5, ] <- results.BF2$marllik.ESS

table_names <- matrix(rownames(table_results), nrow = 5, ncol = 4)
table_sorted <- matrix(table_names[apply(table_results, 2, order, decreasing = TRUE)], nrow = 5, ncol = 4)
table_sorted.ESS <- matrix(table_names[apply(table_results.ESS, 2, order, decreasing = TRUE)], nrow = 5, ncol = 4)
colnames(table_sorted) <- colnames(table_sorted.ESS) <- c("Grade 6", "Grade 7", "Grade 8", "Grade 9")

table_results
talbe_results.ESS

table_sorted # Sorted by log marginal likelihood per group (decreasing order)
table_sorted.ESS