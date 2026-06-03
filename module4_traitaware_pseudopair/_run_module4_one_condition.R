trait <- Sys.getenv("MODULE4_TRAIT", unset = "MCH")
p_drop <- as.numeric(Sys.getenv("MODULE4_P_DROP", unset = "0"))
seed <- as.integer(Sys.getenv("MODULE4_SEED", unset = "1"))
batch_size <- as.integer(Sys.getenv("MODULE4_BATCH_SIZE", unset = "64"))

setwd("D:/RStudio/work/module4_traitaware_pseudopair")

cat("=== module4 one-condition runner ===\n")
cat("trait      : ", trait, "\n", sep = "")
cat("p_drop     : ", p_drop, "\n", sep = "")
cat("seed       : ", seed, "\n", sep = "")
cat("batch_size : ", batch_size, "\n", sep = "")
cat("started_at : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")

MODULE4_FIXED_TRAITS_TO_RUN <- trait
MODULE4_ROBUST_P_DROP_LIST <- p_drop
MODULE4_ROBUST_SEED_LIST <- seed
MODULE4_ROBUST_SPARSE_MAP_BATCH_SIZE <- batch_size

source("RUN_MODULE4_FIXED_FLOW_ROBUSTNESS.R")

cat("\nfinished_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
