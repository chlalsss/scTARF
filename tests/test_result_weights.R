scratch <- tempfile("weight-table-test ")
dir.create(scratch)
Sys.setenv(M3_OUT_DIR = file.path(scratch, "paired"),
           SCTARF_NONPAIRED_ROOT = file.path(scratch, "nonpaired"),
           SCTARF_DELIVERY_ROOT = file.path(scratch, "delivery"))
input <- file.path(scratch, "paired", "04_dual_modality", "eval_v2_auroc.csv")
dir.create(dirname(input), recursive = TRUE)
weights <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
fixture <- data.frame(method = c("baseline_atac", rep("v2", 6)), target = "MONO",
                      atac_w = c(NA_real_, weights), auroc = seq(0.6, 0.9, length.out = 7))
write.csv(fixture, input, row.names = FALSE)
source("scripts/07_figures_tables/02_build_master_result_tables.R")
actual <- read.csv(file.path(scratch, "delivery/final_outputs/01_master_tables/main_result_comparison_table.csv"))
stopifnot(identical(sort(actual$atac_w[actual$method == "v2"]), weights),
          is.na(actual$atac_w[actual$method == "baseline_atac"]))
cat("Paired weight sweep preserved in master results: PASS\n")
