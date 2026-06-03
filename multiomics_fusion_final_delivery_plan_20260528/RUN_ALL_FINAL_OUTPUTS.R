# One-entry final output pipeline.
# Run this file in RStudio:
#   source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/RUN_ALL_FINAL_OUTPUTS.R")
#
# It does not modify original module3/module4 scripts. It only reads existing
# results and writes final organized outputs under:
#   D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/final_outputs

SCRIPT_ROOT <- "D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts"

scripts <- c(
  "00_prepare_delivery_dirs_and_method_flow.R",
  "09_fill_missing_master_table/build_bmmc_paired_scdrs_only_auroc.R",
  "09_fill_missing_master_table/build_nonpaired_trait_screening_scdrs_auroc.R",
  "01_build_master_result_tables.R",
  "02_plot_main_results.R",
  "03_build_ablation_parameter_outputs.R",
  "04_build_efficiency_outputs.R",
  "05_build_rna_reference_dependency_outputs.R",
  "06_prepare_case_study_assets.R",
  "07_compile_manuscript_text_templates.R",
  "08_build_final_outputs.R"
)

run_log <- data.frame(
  script = scripts,
  status = NA_character_,
  started_at = NA_character_,
  finished_at = NA_character_,
  error_message = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(scripts)) {
  script_path <- file.path(SCRIPT_ROOT, scripts[i])
  cat("\n============================================================\n")
  cat("Running: ", script_path, "\n", sep = "")
  cat("============================================================\n")

  run_log$started_at[i] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  res <- tryCatch({
    source(script_path, local = new.env(parent = globalenv()))
    TRUE
  }, error = function(e) {
    run_log$error_message[i] <<- conditionMessage(e)
    FALSE
  })
  run_log$finished_at[i] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  run_log$status[i] <- if (isTRUE(res)) "OK" else "FAILED"

  if (!isTRUE(res)) {
    cat("\nFAILED: ", scripts[i], "\n", sep = "")
    cat("Reason: ", run_log$error_message[i], "\n", sep = "")
    cat("Stopping pipeline. Fix the issue and rerun this file.\n")
    break
  }
}

log_dir <- "D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/final_outputs/99_source_manifest"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(run_log, file.path(log_dir, "final_output_pipeline_run_log.csv"), row.names = FALSE)

cat("\nPipeline log:\n")
print(run_log)

cat("\nHuman-readable final outputs are under:\n")
cat("  D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/final_outputs/最终输出\n")
