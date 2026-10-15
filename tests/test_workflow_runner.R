source("scripts/run_workflow.R")
p <- workflow_plan("nonpaired")
position <- function(name) match(name, basename(p$script))
stopifnot(position("00_prepare_granja_gene_activity.R") < position("11_load_granja.R"),
          position("01_run_granja_atac_baseline.R") < position("54_granja_prepare_inputs.R"),
          position("54_granja_prepare_inputs.R") < position("55_eval_auroc.R"),
          position("02_build_granja_lsi.py") < position("03_run_v3_sparse_mapping_top5significant.R"),
          position("32_build_legacy_dense_mapping.R") < position("64_diagnose_scdrs_neighborhood_fdr.R"),
          position("64_diagnose_scdrs_neighborhood_fdr.R") < position("29_neighborhood_fdr_trait_screening.R"))
a <- workflow_plan("ablation")
stopifnot(match("02_build_master_result_tables.R", basename(a$script)) <
          match("01_build_ablation_parameter_outputs.R", basename(a$script)))
stopifnot(inherits(try(workflow_plan("typo"), silent = TRUE), "try-error"))
test_dir <- tempfile("runner space ")
dir.create(test_dir)
Sys.setenv(SCTARF_RUNNER_TEST_OUT = file.path(test_dir, "child result.txt"))
step <- file.path(test_dir, "step with spaces.R")
writeLines('writeLines("fresh process", Sys.getenv("SCTARF_RUNNER_TEST_OUT"))', step)
workflow_run_step(step, file.path(test_dir, "step.log"))
stopifnot(identical(readLines(Sys.getenv("SCTARF_RUNNER_TEST_OUT")), "fresh process"))
writeLines('stop("expected failure")', step)
stopifnot(inherits(try(workflow_run_step(step, file.path(test_dir, "failure.log")), silent = TRUE), "try-error"))
cat("Workflow order, subprocess paths, and failure propagation: PASS\n")
