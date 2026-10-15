# Run from the repository root. Each manifest step gets a fresh process.
workflow_plan <- function(stage) {
  manifest <- read.delim("scripts/final_reproduction_manifest.tsv", stringsAsFactors = FALSE)
  if (!stage %in% manifest$stage) {
    stop("Unknown stage. Choose: ", paste(unique(manifest$stage), collapse = ", "), call. = FALSE)
  }
  plan <- manifest[manifest$stage == stage & manifest$role == "run", , drop = FALSE]
  plan <- plan[order(plan$order), , drop = FALSE]
  if (anyDuplicated(plan$order)) stop("Duplicate manifest order in stage: ", stage)
  if (!all(file.exists(plan$script))) stop("Manifest references missing scripts.")
  plan
}

workflow_run_step <- function(script, logfile) {
  source("scripts/_paths.R")
  if (grepl("[.]R$", script)) {
    executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    args <- c("--vanilla", "scripts/run_step.R", script, paste0(logfile, ".session.txt"))
  } else if (grepl("[.]py$", script)) {
    executable <- sctarf_python()
    args <- script
  } else stop("Unsupported script type: ", script)
  status <- system2(executable, shQuote(args), stdout = logfile, stderr = logfile)
  if (status != 0L) stop("Step failed (", status, "): ", script, "\nSee: ", logfile, call. = FALSE)
  invisible(status)
}

workflow_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) < 1L || length(args) > 2L ||
      (length(args) == 2L && args[2] != "--dry-run")) {
    stop("Usage: Rscript --vanilla scripts/run_workflow.R STAGE [--dry-run]", call. = FALSE)
  }
  plan <- workflow_plan(args[1])
  print(plan[, c("order", "script", "purpose")], row.names = FALSE)
  if ("--dry-run" %in% args) return(invisible(plan))
  if (file.exists("paths.local.R")) source("paths.local.R")
  source("scripts/_paths.R")
  logdir <- file.path(sctarf_output_root(), "run_logs", paste0(args[1], "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
  write.table(plan, file.path(logdir, "execution_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines(c(paste("Started:", Sys.time()), capture.output(Sys.info()), capture.output(sessionInfo())),
             file.path(logdir, "environment.txt"))
  for (i in seq_len(nrow(plan))) {
    script <- plan$script[i]
    message("[", i, "/", nrow(plan), "] ", script)
    workflow_run_step(script, file.path(logdir, sprintf("%02d_%s.log", i, basename(script))))
  }
  message("Stage completed. Logs: ", logdir)
  invisible(plan)
}

if (sys.nframe() == 0L) workflow_main()
