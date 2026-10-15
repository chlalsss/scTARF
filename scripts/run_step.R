# The workflow runner invokes this wrapper in a fresh R process.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Internal usage: run_step.R SCRIPT SESSION_LOG")
step_path <- args[1]
step_session_log <- args[2]
if (file.exists("paths.local.R")) source("paths.local.R")
tryCatch(source(step_path, local = .GlobalEnv, encoding = "UTF-8"),
         finally = capture.output(sessionInfo(), file = step_session_log))
