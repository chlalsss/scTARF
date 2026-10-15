# Rscript --vanilla paper_reproduction/run_paper.R <new output directory> [repo root]
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(dirname(script), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) stop("Provide a new output directory, optionally followed by repository root")
out <- args[1]
if (dir.exists(out) || file.exists(out)) stop("Output already exists; choose a new directory: ", out)
repo <- if (length(args) == 2L) args[2] else file.path(root, "..")
check_files <- function(manifest, base, source_text = FALSE) {
  paths <- file.path(base, manifest$path)
  if (source_text) {
    actual <- vapply(paths, function(path) {
      bytes <- readBin(path, "raw", n = file.info(path)$size)
      normalized <- charToRaw(gsub("\r\n", "\n", rawToChar(bytes), fixed = TRUE, useBytes = TRUE))
      temporary <- tempfile("source-check-")
      writeBin(normalized, temporary)
      value <- unname(tools::md5sum(temporary))
      unlink(temporary)
      value
    }, character(1), USE.NAMES = FALSE)
  } else actual <- unname(tools::md5sum(paths))
  stopifnot(!anyNA(actual), identical(actual, manifest$md5))
}
check_files(read.delim(file.path(root, "input_checksums.tsv"), stringsAsFactors = FALSE), root)
check_files(read.delim(file.path(root, "required_sources.tsv"), stringsAsFactors = FALSE), repo, source_text = TRUE)
source(file.path(root, "replay_trs.R"), local = TRUE)
source(file.path(root, "evaluate_results.R"), local = TRUE)
cases <- read.delim(file.path(root, "replay_manifest.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(cases) == 30L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
summaries <- list()
for (i in seq_len(nrow(cases))) {
  item <- cases[i, ]
  cat("REPLAY ", item$id, "\n", sep = "")
  result <- replay_trs(repo, file.path(root, item$graph_file), file.path(root, item$seed_file), file.path(root, item$reference_file))
  group <- if (startsWith(item$setting, "nonpaired")) "nonpaired" else item$setting
  directory <- file.path(out, group)
  dir.create(directory, showWarnings = FALSE)
  write.table(result$scores, file.path(directory, paste0(item$id, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  summaries[[i]] <- data.frame(id = item$id, n_cells = result$n_cells,
    max_abs_difference = result$max_abs_difference, matched = result$matched)
}
summary <- do.call(rbind, summaries)
write.csv(summary, file.path(out, "replay_comparison.csv"), row.names = FALSE)
stopifnot(all(summary$matched))
evaluate_paper(file.path(root, "data"), out, file.path(out, "auroc_comparison.csv"))
writeLines(c("Manuscript replay from frozen graphs/seeds and historical baseline scores.", capture.output(sessionInfo())), file.path(out, "session_info.txt"))
cat("PASS: 30 cell-score replays and 61 historical-source AUROC checks\n")
