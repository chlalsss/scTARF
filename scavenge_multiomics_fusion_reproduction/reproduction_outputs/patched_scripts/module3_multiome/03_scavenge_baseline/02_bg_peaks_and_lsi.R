# 03_scavenge_baseline/02_bg_peaks_and_lsi.R
#
# Forked from blood_small.R Step 5 + Step 6. Identical logic, except the
# reference genome is BSgenome.Hsapiens.UCSC.hg38 (BMMC NeurIPS 2021 uses hg38;
# blood_small was hg19). Outputs are reused by 03_run_baseline_traits.R.
#
# Inputs:
#   <out_dir>/03_scavenge_baseline/SE_atac.rds   from 01_build_atac_se.R
#
# Outputs:
#   <out_dir>/03_scavenge_baseline/bg_peaks_n200.rds   chromVAR background peaks
#   <out_dir>/03_scavenge_baseline/atac_lsi.rds        LSI matrix (cells x dims)
#   <out_dir>/03_scavenge_baseline/mutualknn_k30.rds   mutual-KNN sparse adjacency
#
# Run time: addGCBias + getBackgroundPeaks ~3-8 min on 70k cells, LSI/KNN ~5-15 min.

source("D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome/_bootstrap.R")
required_pkgs(c(
  "SummarizedExperiment", "chromVAR", "BSgenome.Hsapiens.UCSC.hg38",
  "SCAVENGE", "Matrix"
))

paths <- m3_paths()
out_dir <- file.path(paths$out_dir, "03_scavenge_baseline")
SE_path <- file.path(out_dir, "SE_atac.rds")
if (!file.exists(SE_path)) stop("Missing ", SE_path, " — run 01_build_atac_se.R first.")

SE <- readRDS(SE_path)
message("Loaded SE: ", nrow(SE), " peaks x ", ncol(SE), " cells")

# --- Step 5: GC bias + background peaks (one-time) -----------------------

t0 <- Sys.time()
message("[Step 5] addGCBias (hg38) ...")
SE <- chromVAR::addGCBias(
  object = SE,
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

message("[Step 5] getBackgroundPeaks niter=200 ...")
bg_peaks <- chromVAR::getBackgroundPeaks(object = SE, niterations = 200)
saveRDS(bg_peaks, file.path(out_dir, "bg_peaks_n200.rds"))

# Save the SE with GC bias annotated, so the trait-loop script doesn't redo it
saveRDS(SE, file.path(out_dir, "SE_atac_gcbiased.rds"))
message("[Step 5] done in ", round(as.numeric(Sys.time() - t0, units = "secs")), " s")

# --- Step 6: LSI + mutual-KNN (one-time) ---------------------------------

t1 <- Sys.time()
message("[Step 6] tfidf -> do_lsi (dims=30) ...")
peak_by_cell <- SummarizedExperiment::assay(SE)
tfidf_mat <- SCAVENGE::tfidf(
  bmat       = peak_by_cell,
  mat_binary = TRUE,
  TF         = TRUE,
  log_TF     = TRUE
)
lsi_mat <- SCAVENGE::do_lsi(mat = tfidf_mat, dims = 30)
saveRDS(lsi_mat, file.path(out_dir, "atac_lsi.rds"))
message("LSI dim: ", nrow(lsi_mat), " x ", ncol(lsi_mat))

message("[Step 6] getmutualknn k=30 ...")
mknn <- SCAVENGE::getmutualknn(lsimat = lsi_mat, num_k = 30)
saveRDS(mknn, file.path(out_dir, "mutualknn_k30.rds"))
message("mutualKNN dim: ", nrow(mknn), " x ", ncol(mknn),
        " | nnz = ", length(mknn@x))
message("[Step 6] done in ", round(as.numeric(Sys.time() - t1, units = "secs")), " s")
