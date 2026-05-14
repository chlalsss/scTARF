library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BiocManager)
library(igraph)

set.seed(9527)

# 初始化总计时
script_start <- Sys.time()

# =========================
# 0) 路径与输出目录
# =========================
data_dir  <- "D:/SCAVENGE-reproducibility-main/blood_data"
se_path   <- file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
gwas_path <- file.path(data_dir, "CS.trans_all.PP01.strong_atac_overlap.tsv")

# 生成当前时间字符串，例如 "2026-03-06_13-45"
time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")

# 将时间戳加入目录名
out_dir <- file.path(data_dir, paste0("scavenge_out_", time_stamp))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1) 读入小数据集（改名：造血 scATAC 小数据集）
# =========================
# ### [不同] getstarted 用 example_data(pbmc5k)，这里改为读取本地 rds
SE_heme_small <- readRDS(se_path)
cat("[SE] cells =", ncol(SE_heme_small), " peaks =", nrow(SE_heme_small), "\n")

# =========================
# 2) 读入 GWAS（改名：BCX2 fine-mapping 汇总）
# =========================
# ### [不同] getstarted 用 example_data(mono.PP001.bed)，这里改为读取你们 tsv
gwas_bcx2 <- data.table::fread(gwas_path)
cat("[GWAS] rows =", nrow(gwas_bcx2), " cols =", ncol(gwas_bcx2), "\n")
print(names(gwas_bcx2))

# =========================
# 3) 固定列名（根据你截图）
# =========================
# ### [新增] getstarted 不需要列名映射（示例 bed 已经准备好）；你们需要指定 tsv 的关键列
col_var   <- "var"
col_trait <- "trait"
col_pp    <- "PP"

stopifnot(all(c(col_var, col_trait, col_pp) %in% names(gwas_bcx2)))

# =========================
# 4) 从 var 解析 chr/start/end
# =========================
# 只解析 var 形如 "1:2164935_G_T"（chr:pos_ref_alt）
# 输出 BED 坐标：chr=chr1, start=pos-1, end=pos
parse_var_to_bed <- function(v) {
  v <- as.character(v)
  
  m <- regexec("^([0-9]+|X|Y|MT|M):(\\d+)_([^_]+)_([^_]+)$", v, ignore.case = TRUE)
  r <- regmatches(v, m)
  
  ok <- lengths(r) == 5
  chr <- start <- end <- rep(NA_character_, length(v))
  
  if (any(ok)) {
    raw_chr <- toupper(vapply(r[ok], `[[`, "", 2))
    raw_chr[raw_chr == "M"] <- "MT"       # 兼容 M 写法
    chr[ok] <- paste0("chr", raw_chr)
    
    pos <- as.integer(vapply(r[ok], `[[`, "", 3))
    start[ok] <- as.character(pmax(pos - 1L, 0L))
    end[ok]   <- as.character(pos)
  }
  
  data.table(
    chr = chr,
    start = as.integer(start),
    end = as.integer(end),
    ok = !is.na(chr) & !is.na(start) & !is.na(end)
  )
}

# =========================
# 4.5) 把 var 解析结果写回 gwas_bcx2（必须！）
# =========================
# ### [新增] getstarted 不需要这一步（因为直接给 bed）；你们必须把 var -> chr/start/end
bed_parsed <- parse_var_to_bed(gwas_bcx2[[col_var]])
bad_n <- sum(!bed_parsed$ok)
cat("[var parse] bad =", bad_n, "\n")
stopifnot(bad_n == 0)

gwas_bcx2[, `:=`(
  chr   = bed_parsed$chr,
  start = bed_parsed$start,
  end   = bed_parsed$end
)]
gwas_bcx2[end <= start, end := start + 1L]

# 保险检查
stopifnot(all(c("chr","start","end") %in% names(gwas_bcx2)))

# =========================
# 5) GC bias + background peaks（只做一次）
# =========================
# ### [不同] getstarted 的对象名是 SE_pbmc5k；这里统一用 SE_heme_small

step5_start <- Sys.time()
cat("[Step 5] 开始计算 GC bias 和 Background Peaks...\n")

cat("[Step] addGCBias...\n")
SE_heme_small <- addGCBias(
  object = SE_heme_small,
  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
)

cat("[Step] getBackgroundPeaks...\n")
heme_bg_peaks <- getBackgroundPeaks(object = SE_heme_small, niterations = 200)
saveRDS(heme_bg_peaks, file.path(out_dir, "bg_peaks_n200.rds"))

step5_end <- Sys.time()
cat(sprintf("[计时] Step 5 耗时: %.2f 秒\n", as.numeric(difftime(step5_end, step5_start, units = "secs"))))


# =========================
# 6) 构建 m-KNN 图（只做一次）
# =========================
# ### [不同] getstarted 单性状示例里也是这样做；这里仍然相同，但对象名不同

step6_start <- Sys.time()
cat("[Step 6] 开始构建 m-KNN 图 (LSI + MutualKNN)...\n")

cat("[Step] build mKNN graph (tfidf -> lsi -> mutualknn)...\n")
peak_by_cell_mat <- SummarizedExperiment::assay(SE_heme_small)

tfidf_mat <- tfidf(
  bmat = peak_by_cell_mat,
  mat_binary = TRUE,
  TF = TRUE,
  log_TF = TRUE
)

lsi_mat <- do_lsi(mat = tfidf_mat, dims = 30)

heme_mknn_k30 <- getmutualknn(lsimat = lsi_mat, num_k = 30)
saveRDS(heme_mknn_k30, file.path(out_dir, "mutualknn_k30.rds"))

step6_end <- Sys.time()
cat(sprintf("[计时] Step 6 耗时: %.2f 秒\n", as.numeric(difftime(step6_end, step6_start, units = "secs"))))


# =========================
# 7) 并行设置
# =========================
# ### [新增] getstarted 没写并行；这里为了跑 22 traits 更快
mycores <- 8
bp <- BiocParallel::MulticoreParam(workers = mycores)

# =========================
# 8) 全性状循环（核心改动）
# =========================
# ### [不同] getstarted 只跑 mono 一个性状；这里跑所有 trait

loop_start <- Sys.time()

traits <- sort(unique(gwas_bcx2[[col_trait]]))
cat("[Traits] n =", length(traits), "\n")

trs_list <- vector("list", length(traits))
names(trs_list) <- traits

cell_ids <- colnames(SE_heme_small)

for (tr in traits) {
  cat("\n============================\n")
  cat("[Trait]", tr, "\n")
  
  dt <- gwas_bcx2[get(col_trait) == tr]
  
  # 生成 bed-like：chr start end name score
  bed <- dt[, .(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end),
    name = tr,
    score = as.numeric(get(col_pp))
  )]
  bed[end <= start, end := start + 1L]
  
  # 写 bed 文件（每 trait 一个）
  # ### [新增] getstarted 的 bed 是 example_data 直接给；这里要自己写
  bed_file <- file.path(out_dir, paste0("trait_", make.names(tr), ".PP001.bed"))
  fwrite(bed, bed_file, sep = "\t", col.names = FALSE)
  
  # import weights
  trait_import <- importBedScore(rowRanges(SE_heme_small), bed_file, colidx = 5)
  
  # gchromVAR deviations
  SE_DEV <- computeWeightedDeviations(
    object = SE_heme_small,
    weights = trait_import,
    background_peaks = heme_bg_peaks
  )
  
  z_score <- t(SummarizedExperiment::assays(SE_DEV)[["z"]]) |> c()
  
  # seed + scale_factor（与 getstarted 相同）
  seed_idx <- seedindex(z_score = z_score, percent_cut = 0.05)
  scale_factor <- cal_scalefactor(z_score = z_score, percent_cut = 0.01)
  
  # random walk（与 getstarted 相同）
  np_score <- randomWalk_sparse(
    intM = heme_mknn_k30,
    queryCells = rownames(heme_mknn_k30)[seed_idx],
    gamma = 0.05
  )
  
  omit_idx <- np_score == 0
  cat("[Trait]", tr, "omit =", sum(omit_idx), "\n")
  
  np_score_use <- np_score[!omit_idx]
  
  TRS <- capOutlierQuantile(x = np_score_use, q_ceiling = 0.95) |> max_min_scale()
  TRS <- TRS * scale_factor
  
  keep_cells <- names(np_score)[!omit_idx]
  
  # 输出表
  out_df <- data.frame(
    cell = keep_cells,
    z_score = z_score[match(keep_cells, cell_ids)],
    seed = seed_idx[match(keep_cells, rownames(heme_mknn_k30))],
    np_score = np_score_use,
    TRS = TRS,
    stringsAsFactors = FALSE
  )
  
  # 拼接 meta（可选）
  meta <- as.data.frame(colData(SE_heme_small))
  meta$cell <- rownames(meta)
  out_df <- dplyr::left_join(out_df, meta, by = "cell")
  
  # 写出每 trait 结果
  out_file <- file.path(out_dir, paste0("SCAVENGE_small_", make.names(tr), ".tsv.gz"))
  fwrite(out_df, out_file, sep = "\t")
  cat("[Trait]", tr, "saved:", out_file, "\n")
  
  # 汇总 TRS
  trs_list[[tr]] <- out_df[, c("cell", "TRS")]
  colnames(trs_list[[tr]]) <- c("cell", make.names(tr))
}

loop_end <- Sys.time()
cat(sprintf("\n[计时] 所有性状 (n=%d) 总耗时: %.2f 秒\n", length(traits), as.numeric(difftime(loop_end, loop_start, units = "secs"))))

# =========================
# 9) 合并 TRS 矩阵并输出（可选但建议）
# =========================
cat("\n[Step] merge TRS matrix...\n")
trs_mat <- Reduce(function(x, y) dplyr::full_join(x, y, by = "cell"), trs_list)

saveRDS(trs_mat, file.path(out_dir, "TRS_matrix_small_all_traits.rds"))
fwrite(trs_mat, file.path(out_dir, "TRS_matrix_small_all_traits.tsv.gz"), sep = "\t")

cat("[DONE] Outputs in:", out_dir, "\n")




# 结果文件数（应该 >= 15）
length(list.files(out_dir, pattern = "^SCAVENGE_small_.*\\.tsv\\.gz$"))

# TRS 矩阵维度（行数≈细胞数，列数=1+trait数）
dim(trs_mat)

# 看一个 trait 的 TRS 是否不是全 0
summary(trs_mat[[2]])




# 全部细胞 ID（来自 SE）
all_cells <- colnames(SE_heme_small)

# 以 all_cells 为主键补齐
trs_full <- data.frame(cell = all_cells, stringsAsFactors = FALSE) |>
  dplyr::left_join(trs_mat, by = "cell")

# 把所有 trait 列的 NA 变成 0
trait_cols <- setdiff(colnames(trs_full), "cell")
for (cc in trait_cols) trs_full[[cc]][is.na(trs_full[[cc]])] <- 0

# 保存“全细胞版”
saveRDS(trs_full, file.path(out_dir, "TRS_matrix_small_all_traits_fullcells.rds"))
data.table::fwrite(trs_full, file.path(out_dir, "TRS_matrix_small_all_traits_fullcells.tsv.gz"), sep="\t")

dim(trs_full)  # 期望 35038 × 16

script_end <- Sys.time()
cat("\n==============================================\n")
cat(sprintf("整个脚本运行结束，累计总耗时: %.2f 秒\n", as.numeric(difftime(script_end, script_start, units = "secs"))))
cat("==============================================\n")
