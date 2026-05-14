# multiomics_phase3to6.R
# 阶段 3-6：跨模态整合 → Label Transfer → TRS 调控解释 → 交叉验证
#
# 前置条件：阶段 2 的产出文件已存在于 multiomics_output/
# 运行方式：source("D:/RStudio/work/multiomics_phase3to6.R")

library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)
library(harmony)

set.seed(9527)
script_start <- Sys.time()

# =========================
# 0) 路径与加载
# =========================
data_dir <- "D:/SCAVENGE-reproducibility-main/blood_data"
out_dir  <- file.path(data_dir, "multiomics_output_v2")

cat("========== 清理内存 + 加载数据 ==========\n")
rm(list = ls()[!ls() %in% c("data_dir", "out_dir", "script_start")]); gc()

# 第一步：获取 ATAC metadata（首次从 SE 提取并缓存，之后直接读小文件）
atac_meta_path <- file.path(out_dir, "atac_metadata.rds")
if (!file.exists(atac_meta_path)) {
  cat("[1] 首次提取 ATAC metadata...\n")
  SE_atac <- readRDS(file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds"))
  atac_meta <- as.data.frame(colData(SE_atac))
  saveRDS(atac_meta, atac_meta_path, compress = FALSE)
  rm(SE_atac); gc()
} else {
  cat("[1] 加载已缓存的 ATAC metadata...\n")
  atac_meta <- readRDS(atac_meta_path)
}
cat("  ATAC metadata:", nrow(atac_meta), "cells\n")

# 第二步：加载 Gene Activity + scRNA
cat("[2] 加载 Gene Activity 矩阵...\n")
gene_activities <- readRDS(file.path(out_dir, "gene_activity_matrix.rds"))
cat("  Gene Activity:", nrow(gene_activities), "genes x", ncol(gene_activities), "cells\n")

cat("[3] 加载 scRNA Seurat...\n")
seurat_rna <- readRDS(file.path(out_dir, "seurat_rna_normalized.rds"))
cat("  scRNA Seurat:", ncol(seurat_rna), "cells x", nrow(seurat_rna), "genes\n")

# 第三步：构建轻量 ATAC Seurat
seurat_atac <- CreateSeuratObject(
  counts    = gene_activities,
  assay     = "RNA_from_ATAC",
  meta.data = atac_meta
)
seurat_atac <- NormalizeData(seurat_atac)
cat("  ATAC Seurat (轻量GA版):", ncol(seurat_atac), "cells\n")


# ╔═══════════════════════════════════════════════════════════╗
# ║  阶段 3：跨模态整合（pcaproject + Anchor）               ║
# ╚═══════════════════════════════════════════════════════════╝
cat("\n\n")
cat("############################################\n")
cat("#  阶段 3：跨模态整合 (pcaproject + Anchor)#\n")
cat("############################################\n")

step3_start <- Sys.time()

# 3a) 以 scRNA 的高变基因为 anchor 特征，只保留两边都有的基因
hvg <- VariableFeatures(seurat_rna)
shared_genes <- intersect(rownames(seurat_atac), rownames(seurat_rna))
hvg_use <- intersect(hvg, shared_genes)
cat("用于 anchor 的高变基因:", length(hvg_use), "\n")

# 3b) FindTransferAnchors
#     用 pcaproject（比 CCA 快 3-5 倍，同体系数据精度接近）
cat("[Step 3b] FindTransferAnchors (pcaproject)...\n")
transfer_anchors <- FindTransferAnchors(
  reference       = seurat_rna,
  query           = seurat_atac,
  features        = hvg_use,
  reference.assay = "RNA",
  query.assay     = "RNA_from_ATAC",
  reduction       = "pcaproject",
  dims            = 1:30
)
cat("Anchor 数量:", nrow(transfer_anchors@anchors), "\n")

# 3c) TransferData: 将 RNA 的细胞类型标签 transfer 到 ATAC
#     关键：结果存到新变量 transfer_result，不要覆盖 seurat_atac
cat("[Step 3c] TransferData (Label Transfer)...\n")
transfer_result <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = seurat_rna$BioClassification,
  weight.reduction = "pcaproject",
  dims             = 1:30
)

# 从返回的 Seurat 对象中提取预测标签和置信度
seurat_atac$predicted.id         <- transfer_result$predicted.id
seurat_atac$prediction.score.max <- transfer_result$prediction.score.max
cat("Label Transfer 完成\n")

# 释放 transfer_result 节省内存
rm(transfer_result); gc()

# 3d) Co-embedding: 采样后合并做联合降维（纯可视化用途，不影响 DE 分析）
cat("[Step 3d] Co-embedding (subsample + merge + PCA + UMAP)...\n")

# 释放大对象腾内存
rm(gene_activities, transfer_anchors); gc()

# 给两边加上标记
seurat_rna$predicted.id <- seurat_rna$BioClassification
seurat_rna$modality <- "scRNA"

atac_merge_meta <- seurat_atac@meta.data
atac_merge_meta$modality <- "scATAC"

ga_counts <- GetAssayData(seurat_atac, assay = "RNA_from_ATAC", layer = "counts")

# 每模态采样 10000 个细胞（70k → 20k，内存降至约 1/3）
n_sub <- 10000
set.seed(9527)
rna_cells_sub  <- sample(colnames(seurat_rna), min(n_sub, ncol(seurat_rna)))
atac_cells_sub <- sample(colnames(seurat_atac), min(n_sub, ncol(seurat_atac)))
cat("  Co-embedding 采样: scRNA", length(rna_cells_sub),
    " scATAC", length(atac_cells_sub), "\n")

seurat_atac_for_merge <- CreateSeuratObject(
  counts    = ga_counts[, atac_cells_sub],
  assay     = "RNA",
  meta.data = atac_merge_meta[atac_cells_sub, ]
)
seurat_atac_for_merge <- NormalizeData(seurat_atac_for_merge)

seurat_rna_sub <- subset(seurat_rna, cells = rna_cells_sub)

# 合并 RNA + ATAC
co_embed <- merge(seurat_rna_sub, seurat_atac_for_merge)
rm(seurat_atac_for_merge, seurat_rna_sub, ga_counts); gc()

DefaultAssay(co_embed) <- "RNA"
co_embed <- JoinLayers(co_embed)
co_embed <- FindVariableFeatures(co_embed, nfeatures = 3000)
co_embed <- ScaleData(co_embed)
co_embed <- RunPCA(co_embed, npcs = 30)

# Harmony 批次校正：消除 scATAC / scRNA 模态间的系统性差异
cat("[Step 3d.2] Harmony batch correction on modality...\n")
co_embed <- RunHarmony(co_embed, group.by.vars = "modality",
                       reduction = "pca", dims.use = 1:30,
                       reduction.save = "harmony")

# UMAP 基于 Harmony 校正后的嵌入
co_embed <- RunUMAP(co_embed, reduction = "harmony", dims = 1:30)

step3_end <- Sys.time()
cat(sprintf("[计时] 阶段 3 耗时: %.1f 秒\n",
            as.numeric(difftime(step3_end, step3_start, units = "secs"))))

# 3e) 联合 UMAP 可视化
p_coembed_modality <- DimPlot(co_embed, group.by = "modality",
                               cols = c("scATAC" = "#E64B35", "scRNA" = "#4DBBD5"),
                               pt.size = 0.3, shuffle = TRUE) +
  ggtitle("Co-embedding (Harmony): scATAC vs scRNA")

p_coembed_celltype <- DimPlot(co_embed, group.by = "predicted.id",
                               pt.size = 0.3, label = TRUE, repel = TRUE,
                               label.size = 2.5) +
  ggtitle("Co-embedding (Harmony): Cell Types") + NoLegend()

ggsave(file.path(out_dir, "Fig_phase3_coembedding_modality.png"),
       plot = p_coembed_modality, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "Fig_phase3_coembedding_celltype.png"),
       plot = p_coembed_celltype, width = 8, height = 6, dpi = 300)
cat("已保存: Fig_phase3 co-embedding 图 (Harmony校正)\n")


# ╔═══════════════════════════════════════════════════════════╗
# ║  阶段 4：Label Transfer 评估与细胞类型校准               ║
# ╚═══════════════════════════════════════════════════════════╝
cat("\n\n")
cat("############################################\n")
cat("#  阶段 4：Label Transfer 评估             #\n")
cat("############################################\n")

# 4a) 对比：ATAC 原有 BioClassification vs RNA transfer 的 predicted.id
original_label    <- seurat_atac$BioClassification
transferred_label <- seurat_atac$predicted.id
prediction_score  <- seurat_atac$prediction.score.max

cat("[Transfer 置信度] median =", round(median(prediction_score), 3),
    " mean =", round(mean(prediction_score), 3),
    " min =", round(min(prediction_score), 3), "\n")

# 4b) 混淆矩阵
conf_mat <- table(Original = original_label, Transferred = transferred_label)
cat("\n混淆矩阵维度:", nrow(conf_mat), "x", ncol(conf_mat), "\n")

# 逐类型一致率（正确处理非方阵：按行名逐个查找对角元素）
orig_types <- rownames(conf_mat)
trans_types <- colnames(conf_mat)
shared_labels <- intersect(orig_types, trans_types)

agreement_per_type <- sapply(orig_types, function(ct) {
  if (ct %in% trans_types) {
    conf_mat[ct, ct] / sum(conf_mat[ct, ])
  } else {
    0  # 该类型在 RNA 中不存在，一致率为 0
  }
})

# 总体一致率：只统计行列标签匹配的对角元素之和
diag_sum <- sum(sapply(shared_labels, function(ct) conf_mat[ct, ct]))
overall_agreement <- diag_sum / sum(conf_mat)

cat("总体一致率:", round(overall_agreement * 100, 1), "%\n")
cat("(共有标签:", length(shared_labels), "/", length(orig_types), "原始类型)\n")

agreement_df <- data.frame(
  celltype   = orig_types,
  n_cells    = as.integer(rowSums(conf_mat)),
  agreement  = round(agreement_per_type * 100, 1)
) %>% arrange(desc(agreement))

cat("\n逐类型一致率:\n")
print(agreement_df, row.names = FALSE)

fwrite(agreement_df, file.path(out_dir, "label_transfer_agreement.csv"))
cat("已保存: label_transfer_agreement.csv\n")

# 4b.2) 高置信度子集一致率（prediction.score.max >= 0.5）
conf_threshold <- 0.5
high_conf_mask <- prediction_score >= conf_threshold
n_high_conf <- sum(high_conf_mask)
cat("\n[高置信度分析] threshold =", conf_threshold, "\n")
cat("  高置信度细胞:", n_high_conf, "/", length(high_conf_mask),
    "(", round(n_high_conf / length(high_conf_mask) * 100, 1), "%)\n")

if (n_high_conf > 0) {
  high_conf_agreement <- mean(original_label[high_conf_mask] ==
                              transferred_label[high_conf_mask])
  cat("  高置信度一致率:", round(high_conf_agreement * 100, 1), "%\n")
  cat("  （对比整体一致率:", round(overall_agreement * 100, 1), "%）\n")

  # 高置信度子集的逐类型一致率
  conf_mat_hc <- table(
    Original    = original_label[high_conf_mask],
    Transferred = transferred_label[high_conf_mask]
  )
  orig_types_hc <- rownames(conf_mat_hc)
  trans_types_hc <- colnames(conf_mat_hc)

  agreement_hc <- sapply(orig_types_hc, function(ct) {
    if (ct %in% trans_types_hc) conf_mat_hc[ct, ct] / sum(conf_mat_hc[ct, ]) else 0
  })

  agreement_hc_df <- data.frame(
    celltype       = orig_types_hc,
    n_cells_hc     = as.integer(rowSums(conf_mat_hc)),
    agreement_hc   = round(agreement_hc * 100, 1)
  ) %>% arrange(desc(agreement_hc))

  cat("\n  高置信度逐类型一致率:\n")
  print(agreement_hc_df, row.names = FALSE)
  fwrite(agreement_hc_df, file.path(out_dir, "label_transfer_agreement_highconf.csv"))
  cat("  已保存: label_transfer_agreement_highconf.csv\n")
}

# 4c) Transfer 置信度分布
p_score <- ggplot(data.frame(score = prediction_score), aes(score)) +
  geom_histogram(bins = 50, fill = "#3C5488", color = "white") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "Prediction Score Distribution (Label Transfer Confidence)",
       x = "Max Prediction Score", y = "Cells") +
  theme_classic()

ggsave(file.path(out_dir, "Fig_phase4_prediction_score.png"),
       plot = p_score, width = 6, height = 4, dpi = 300)
cat("已保存: Fig_phase4_prediction_score.png\n")


# ╔═══════════════════════════════════════════════════════════╗
# ║  阶段 5：TRS 调控解释（差异表达 + 通路富集）             ║
# ╚═══════════════════════════════════════════════════════════╝
cat("\n\n")
cat("############################################\n")
cat("#  阶段 5：TRS 调控解释                    #\n")
cat("############################################\n")

step5_start <- Sys.time()

# 5a) 加载 TRS 矩阵
trs_path <- file.path(data_dir, "scavenge_out_small_all_traits",
                       "TRS_matrix_small_all_traits_fullcells.rds")
trs_full <- readRDS(trs_path)
cat("TRS 矩阵:", nrow(trs_full), "cells x", ncol(trs_full) - 1, "traits\n")

# 5b) 将 TRS 合并到 ATAC Seurat metadata
trait_cols <- setdiff(colnames(trs_full), "cell")
trs_df <- as.data.frame(trs_full)
rownames(trs_df) <- trs_df$cell
trs_df$cell <- NULL

common_cells <- intersect(colnames(seurat_atac), rownames(trs_df))
cat("TRS 可匹配细胞:", length(common_cells), "/", ncol(seurat_atac), "\n")
seurat_atac <- AddMetaData(seurat_atac, metadata = trs_df[common_cells, , drop = FALSE])

# 5c) 对选定 trait 做 TRS-high vs TRS-low 差异表达
traits_to_analyze <- c("MONO", "LYM", "PLT")
traits_to_analyze <- intersect(traits_to_analyze, trait_cols)
cat("将分析的 traits:", paste(traits_to_analyze, collapse = ", "), "\n")

# DE 用 RNA_from_ATAC assay（即 Gene Activity，已归一化）
DefaultAssay(seurat_atac) <- "RNA_from_ATAC"

de_results_all <- list()

for (tr in traits_to_analyze) {
  cat("\n--- Trait:", tr, "---\n")

  trs_col <- make.names(tr)
  if (!(trs_col %in% colnames(seurat_atac@meta.data))) {
    cat("  跳过: TRS 列不存在\n")
    next
  }

  trs_values <- seurat_atac@meta.data[[trs_col]]
  q75 <- quantile(trs_values, 0.75, na.rm = TRUE)
  q25 <- quantile(trs_values, 0.25, na.rm = TRUE)

  seurat_atac@meta.data$trs_group <- "mid"
  seurat_atac@meta.data$trs_group[trs_values >= q75] <- "TRS_high"
  seurat_atac@meta.data$trs_group[trs_values <= q25] <- "TRS_low"

  n_high <- sum(seurat_atac$trs_group == "TRS_high")
  n_low  <- sum(seurat_atac$trs_group == "TRS_low")
  cat("  TRS_high:", n_high, " TRS_low:", n_low, "\n")

  Idents(seurat_atac) <- "trs_group"

  de <- FindMarkers(seurat_atac, ident.1 = "TRS_high", ident.2 = "TRS_low",
                    logfc.threshold = 0.1, min.pct = 0.05)

  de$gene <- rownames(de)
  de$trait <- tr
  de_results_all[[tr]] <- de

  n_up   <- sum(de$avg_log2FC > 0.25 & de$p_val_adj < 0.05)
  n_down <- sum(de$avg_log2FC < -0.25 & de$p_val_adj < 0.05)
  cat("  差异基因: up =", n_up, " down =", n_down, "(|log2FC|>0.25, padj<0.05)\n")

  fwrite(de, file.path(out_dir, paste0("DE_TRS_high_vs_low_", tr, ".csv")))
}

de_combined <- rbindlist(de_results_all, use.names = TRUE)
fwrite(de_combined, file.path(out_dir, "DE_TRS_all_traits_combined.csv"))
cat("\n已保存: DE_TRS_all_traits_combined.csv\n")

# 5d) GO 通路富集（如果 clusterProfiler 已安装）
run_go <- requireNamespace("clusterProfiler", quietly = TRUE) &&
          requireNamespace("org.Hs.eg.db", quietly = TRUE)

if (run_go) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  cat("\n[GO 富集分析]\n")

  go_results_all <- list()

  for (tr in traits_to_analyze) {
    de <- de_results_all[[tr]]
    if (is.null(de)) next

    up_genes <- de$gene[de$avg_log2FC > 0.25 & de$p_val_adj < 0.05]

    if (length(up_genes) < 5) {
      cat("  [", tr, "] 上调基因不足5个，跳过GO\n")
      next
    }

    cat("  [", tr, "] GO enrichment for", length(up_genes), "upregulated genes...\n")

    ego <- enrichGO(
      gene     = up_genes,
      OrgDb    = org.Hs.eg.db,
      keyType  = "SYMBOL",
      ont      = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1
    )

    if (nrow(ego) > 0) {
      ego_df <- as.data.frame(ego)
      ego_df$trait <- tr
      go_results_all[[tr]] <- ego_df
      fwrite(ego_df, file.path(out_dir, paste0("GO_BP_TRS_high_", tr, ".csv")))
      cat("    显著 GO terms:", nrow(ego_df), "\n")

      p_go <- barplot(ego, showCategory = 15) +
        ggtitle(paste0(tr, ": GO (BP) enrichment in TRS-high cells"))
      ggsave(file.path(out_dir, paste0("Fig_phase5_GO_", tr, ".png")),
             plot = p_go, width = 10, height = 7, dpi = 300)
    } else {
      cat("    无显著 GO term\n")
    }
  }

  if (length(go_results_all) > 0) {
    go_combined <- rbindlist(go_results_all, use.names = TRUE, fill = TRUE)
    fwrite(go_combined, file.path(out_dir, "GO_BP_all_traits_combined.csv"))
    cat("已保存: GO_BP_all_traits_combined.csv\n")
  }
} else {
  cat("\n[跳过 GO 富集] 请安装: BiocManager::install(c('clusterProfiler','org.Hs.eg.db'))\n")
}

step5_end <- Sys.time()
cat(sprintf("\n[计时] 阶段 5 耗时: %.1f 秒\n",
            as.numeric(difftime(step5_end, step5_start, units = "secs"))))


# ╔═══════════════════════════════════════════════════════════╗
# ║  阶段 6：交叉验证与汇总可视化                            ║
# ╚═══════════════════════════════════════════════════════════╝
cat("\n\n")
cat("############################################\n")
cat("#  阶段 6：交叉验证与汇总可视化            #\n")
cat("############################################\n")

# 6a) 已知 marker 基因验证
marker_sets <- list(
  MONO = c("CD14", "LYZ", "S100A8", "S100A9", "VCAN", "FCN1"),
  LYM  = c("CD3D", "CD3E", "CD79A", "CD79B", "NKG7", "GNLY"),
  PLT  = c("PF4", "PPBP", "ITGA2B", "GP9", "TUBB1", "TREML1")
)

cat("\n[Marker 基因验证]\n")
validation_results <- list()

for (tr in traits_to_analyze) {
  markers <- marker_sets[[tr]]
  if (is.null(markers)) next

  de <- de_results_all[[tr]]
  if (is.null(de)) next

  markers_in_de <- de[de$gene %in% markers, ]
  if (nrow(markers_in_de) > 0) {
    cat("\n  [", tr, "] Marker 基因在 DE 结果中:\n")
    print_df <- markers_in_de[, c("gene", "avg_log2FC", "p_val_adj")]
    print_df$avg_log2FC <- round(print_df$avg_log2FC, 3)
    print_df$p_val_adj  <- signif(print_df$p_val_adj, 3)
    print_df$direction  <- ifelse(print_df$avg_log2FC > 0, "UP in TRS-high", "DOWN in TRS-high")
    print(print_df, row.names = FALSE)

    n_correct <- sum(markers_in_de$avg_log2FC > 0)
    cat("    方向正确:", n_correct, "/", nrow(markers_in_de), "\n")

    validation_results[[tr]] <- data.frame(
      trait = tr,
      n_markers_tested = nrow(markers_in_de),
      n_correct_direction = n_correct,
      pct_correct = round(n_correct / nrow(markers_in_de) * 100, 1)
    )
  } else {
    cat("  [", tr, "] 无 marker 基因在 DE 结果中\n")
  }
}

if (length(validation_results) > 0) {
  validation_df <- rbindlist(validation_results)
  cat("\n[验证汇总]\n")
  print(validation_df, row.names = FALSE)
  fwrite(validation_df, file.path(out_dir, "marker_validation_summary.csv"))
  cat("已保存: marker_validation_summary.csv\n")
}

# 6b) Volcano plots
cat("\n[Volcano Plots]\n")
for (tr in traits_to_analyze) {
  de <- de_results_all[[tr]]
  if (is.null(de)) next

  de$neg_log10_padj <- -log10(de$p_val_adj + 1e-300)
  de$significance <- "NS"
  de$significance[de$avg_log2FC > 0.25 & de$p_val_adj < 0.05] <- "Up"
  de$significance[de$avg_log2FC < -0.25 & de$p_val_adj < 0.05] <- "Down"

  markers <- marker_sets[[tr]]
  de$label <- ifelse(de$gene %in% markers, de$gene, "")

  p_vol <- ggplot(de, aes(x = avg_log2FC, y = neg_log10_padj, color = significance)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey70")) +
    geom_text(aes(label = label), size = 3, color = "black", vjust = -0.5, hjust = 0.5) +
    labs(title = paste0(tr, ": TRS-high vs TRS-low (Gene Activity)"),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
    theme_classic() +
    theme(legend.position = "top")

  ggsave(file.path(out_dir, paste0("Fig_phase6_volcano_", tr, ".png")),
         plot = p_vol, width = 7, height = 6, dpi = 300)
}
cat("已保存: Volcano plots\n")

# 6c) TRS UMAP 可视化（在 co-embedding 上，只画 ATAC 细胞）
cat("\n[TRS UMAP]\n")
atac_cells_in_co <- intersect(colnames(seurat_atac), colnames(co_embed))

for (tr in traits_to_analyze) {
  trs_col <- make.names(tr)
  if (!(trs_col %in% colnames(seurat_atac@meta.data))) next

  trs_vec <- rep(NA_real_, ncol(co_embed))
  names(trs_vec) <- colnames(co_embed)
  trs_vec[atac_cells_in_co] <- seurat_atac@meta.data[atac_cells_in_co, trs_col]
  co_embed[[paste0("TRS_", tr)]] <- trs_vec

  p_trs <- FeaturePlot(co_embed, features = paste0("TRS_", tr), pt.size = 0.3,
                        order = TRUE) +
    scale_color_viridis_c(na.value = "grey90") +
    ggtitle(paste0("TRS (", tr, ") on Co-embedding"))

  ggsave(file.path(out_dir, paste0("Fig_phase6_TRS_umap_", tr, ".png")),
         plot = p_trs, width = 7, height = 6, dpi = 300)
}
cat("已保存: TRS UMAP plots\n")


# =========================
# 保存最终对象（不压缩，避免卡住）
# =========================
cat("\n========== 保存最终结果 ==========\n")
saveRDS(seurat_atac, file.path(out_dir, "seurat_atac_final.rds"), compress = FALSE)
saveRDS(co_embed,    file.path(out_dir, "seurat_coembedding.rds"), compress = FALSE)
cat("已保存: seurat_atac_final.rds\n")
cat("已保存: seurat_coembedding.rds\n")


# =========================
# 全局总结
# =========================
script_end <- Sys.time()
cat("\n")
cat("══════════════════════════════════════════════════════\n")
cat("  多组学融合阶段 3-6 全部完成\n")
cat("══════════════════════════════════════════════════════\n")
cat("  阶段 3: 跨模态整合 (pcaproject + co-embedding) done\n")
cat("  阶段 4: Label Transfer 一致率 =",
    round(overall_agreement * 100, 1), "%\n")
cat("  阶段 5: TRS 差异表达 + GO 通路富集            done\n")
cat("  阶段 6: Marker 验证 + Volcano + TRS UMAP      done\n")
cat("  总耗时:", round(as.numeric(difftime(script_end, script_start, units = "mins")), 1), "分钟\n")
cat("  输出目录:", out_dir, "\n")
cat("══════════════════════════════════════════════════════\n")
