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

trait_import <- example_data(name="mono.PP001.bed")
SE_pbmc5k <- example_data(name="pbmc5k_SE.rda")

#gchromVAR analysis
SE_pbmc5k <- addGCBias(object = SE_pbmc5k, 
                       genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
SE_pbmc5k_bg <- getBackgroundPeaks(object = SE_pbmc5k, 
                                   niterations = 200)
SE_pbmc5k_DEV <- computeWeightedDeviations(object = SE_pbmc5k, 
                                           weights = trait_import, 
                                           background_peaks = SE_pbmc5k_bg)

z_score_mat <- data.frame(SummarizedExperiment::colData(SE_pbmc5k), 
                          z_score=t(SummarizedExperiment::assays(SE_pbmc5k_DEV)[["z"]]) |> c())
head(z_score_mat)

#生成种子细胞索引（如果太多就只使用前5%）
seed_idx <- seedindex(z_score = z_score_mat$z_score,
                      percent_cut = 0.05)

#计算规模因子
scale_factor <- cal_scalefactor(z_score = z_score_mat$z_score, 
                                percent_cut = 0.01)


#构建m—knn图
#Calculate tfidf-mat
peak_by_cell_mat <- SummarizedExperiment::assay(SE_pbmc5k)
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, 
                   mat_binary=TRUE, 
                   TF=TRUE, 
                   log_TF=TRUE)

#Calculate lsi-mat
lsi_mat <- do_lsi(mat = tfidf_mat, 
                  dims = 30)

#Calculate m-knn graph
mutualknn30 <- getmutualknn(lsimat = lsi_mat, 
                            num_k = 30)


#Network propagation网络传播
np_score <- randomWalk_sparse(intM = mutualknn30, 
                              queryCells = rownames(mutualknn30)[seed_idx], 
                              gamma = 0.05)

omit_idx <- np_score==0
sum(omit_idx)

mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- capOutlierQuantile(x = np_score, 
                          q_ceiling = 0.95) |> max_min_scale()
TRS <- TRS * scale_factor
mono_mat <- data.frame(z_score_mat[!omit_idx, ], 
                       seed_idx[!omit_idx], 
                       np_score, 
                       TRS)
head(mono_mat)


#UMAP plots of cell type annotation and cell-to-cell graph单元格类型注释绘图
#Cell type annotation单元格类型注释
p <- ggplot(data=mono_mat, aes(x, y, color=color)) + 
  geom_point(size=1, na.rm = TRUE) + 
  pretty_plot() + 
  theme(legend.title = element_blank()) + 
  labs(x="UMAP 1",y="UMAP 2")
p

#Visualize cell-to-cell graph if you have low-dimensional coordinates such as UMAP1 and UMAP2
mutualknn30_graph <- graph_from_adjacency_matrix(adjmatrix = mutualknn30, 
                                                 mode = "undirected", 
                                                 diag = FALSE)
igraph::plot.igraph(x = mutualknn30_graph, 
                    vertex.size=0.8, 
                    vertex.label=NA, 
                    vertex.color=adjustcolor("#c7ce3d", alpha.f = 1), 
                    vertex.frame.color=NA, 
                    edge.color=adjustcolor("#443dce", alpha.f = 1), 
                    edge.width=0.3, edge.curved=.5, 
                    layout=as.matrix(data.frame(mono_mat$x, mono_mat$y)))


#Comparsion before and after SCAVENGE analysis
#Z score based visualization Z评分可视化
p1 <- ggplot(data=mono_mat, aes(x, y, color=z_score)) + 
  geom_point(size=1, na.rm = TRUE, alpha = 0.6) + 
  scale_color_viridis_c() + 
  scale_alpha() +
  pretty_plot() + 
  theme(legend.title = element_blank()) + 
  labs(x="UMAP 1", y="UMAP 2")
p1

#柱状图
pp1 <- ggplot(data=mono_mat,  aes(x=color, y=z_score))  +
  geom_boxplot(aes(fill=color, color=color), outlier.shape=NA) + 
  guides(fill=FALSE) + 
  pretty_plot(fontsize = 10) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", 
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp1

#SCAVENGE TRS based visualization TRS可视化
p2 <- ggplot(data=mono_mat, aes(x, y, color=TRS)) + 
  geom_point(size=1, na.rm = TRUE, alpha = 0.6) + 
  scale_color_viridis_c() + 
  scale_alpha() + pretty_plot() + 
  theme(legend.title = element_blank()) + 
  labs(x="UMAP 1", y="UMAP 2")
p2
#柱状图
pp2 <- ggplot(data=mono_mat,  aes(x=color, y=TRS))  +
  geom_boxplot(aes(fill=color, color=color), outlier.shape=NA) + 
  guides(fill=FALSE) + 
  pretty_plot(fontsize = 10) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + theme(legend.position = "none")
pp2


#Trait relevant cell determination from permutation test从排列结果中确定细胞特性
#About 2 mins
mono_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, 
                                 seed_idx=mono_mat$seed_idx, 
                                 topseed_npscore=mono_mat$np_score, 
                                 permutation_times=100, # Increase to >=1000 in practice
                                 true_cell_significance=0.05, 
                                 rda_output=FALSE, 
                                 # mycores=8,# Increase in practice
                                 rw_gamma=0.05)
mono_mat2 <- data.frame(mono_mat, mono_permu)


#Look at the distribution of statistically significant phenotypically enriched and depleted cells
#观察统计学显著的表型富集和耗竭细胞的分布
#Enriched cells
mono_mat2 |>
  dplyr::group_by(color) |> 
  dplyr::summarise(enriched_cell=sum(true_cell_top_idx)) |> 
  ggplot(aes(x=color, y=enriched_cell, fill=color)) + 
  geom_bar(stat="identity") + 
  theme_classic()
#Depleted cells
mono_mat2$rev_true_cell_top_idx <- !mono_mat2$true_cell_top_idx
mono_mat2 |>
  dplyr::group_by(color) |> 
  dplyr::summarise(depleted_cell=sum(rev_true_cell_top_idx)) |> 
  ggplot(aes(x=color, y=depleted_cell, fill=color)) + 
  geom_bar(stat="identity") + 
  theme_classic()

