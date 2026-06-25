#########################################
### Sample integration and clustering ###
#########################################

# load the merged filtered RNA object, integrate RNA across samples,
# inspect clustering stability with clustree, choose one Leiden resolution,
# and save the integrated object.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(clustree)
  library(future)
})

setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")

# directories
input_object <- file.path(filtered_data_dir, "filtered_samples.rds")
out_dir <- integrated_data_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(figures_dir, "integration")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
results_int_dir <- file.path(results_dir, "integration")
dir.create(results_int_dir, recursive = TRUE, showWarnings = FALSE)
output_object <- file.path(out_dir, "integrated.rds")

# reproducibility
set.seed(1234)
plan("sequential") # avoid SCTransform failures 
options(future.globals.maxSize = 8 * 1024^3)

# load object
cat("Loading filtered RNA object...\n")
obj <- readRDS(input_object)
obj[[batch_var]] <- factor(obj[[batch_var]][, 1])

# check object contents
cat("Cells:", ncol(obj), "\n")
cat("Samples:", paste(levels(obj[[batch_var]][, 1]), collapse = ", "), "\n")
cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")

# run SCTransform and PCA while excluding ribosomal genes only from the PCA input
cat("\nRNA preprocessing (SCTransform + PCA)...\n")
DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE)
ribo_genes <- grep("^(RPS|RPL)", rownames(obj), value = TRUE)
pca_features_before <- VariableFeatures(obj)
pca_features <- setdiff(pca_features_before, ribo_genes)

cat("Variable features before ribosomal exclusion:", length(pca_features_before), "\n")
cat("Ribosomal genes among variable features:", length(intersect(pca_features_before, ribo_genes)), "\n")
cat("Variable features after ribosomal exclusion:", length(pca_features), "\n")

# run PCA and UMAP on the SCTransform assay, excluding ribosomal genes from the PCA input
obj <- RunPCA(
  object = obj,
  assay = "SCT",
  features = pca_features,
  npcs = n_pcs,
  reduction.name = "pca",
  verbose = FALSE
)

obj <- RunUMAP(
  object = obj,
  reduction = "pca",
  dims = harmony_dims,
  reduction.name = "rna.umap",
  reduction.key = "rnaUMAP_",
  umap.method = "uwot",
  metric = "cosine",
  n.neighbors = umap_neighbors,
  min.dist = umap_min_dist,
  spread = umap_spread,
  verbose = FALSE
)

p_rna_pre <- DimPlot(
  object = obj,
  reduction = "rna.umap",
  group.by = batch_var,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("RNA pre-Harmony batch correction")

ggsave(
  filename = file.path(fig_dir, "RNA_preHarmony_by_sample.png"),
  plot = p_rna_pre,
  width = 8,
  height = 8,
  dpi = 600
)

# correct RNA embeddings across samples with Harmony
cat("Harmony on RNA (PCA)...\n")
obj <- RunHarmony(
  object = obj,
  group.by.vars = batch_var,
  reduction = "pca",
  dims.use = harmony_dims,
  reduction.save = "harmony",
  verbose = FALSE
)

obj <- RunUMAP(
  object = obj,
  reduction = "harmony",
  dims = harmony_dims,
  reduction.name = "umap.harmony.rna",
  reduction.key = "hRNAUMAP_",
  umap.method = "uwot",
  metric = "cosine",
  n.neighbors = umap_neighbors,
  min.dist = umap_min_dist,
  spread = umap_spread,
  verbose = FALSE
)

p_rna_post <- DimPlot(
  object = obj,
  reduction = "umap.harmony.rna",
  group.by = batch_var,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("RNA post-Harmony batch correction")

ggsave(
  filename = file.path(fig_dir, "RNA_postHarmony_by_sample.png"),
  plot = p_rna_post,
  width = 8,
  height = 8,
  dpi = 600
)

# build the shared nearest-neighbor graph on the Harmony RNA space
cat("\nBuilding RNA graph on Harmony embeddings...\n")
DefaultAssay(obj) <- "SCT"

obj <- FindNeighbors(
  object = obj,
  reduction = "harmony",
  dims = harmony_dims,
  k.param = graph_k_param,
  graph.name = c("RNA_nn", "RNA_snn"),
  verbose = FALSE
)

# run Leiden clustering on a short grid of resolutions.
# use clustree to find stable splits then keep one final resolution 
cat("\nScreening Leiden resolutions on the RNA SNN graph...\n")

cluster_summary <- vector("list", length(cluster_res_grid))

# run FindClusters across the grid of resolutions
for (i in seq_along(cluster_res_grid)) {
  resolution <- cluster_res_grid[[i]]
  cluster_col <- res_to_colname(resolution) # helper in utils.R

  obj <- FindClusters(
    object = obj,
    graph.name = "RNA_snn",
    algorithm = 4,
    resolution = resolution,
    cluster.name = cluster_col,
    random.seed = 1,
    verbose = FALSE
  )

  n_clusters <- length(unique(obj@meta.data[[cluster_col]]))
  cat(sprintf("resolution = %.1f -> %d clusters\n", resolution, n_clusters))

  cluster_summary[[i]] <- data.frame(
    resolution = resolution,
    cluster_col = cluster_col,
    n_clusters = n_clusters
  )
}

# save the cluster summary table
cluster_summary <- do.call(rbind, cluster_summary)
write.csv(cluster_summary, file.path(results_int_dir, "cluster_resolution_summary.csv"), row.names = FALSE)

# clustree used after running FindClusters across multiple resolutions
p_tree <- clustree(
  obj@meta.data,
  prefix = "snn_res."
) +
  ggtitle("Clustering tree across RNA Leiden resolutions")

ggsave(
  filename = file.path(fig_dir, "RNA_clustree.png"),
  plot = p_tree,
  width = 10,
  height = 8,
  dpi = 600
)

cl_col <- res_to_colname(selected_cluster_res)
obj$cluster <- obj@meta.data[[cl_col]]
Idents(obj) <- "cluster"

n_clust <- length(unique(obj$cluster))
cat(sprintf("\nSelected resolution = %.1f -> %d clusters\n", selected_cluster_res, n_clust))

p_harmony_cluster <- DimPlot(
  object = obj,
  reduction = "umap.harmony.rna",
  group.by = "cluster",
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("Harmony RNA Leiden clusters")

ggsave(
  filename = file.path(fig_dir, "Harmony_RNA_leiden_clusters.png"),
  plot = p_harmony_cluster,
  width = 8,
  height = 8,
  dpi = 600
)

p_panel <- (p_rna_pre + p_rna_post + p_harmony_cluster) &
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(fig_dir, "RNA_integration_comparison.png"),
  plot = p_panel,
  width = 18,
  height = 6,
  dpi = 600
)

# save QC diagnostics on the integrated embedding
cat("Saving QC plots...\n")

qc_feats <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "soupx_removed_fraction")

p_qc <- FeaturePlot(
  object = obj,
  features = qc_feats,
  reduction = "umap.harmony.rna",
  ncol = 2,
  order = TRUE
)

ggsave(
  filename = file.path(fig_dir, "QC_on_Harmony_RNA.png"),
  plot = p_qc,
  width = 12,
  height = 10,
  dpi = 600
)

# save the integrated object with all screened resolutions with 'cluster' column
cat("\nSaving integrated object...\n")
saveRDS(obj, output_object)

cat("\n============================================================\n")
cat("Integration complete.\n")
cat("Object saved        : ", output_object, "\n")
cat("Figures dir         : ", fig_dir, "\n")
cat("Resolution summary  : ", file.path(results_int_dir, "cluster_resolution_summary.csv"), "\n")
cat("Selected resolution : ", selected_cluster_res, "\n")
cat("Final clusters col  : cluster\n")
cat("N clusters          : ", n_clust, "\n")
cat("============================================================\n")
