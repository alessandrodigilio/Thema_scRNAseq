# ---------------------------------------
# 02_GEX_integration.R
# Integrate filtered GEX samples
# ---------------------------------------

# Load the merged filtered GEX object, integrate RNA across samples,
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

setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
INPUT_OBJECT <- file.path(FILTERED_DATA_DIR, "gex_filtered_samples.rds")

OUT_DIR <- INTEGRATED_DATA_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(FIGURES_DIR, "integration")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RESULTS_INT_DIR <- file.path(RESULTS_DIR, "integration")
dir.create(RESULTS_INT_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_OBJECT <- file.path(OUT_DIR, "gex_integrated.rds")

# Set seed for reproducibility
set.seed(1234)

# Avoid SCTransform failures caused by future exporting large globals.
future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

# Convert numeric resolution to a stable metadata column name
res_to_colname <- function(resolution) {
  paste0("snn_res.", format(resolution, nsmall = 1, trim = TRUE))
}

# Load the filtered GEX object produced by the previous step
if (!file.exists(INPUT_OBJECT)) stop("Missing filtered object: ", INPUT_OBJECT)

cat("Loading filtered GEX object...\n")
obj <- readRDS(INPUT_OBJECT)

if (!BATCH_VAR %in% colnames(obj@meta.data)) {
  stop("Batch variable not found in metadata: ", BATCH_VAR)
}

obj[[BATCH_VAR]] <- factor(obj[[BATCH_VAR]][, 1])

cat("Cells:", ncol(obj), "\n")
cat("Samples:", paste(levels(obj[[BATCH_VAR]][, 1]), collapse = ", "), "\n")
cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")

#################
# ----- RNA -----#
#################

# Run SCTransform and PCA while excluding ribosomal genes only from
# the PCA input, keeping the assay otherwise unchanged.
cat("\nRNA preprocessing (SCTransform + PCA)...\n")
DefaultAssay(obj) <- "RNA"

if (!"percent.mt" %in% colnames(obj@meta.data)) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
}

obj <- SCTransform(obj, verbose = FALSE)

ribo_genes <- grep("^(RPS|RPL)", rownames(obj), value = TRUE)
pca_features_before <- VariableFeatures(obj)
pca_features <- setdiff(pca_features_before, ribo_genes)

cat("Variable features before ribosomal exclusion:", length(pca_features_before), "\n")
cat("Ribosomal genes among variable features:", length(intersect(pca_features_before, ribo_genes)), "\n")
cat("Variable features after ribosomal exclusion:", length(pca_features), "\n")

if (length(pca_features) < 50) {
  stop("Too few PCA features remain after excluding ribosomal genes")
}

obj <- RunPCA(
  object = obj,
  assay = "SCT",
  features = pca_features,
  npcs = N_PCS,
  reduction.name = "pca",
  verbose = FALSE
)

obj <- RunUMAP(
  object = obj,
  reduction = "pca",
  dims = HARMONY_DIMS,
  reduction.name = "rna.umap",
  reduction.key = "rnaUMAP_",
  umap.method = "uwot",
  metric = "cosine",
  n.neighbors = UMAP_NEIGHBORS,
  min.dist = UMAP_MIN_DIST,
  spread = UMAP_SPREAD,
  verbose = FALSE
)

p_rna_pre <- DimPlot(
  object = obj,
  reduction = "rna.umap",
  group.by = BATCH_VAR,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("RNA pre-Harmony batch correction")

ggsave(
  filename = file.path(FIG_DIR, "RNA_preHarmony_by_sample.png"),
  plot = p_rna_pre,
  width = 8,
  height = 8,
  dpi = 600
)

# Correct RNA embeddings across samples with Harmony
cat("Harmony on RNA (PCA)...\n")
obj <- harmony::RunHarmony(
  object = obj,
  group.by.vars = BATCH_VAR,
  reduction = "pca",
  dims.use = HARMONY_DIMS,
  reduction.save = "harmony",
  verbose = FALSE
)

obj <- RunUMAP(
  object = obj,
  reduction = "harmony",
  dims = HARMONY_DIMS,
  reduction.name = "umap.harmony.rna",
  reduction.key = "hRNAUMAP_",
  umap.method = "uwot",
  metric = "cosine",
  n.neighbors = UMAP_NEIGHBORS,
  min.dist = UMAP_MIN_DIST,
  spread = UMAP_SPREAD,
  verbose = FALSE
)

p_rna_post <- DimPlot(
  object = obj,
  reduction = "umap.harmony.rna",
  group.by = BATCH_VAR,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("RNA post-Harmony batch correction")

ggsave(
  filename = file.path(FIG_DIR, "RNA_postHarmony_by_sample.png"),
  plot = p_rna_post,
  width = 8,
  height = 8,
  dpi = 600
)

##################################
# ----- Resolution screening -----#
##################################

# Build the shared nearest-neighbor graph on the Harmony RNA space.
cat("\nBuilding RNA graph on Harmony embeddings...\n")
DefaultAssay(obj) <- "SCT"

obj <- FindNeighbors(
  object = obj,
  reduction = "harmony",
  dims = HARMONY_DIMS,
  k.param = GRAPH_K_PARAM,
  graph.name = c("RNA_nn", "RNA_snn"),
  verbose = FALSE
)

# Run Leiden clustering on a short grid of resolutions. The idea is:
# first inspect the clustree to find stable splits, then keep one
# final resolution for the rest of the workflow.
cat("\nScreening Leiden resolutions on the RNA SNN graph...\n")

cluster_summary <- vector("list", length(CLUSTER_RES_GRID))

for (i in seq_along(CLUSTER_RES_GRID)) {
  resolution <- CLUSTER_RES_GRID[[i]]
  cluster_col <- res_to_colname(resolution)

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

cluster_summary <- do.call(rbind, cluster_summary)
write.csv(
  cluster_summary,
  file.path(RESULTS_INT_DIR, "cluster_resolution_summary.csv"),
  row.names = FALSE
)

# clustree is used after the graph is built and after running
# FindClusters across multiple resolutions.
p_tree <- clustree::clustree(
  obj@meta.data,
  prefix = "snn_res."
) +
  ggplot2::ggtitle("Clustering tree across RNA Leiden resolutions")

ggplot2::ggsave(
  filename = file.path(FIG_DIR, "RNA_clustree.png"),
  plot = p_tree,
  width = 10,
  height = 8,
  dpi = 600
)

cl_col <- res_to_colname(SELECTED_CLUSTER_RES)
if (!cl_col %in% colnames(obj@meta.data)) {
  stop("Selected clustering column not found in metadata: ", cl_col)
}

obj$gex_cluster <- obj@meta.data[[cl_col]]
Idents(obj) <- "gex_cluster"

n_clust <- length(unique(obj$gex_cluster))
cat(sprintf("\nSelected resolution = %.1f -> %d clusters\n", SELECTED_CLUSTER_RES, n_clust))

#############################
# ----- Final embedding -----#
#############################

p_harmony_cluster <- DimPlot(
  object = obj,
  reduction = "umap.harmony.rna",
  group.by = "gex_cluster",
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) + ggtitle("Harmony RNA Leiden clusters")

ggsave(
  filename = file.path(FIG_DIR, "Harmony_RNA_leiden_clusters.png"),
  plot = p_harmony_cluster,
  width = 8,
  height = 8,
  dpi = 600
)

########################################
# Comparison pre/post Harmony clusters #
########################################

p_panel <- (p_rna_pre + p_rna_post + p_harmony_cluster) &
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(FIG_DIR, "RNA_integration_comparison.png"),
  plot = p_panel,
  width = 18,
  height = 6,
  dpi = 600
)

# Save QC diagnostics on the integrated embedding.
cat("Saving QC plots...\n")

qc_feats <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "soupx_removed_fraction"),
  colnames(obj@meta.data)
)

if (length(qc_feats) > 0) {
  p_qc <- FeaturePlot(
    object = obj,
    features = qc_feats,
    reduction = "umap.harmony.rna",
    ncol = 2,
    order = TRUE
  )

  ggsave(
    filename = file.path(FIG_DIR, "QC_on_Harmony_RNA.png"),
    plot = p_qc,
    width = 12,
    height = 10,
    dpi = 600
  )
}

# Save the integrated object with all screened resolutions plus the
# final column 'gex_cluster' used downstream.
cat("\nSaving integrated object...\n")
saveRDS(obj, OUTPUT_OBJECT)

cat("\n============================================================\n")
cat("Integration complete.\n")
cat("Object saved        : ", OUTPUT_OBJECT, "\n")
cat("Figures dir         : ", FIG_DIR, "\n")
cat("Resolution summary  : ", file.path(RESULTS_INT_DIR, "cluster_resolution_summary.csv"), "\n")
cat("Selected resolution : ", SELECTED_CLUSTER_RES, "\n")
cat("Final clusters col  : gex_cluster\n")
cat("N clusters          : ", n_clust, "\n")
cat("============================================================\n")
