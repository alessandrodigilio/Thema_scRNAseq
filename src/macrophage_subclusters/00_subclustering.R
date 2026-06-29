################################
### Macrophage subclustering ###
################################

# subset the annotated object to the main macrophage populations

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(harmony)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")
source("src/macrophage_subclusters/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "annotated.rds")
out_dir <- file.path(data_dir, "integrated_object")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(figures_dir, "macrophage_subclustering")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
res_dir <- file.path(results_dir, "macrophage_subclustering")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

output_subset_object <- file.path(out_dir, "macrophages_subclustered.rds")
output_full_object <- file.path(out_dir, "annotated_macrophage_subclusters.rds")

# set parameters
parent_celltypes <- c(
  "Inflammatory macrophages (IL1B+)",
  "Resident macrophages (C1QC+)"
)
group_col <- "condition"
batch_col <- "sample_id"
reduction_name <- "umap.harmony.macrophage"
cluster_col <- "macrophage_subcluster"

set.seed(1234)

# load annotated object
if (!file.exists(input_object)) {
  stop("Missing annotated object: ", input_object)
}

cat("Loading annotated object...\n")
obj_full <- readRDS(input_object)
cat("Cells:", ncol(obj_full), "\n")

if (!"cell_type" %in% colnames(obj_full@meta.data)) {
  stop("cell_type not found in metadata. Run final annotation first.")
}

DefaultAssay(obj_full) <- "RNA"
obj_full <- JoinLayers(obj_full, assay = "RNA")

rna_data_layer <- tryCatch(
  LayerData(obj_full[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  cat("Normalizing RNA assay...\n")
  obj_full <- NormalizeData(obj_full, assay = "RNA", verbose = FALSE)
}

# subset to the main macrophage populations only
cells_use <- rownames(obj_full@meta.data)[obj_full$cell_type %in% parent_celltypes]

if (length(cells_use) == 0) {
  stop("No cells found for the selected parent macrophage populations")
}

obj <- subset(obj_full, cells = cells_use)
cat("Macrophage subset cells:", ncol(obj), "\n")
print(table(obj$cell_type))

# recompute subset embedding and clustering
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 20, verbose = FALSE)

if (!batch_col %in% colnames(obj@meta.data)) {
  stop("Batch column not found in metadata: ", batch_col)
}

cat("\nRunning Harmony on macrophage subset...\n")
obj <- RunHarmony(
  object = obj,
  group.by.vars = batch_col,
  reduction = "pca",
  dims.use = 1:20,
  reduction.save = "harmony",
  verbose = FALSE
)

obj <- RunUMAP(
  object = obj,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = reduction_name,
  reduction.key = "UMAPMAC_",
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  verbose = FALSE
)

obj <- FindNeighbors(
  object = obj,
  reduction = "harmony",
  dims = 1:20,
  k.param = 20,
  verbose = FALSE
)

cluster_summary <- data.frame(
  resolution = seq(0.2, 0.8, by = 0.1),
  n_clusters = NA_integer_,
  stringsAsFactors = FALSE
)

# run clustering for each resolution and summarize the number of clusters
for (i in seq_along(cluster_summary$resolution)) {
  res_here <- cluster_summary$resolution[i]
  obj <- FindClusters(
    object = obj,
    resolution = res_here,
    algorithm = 1,
    random.seed = 1234,
    verbose = FALSE
  )
  cluster_col_here <- paste0("RNA_snn_res.", res_here)
  cluster_summary$n_clusters[i] <- length(unique(obj@meta.data[[cluster_col_here]]))
}
# save res summary to decide
write.csv(cluster_summary, file.path(res_dir, "macrophage_cluster_resolution_summary.csv"), row.names = FALSE)

# set the selected clustering resolution and update the Seurat object
selected_cluster_col <- paste0("RNA_snn_res.", 0.3)
if (!selected_cluster_col %in% colnames(obj@meta.data)) {
  stop("Selected clustering column not found: ", selected_cluster_col)
}

# update the cluster column
obj[[cluster_col]] <- as.character(obj@meta.data[[selected_cluster_col]])
cluster_levels <- sort_cluster_levels(obj@meta.data[[cluster_col]])
Idents(obj) <- factor(obj@meta.data[[cluster_col]], levels = cluster_levels)

# marker finding for macrophage subclusters
cat("\nRunning FindAllMarkers on macrophage subclusters...\n")
markers <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

# print summary of marker genes found
cat("Total marker genes found:", nrow(markers), "\n")
# save all markers
write.csv(markers, file.path(res_dir, "all_markers_macrophage_subclusters.csv"), row.names = FALSE)
# top markers
top_markers_save <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>% # top 5 markers per cluster
  ungroup()
# save top 5 markers 
write.csv(top_markers_save, file.path(res_dir, "top_markers_per_macrophage_subcluster.csv"), row.names = FALSE)
top_markers_review <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) %>%
  summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
top_markers_review$cluster <- as.character(top_markers_review$cluster)
cat("\nTop markers used for macrophage subcluster review:\n")
print(top_markers_review)

# save subset object
saveRDS(obj, output_subset_object)

# add subcluster layer back to the full annotated object
obj_full$macrophage_subcluster <- NA_character_
obj_full$macrophage_subcluster[Cells(obj)] <- as.character(obj$macrophage_subcluster)
obj_full$macrophage_parent_type <- NA_character_
obj_full$macrophage_parent_type[Cells(obj)] <- as.character(obj$cell_type)
saveRDS(obj_full, output_full_object) # updated full object with macrophage subcluster layer

# dotplot of top subcluster markers
review_features <- top_markers_review$top_markers
review_features <- unique(unlist(strsplit(review_features, ", ", fixed = TRUE)))
review_features <- review_features[review_features %in% rownames(obj)]

if (length(review_features) > 0) {
  p_marker_dot <- DotPlot(
    object = obj,
    features = review_features,
    cols = c("white", "#d84b1c"),
    dot.scale = 5,
    col.min = 0,
    col.max = 3
  ) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11),
      legend.position = "top",
      legend.key.width = grid::unit(0.25, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.08, "cm"),
      plot.margin = margin(20, 35, 35, 25)
    )

  ggsave(
    file.path(fig_dir, "dotplot_top_markers_macrophage_subclusters.png"),
    p_marker_dot,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# dotplot focused on iron-related genes
iron_features <- macrophage_iron_features[macrophage_iron_features %in% rownames(obj)] # defined in global_config.R

if (length(iron_features) > 0) {
  p_iron_dot <- DotPlot(
    object = obj,
    features = iron_features,
    cols = c("white", "#7A1F2B"),
    dot.scale = 5,
    col.min = 0,
    col.max = 3
  ) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11),
      legend.position = "top",
      legend.key.width = grid::unit(0.25, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.08, "cm"),
      plot.margin = margin(20, 35, 35, 25)
    )

  ggsave(
    file.path(fig_dir, "dotplot_iron_related_markers_macrophage_subclusters.png"),
    p_iron_dot,
    width = 12,
    height = 6,
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Macrophage subclustering complete.\n")
cat("Subset object : ", output_subset_object, "\n")
cat("Updated object: ", output_full_object, "\n")
cat("Results dir   : ", res_dir, "\n")
cat("Figures dir   : ", fig_dir, "\n")
cat("============================================================\n")
