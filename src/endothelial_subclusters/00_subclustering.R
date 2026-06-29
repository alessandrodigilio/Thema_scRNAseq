#######################################
### Endothelial cell subclustering ###
#######################################

# subset the annotated object to activated endothelial cells

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
source("src/endothelial_subclusters/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "annotated.rds")
out_dir <- file.path(data_dir, "integrated_object")
fig_dir <- file.path(figures_dir, "activated_endothelial_subclustering")
res_dir <- file.path(results_dir, "activated_endothelial_subclustering")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
output_subset_object <- file.path(out_dir, "activated_endothelial_subclustered.rds")
output_full_object <- file.path(out_dir, "annotated_endothelial_subclusters.rds")

# set parameters
parent_celltype <- "Activated endothelial cells"
batch_col <- "sample_id"
reduction_name <- "umap.harmony.endothelial"
cluster_col <- "endothelial_subcluster"

set.seed(1234)

# load annotated object
cat("Loading annotated object...\n")
obj_full <- readRDS(input_object)
cat("Cells:", ncol(obj_full), "\n")

DefaultAssay(obj_full) <- "RNA"
obj_full <- JoinLayers(obj_full, assay = "RNA")
obj_full <- NormalizeData(obj_full, assay = "RNA", verbose = FALSE)

# subset to activated endothelial cells
cells_use <- rownames(obj_full@meta.data)[obj_full$cell_type == parent_celltype]
obj <- subset(obj_full, cells = cells_use)
cat("Activated endothelial subset cells:", ncol(obj), "\n")
print(table(obj$cell_type))

# recompute subset embedding and clustering
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 20, verbose = FALSE)

# run Harmony integration to correct for batch effects
cat("\nRunning Harmony on endothelial subset...\n")
obj <- RunHarmony(
  object = obj,
  group.by.vars = batch_col,
  reduction = "pca",
  dims.use = 1:20,
  reduction.save = "harmony",
  verbose = FALSE
)

# run UMAP on the integrated data
obj <- RunUMAP(
  object = obj,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = reduction_name,
  reduction.key = "UMAPENDO_",
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  verbose = FALSE
)

# compute KNN and clustering (range of res)
obj <- FindNeighbors(
  object = obj,
  reduction = "harmony",
  dims = 1:20,
  k.param = 20,
  verbose = FALSE
)
# summary
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
# save
write.csv(cluster_summary, file.path(res_dir, "activated_endothelial_cluster_resolution_summary.csv"), row.names = FALSE)

# set the 0.3 res and update the Seurat object
selected_cluster_col <- paste0("RNA_snn_res.", 0.3)
obj[[cluster_col]] <- as.character(obj@meta.data[[selected_cluster_col]])
cluster_levels <- sort_cluster_levels(obj@meta.data[[cluster_col]])
Idents(obj) <- factor(obj@meta.data[[cluster_col]], levels = cluster_levels)

# marker finding for endothelial subclusters
cat("\nRunning FindAllMarkers on endothelial subclusters...\n")
markers <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)
# save
cat("Total marker genes found:", nrow(markers), "\n")
write.csv(markers, file.path(res_dir, "all_markers_activated_endothelial_subclusters.csv"), row.names = FALSE)

# top markers
top_markers_save <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()
# save
write.csv(top_markers_save, file.path(res_dir, "top_markers_per_activated_endothelial_subcluster.csv"), row.names = FALSE)
top_markers_review <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) %>%
  summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
top_markers_review$cluster <- as.character(top_markers_review$cluster)
cat("\nTop markers used for endothelial subcluster review:\n")
print(top_markers_review)

# save subset object
saveRDS(obj, output_subset_object)

# add endothelial subcluster layer back to the full annotated object
obj_full$endothelial_subcluster <- NA_character_
obj_full$endothelial_subcluster[Cells(obj)] <- as.character(obj$endothelial_subcluster)
obj_full$endothelial_parent_type <- NA_character_
obj_full$endothelial_parent_type[Cells(obj)] <- as.character(obj$cell_type)
saveRDS(obj_full, output_full_object)

# dotplot of top subcluster markers
review_features <- top_markers_review$top_markers
review_features <- unique(unlist(strsplit(review_features, ", ", fixed = TRUE)))
review_features <- review_features[review_features %in% rownames(obj)]

# plot dotplot
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
    plot.margin = margin(20, 35, 35, 25)
  )
# save
ggsave(file.path(fig_dir, "dotplot_top_markers_activated_endothelial_subclusters.png"), p_marker_dot, width = 12, height = 6, dpi = 600)

# dotplot focused on iron-related genes
iron_features <- endothelial_iron_features[endothelial_iron_features %in% rownames(obj)]
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
    plot.margin = margin(20, 35, 35, 25)
  )
# save
ggsave(file.path(fig_dir, "dotplot_iron_related_markers_activated_endothelial_subclusters.png"), p_iron_dot, width = 12, height = 6, dpi = 600)

cat("\n============================================================\n")
cat("Endothelial subclustering complete.\n")
cat("Subset object : ", output_subset_object, "\n", sep = "")
cat("Updated object: ", output_full_object, "\n", sep = "")
cat("Results dir   : ", res_dir, "\n", sep = "")
cat("Figures dir   : ", fig_dir, "\n", sep = "")
cat("============================================================\n")
