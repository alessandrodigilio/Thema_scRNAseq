###############################################
### Cluster marker review and final labels ###
###############################################

# compute cluster markers, assign final cell type labels and save the
# annotated object used by the downstream analyses.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "integrated.rds")
output_object <- file.path(data_dir, "integrated_object", "annotated.rds")

fig_dir <- file.path(figures_dir, "annotation")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
res_dir <- file.path(results_dir, "annotation")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
cluster_col <- "cluster"
reduction_name <- "umap.harmony.rna"
min_pct <- 0.25
logfc_thr <- 0.25
n_save <- 5
n_review_markers <- 3
pt_size <- 0.01

# this list of markers will be helpful to better characterize
# the final cell type labels in the feature plot panel
feature_markers <- c(
  "PRG4", "CLIC5", "HBEGF",
  "MFAP5", "CXCL12", "SFRP4",
  "MMP3", "MMP1", "INHBA",
  "IL1B", "C1QC", "FOLR2",
  "SELE", "ACKR1", "EMCN",
  "RGS5", "ACTA2", "MYH11",
  "CD1C", "FCER1A", "CLEC10A",
  "HMOX1", "SLC40A1"
)
feature_ncol <- 4 # number of columns in the featureplot

set.seed(1234)

# load the integrated object
cat("Loading integrated object...\n")
obj <- readRDS(input_object)
cat("Cells:", ncol(obj), "\n")

cluster_levels <- sort_cluster_levels(obj@meta.data[[cluster_col]])
cluster_ids <- as.character(obj@meta.data[[cluster_col]])
cat("Clusters:", length(cluster_levels), "\n")

Idents(obj) <- factor(cluster_ids, levels = cluster_levels)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")

# normalize RNA for marker plots and marker detection
cat("Normalizing RNA assay...\n")
obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)

# find positive markers for each cluster
cat("\nRunning FindAllMarkers...\n")
markers <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = min_pct,
  logfc.threshold = logfc_thr,
  test.use = "wilcox",
  verbose = FALSE
)

cat("Total marker genes found:", nrow(markers), "\n")
write.csv(markers, file.path(res_dir, "all_markers.csv"), row.names = FALSE)

# save the strongest markers per cluster for manual review
top_markers_save <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = n_save, with_ties = FALSE) %>%
  ungroup()
write.csv(top_markers_save, file.path(res_dir, "top_markers_per_cluster.csv"), row.names = FALSE)

top_markers_review <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = n_review_markers, with_ties = FALSE) %>%
  summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")

# dotplot of the strongest marker genes found per cluster
review_features <- top_markers_review$top_markers
review_features <- unique(unlist(strsplit(review_features, ", ", fixed = TRUE)))
review_features <- review_features[review_features %in% rownames(obj)]

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
  ) +
  guides(
    color = guide_colorbar(title = "Average Expression", title.position = "top", barwidth = 5, barheight = 0.5),
    size = guide_legend(title = "Percent Expressed", title.position = "top", nrow = 1, byrow = TRUE)
  )

ggsave(file.path(fig_dir, "dotplot_top_markers.png"), p_marker_dot, width = 12, height = 6, dpi = 600)

# assign final labels from global_config.R
cat("\nAssigning final cell type labels...\n")
obj$cell_type <- unname(cluster_celltype[cluster_ids])

cluster_annotation_table <- data.frame(
  cluster_id = cluster_levels,
  cell_type = unname(cluster_celltype[cluster_levels]),
  n_cells = as.integer(table(factor(cluster_ids, levels = cluster_levels))),
  stringsAsFactors = FALSE
)
write.csv(cluster_annotation_table, file.path(res_dir, "cluster_annotation_table.csv"), row.names = FALSE)

# save the canonical marker table used for the final annotation dotplot
canonical_marker_table <- build_marker_table(marker_genes, rownames(obj))
write.csv(canonical_marker_table, file.path(res_dir, "canonical_marker_table.csv"), row.names = FALSE)

# plot the final cell type UMAP
cat("\nPlotting final cell type UMAP...\n")
celltype_levels <- unique(cluster_annotation_table$cell_type)
final_colors <- build_celltype_colors(celltype_levels, cluster_name_colors)

p_umap <- DimPlot(
  object = obj,
  reduction = reduction_name,
  group.by = "cell_type",
  cols = final_colors,
  raster = FALSE,
  pt.size = pt_size
) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "right",
    legend.key.size = grid::unit(0.45, "cm"),
    legend.spacing.y = grid::unit(0.12, "cm"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

ggsave(file.path(fig_dir, "umap_by_celltype.png"), p_umap, width = 10, height = 6, dpi = 600)

# split the final UMAP by disease group
for (split_col in c("condition", "condition_all")) {
  split_values <- as.character(obj@meta.data[[split_col]])
  split_values[is.na(split_values) | trimws(split_values) == ""] <- "NA"
  obj@meta.data[[split_col]] <- factor(split_values, levels = unique(split_values))

  p_umap_split <- DimPlot(
    object = obj,
    reduction = reduction_name,
    group.by = "cell_type",
    split.by = split_col,
    cols = final_colors,
    raster = FALSE,
    pt.size = pt_size
  ) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 16) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 20, color = "black"),
      plot.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 20, face = "bold", color = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      legend.position = "right",
      legend.key.size = grid::unit(0.45, "cm"),
      legend.spacing.y = grid::unit(0.12, "cm"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

  ggsave(
    file.path(fig_dir, paste0("umap_by_celltype_split_", split_col, ".png")),
    p_umap_split,
    width = 7 * length(levels(obj@meta.data[[split_col]])),
    height = 6,
    dpi = 600
  )
}

# plot cell type ratios by sample and condition
meta_df <- obj@meta.data
ratio_colors <- build_celltype_colors(obj$cell_type, cluster_name_colors)

for (group_col in c("sample_id", "condition", "condition_all")) {
  group_values <- as.character(meta_df[[group_col]])
  group_values[is.na(group_values) | trimws(group_values) == ""] <- "NA"
  meta_df[[group_col]] <- group_values

  ratio_df <- build_ratio_plot_data(meta_df, group_col, "cell_type")
  p_ratio <- make_ratio_plot(ratio_df, ratio_colors)
  plot_width <- ifelse(group_col == "sample_id", max(8, 1.2 * length(unique(ratio_df$group))), max(12, 3.2 * length(unique(ratio_df$group))))

  ggsave(file.path(fig_dir, paste0("celltype_ratio_by_", group_col, ".png")), p_ratio, width = plot_width, height = 10, dpi = 600)
}

# dotplot of canonical markers used for the final cell type labels
cat("\nPlotting canonical marker dotplot...\n")
Idents(obj) <- factor(obj$cell_type, levels = rev(celltype_levels))
dot_features <- unique(as.character(canonical_marker_table$gene))
dot_features <- dot_features[dot_features %in% rownames(obj)]

p_dot <- DotPlot(
  object = obj,
  features = dot_features,
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
  ) +
  guides(
    color = guide_colorbar(title = "Average Expression", title.position = "top", barwidth = 5, barheight = 0.5),
    size = guide_legend(title = "Percent Expressed", title.position = "top", nrow = 1, byrow = TRUE)
  )

ggsave(file.path(fig_dir, "dotplot_canonical_markers_by_celltype.png"), p_dot, width = 12, height = 6, dpi = 600)

# feature plots for selected annotation genes
feature_markers <- unique(trimws(feature_markers))
feature_markers <- feature_markers[feature_markers %in% rownames(obj)]

p_feature_panel <- FeaturePlot(
  object = obj,
  features = feature_markers,
  reduction = reduction_name,
  raster = FALSE,
  order = TRUE,
  ncol = feature_ncol
)

ggsave(
  file.path(fig_dir, "featureplot_selected_markers.png"),
  p_feature_panel,
  width = 4 * min(feature_ncol, length(feature_markers)),
  height = 3.8 * ceiling(length(feature_markers) / feature_ncol),
  dpi = 600
)

# save the object
cat("\nSaving final annotated object...\n")
saveRDS(obj, output_object)

cat("\n============================================================\n")
cat("Cluster marker review and final annotation complete.\n")
cat("Annotated object: ", output_object, "\n", sep = "")
cat("Results dir     : ", res_dir, "\n", sep = "")
cat("Figures dir     : ", fig_dir, "\n", sep = "")
cat("============================================================\n")
