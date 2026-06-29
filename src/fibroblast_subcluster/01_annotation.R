######################################################
### Final destructive lining fibroblast annotation ###
######################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/fibroblast_subcluster/utils.R")

# directories
input_subset_object <- file.path(data_dir, "integrated_object", "destructive_lining_fibroblasts_subclustered.rds")
input_full_object <- file.path(data_dir, "integrated_object", "annotated_destructive_lining_fibroblast_subclusters.rds")
out_dir <- file.path(data_dir, "integrated_object")
fig_dir <- file.path(figures_dir, "destructive_lining_fibroblast_final_annotation")
res_dir <- file.path(results_dir, "destructive_lining_fibroblast_final_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
output_full_object <- file.path(out_dir, "annotated_destructive_lining_fibroblast_states.rds")

# set parameters
cluster_col <- "destructive_lining_fibroblast_subcluster"
label_col <- "destructive_lining_fibroblast_subtype"
reduction_name <- "umap.harmony.destructive.lining.fibroblast"
pt_size <- 0.18
min_cells_per_sample_ratio_plot <- 5
low_cell_count_label <- "Low cell count (<5 cells)"

condition_colors <- c(
  "other" = "#B65A5A",
  "HA" = "#5B8DB8"
)

# load subset and updated full object
cat("Loading destructive lining fibroblast subset object...\n")
obj_fibro <- readRDS(input_subset_object)
cat("Destructive lining fibroblast cells:", ncol(obj_fibro), "\n")

cat("Loading full annotated object...\n")
obj_full <- readRDS(input_full_object)
cat("Full object cells:", ncol(obj_full), "\n")

DefaultAssay(obj_fibro) <- "RNA"
obj_fibro <- JoinLayers(obj_fibro, assay = "RNA")
obj_fibro <- NormalizeData(obj_fibro, assay = "RNA", verbose = FALSE)

# apply the final subcluster-to-state annotation
cluster_ids <- as.character(obj_fibro@meta.data[[cluster_col]])
obj_fibro <- add_destructive_lining_fibroblast_subtype_labels(obj_fibro, cluster_col = cluster_col, label_col = label_col)

cat("Destructive lining fibroblast subtype labels after mapping:\n")
print(table(obj_fibro@meta.data[[label_col]], useNA = "ifany"))

obj_full@meta.data[[label_col]] <- NA_character_
obj_full@meta.data[Cells(obj_fibro), label_col] <- as.character(obj_fibro@meta.data[[label_col]])

# save the final annotation table
cluster_levels <- unique(cluster_ids)
cluster_levels <- sort_cluster_levels(cluster_levels)
cluster_counts <- table(cluster_ids)

cluster_annotation_table <- data.frame(
  destructive_lining_fibroblast_subcluster = cluster_levels,
  destructive_lining_fibroblast_subtype = unname(destructive_lining_fibroblast_subcluster_labels[cluster_levels]),
  n_cells = as.integer(cluster_counts[cluster_levels]),
  stringsAsFactors = FALSE
)

write.csv(cluster_annotation_table, file.path(res_dir, "destructive_lining_fibroblast_subcluster_annotation_table.csv"), row.names = FALSE)

canonical_marker_table <- make_destructive_lining_fibroblast_marker_table(obj_fibro)
write.csv(canonical_marker_table, file.path(res_dir, "destructive_lining_fibroblast_canonical_marker_table.csv"), row.names = FALSE)

# save the updated full object with final fibroblast state labels
saveRDS(obj_full, output_full_object)

cat("\nPlotting destructive lining fibroblast UMAPs...\n")
fibro_levels <- levels(obj_fibro@meta.data[[label_col]])
fibro_color_map <- get_destructive_lining_fibroblast_colors(fibro_levels)

p_umap <- DimPlot(
  object = obj_fibro,
  reduction = reduction_name,
  group.by = label_col,
  raster = FALSE,
  pt.size = pt_size
) +
  scale_color_manual(values = fibro_color_map, drop = FALSE) +
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
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.key.size = grid::unit(0.45, "cm"),
    legend.spacing.y = grid::unit(0.12, "cm"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

ggsave(
  filename = file.path(fig_dir, "umap_by_destructive_lining_fibroblast_subtype.png"),
  plot = p_umap,
  width = 12,
  height = 6.5,
  dpi = 600
)

# condition-colored UMAP view by HA vs other condition
condition_values <- ifelse(as.character(obj_fibro@meta.data$condition) == "HA", "HA", "other")
obj_fibro@meta.data$condition <- factor(condition_values, levels = c("HA", "other"))
condition_color_map <- condition_colors[levels(obj_fibro@meta.data$condition)]

p_umap_condition <- DimPlot(
  object = obj_fibro,
  reduction = reduction_name,
  group.by = "condition",
  raster = FALSE,
  pt.size = pt_size
) +
  scale_color_manual(values = condition_color_map, drop = FALSE) +
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
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

ggsave(
  filename = file.path(fig_dir, "umap_by_condition.png"),
  plot = p_umap_condition,
  width = 9,
  height = 6.5,
  dpi = 600
)

# plot split UMAPs by condition and condition_all
for (split_col in c("condition", "condition_all")) {
  split_values <- as.character(obj_fibro@meta.data[[split_col]])
  split_values[is.na(split_values) | trimws(split_values) == ""] <- "NA"
  obj_fibro@meta.data[[split_col]] <- factor(split_values, levels = unique(split_values))
  obj_fibro@meta.data[[split_col]] <- droplevels(obj_fibro@meta.data[[split_col]])

  p_umap_split <- DimPlot(
    object = obj_fibro,
    reduction = reduction_name,
    group.by = label_col,
    split.by = split_col,
    raster = FALSE,
    pt.size = pt_size
  ) +
    scale_color_manual(values = fibro_color_map, drop = FALSE) +
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
      legend.text = element_text(size = 12),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(size = 13, face = "bold", color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

  ggsave(
    filename = file.path(fig_dir, paste0("umap_by_destructive_lining_fibroblast_subtype_split_", split_col, ".png")),
    plot = p_umap_split,
    width = 13,
    height = 6.5,
    dpi = 600
  )
}

# plot the two usual composition barplots by condition and condition_all
for (x_col in c("condition", "condition_all")) {
  plot_df <- build_fibroblast_ratio_plot_data(
    meta_df = obj_fibro@meta.data,
    x_col = x_col,
    fill_col = label_col
  )

  p_ratio <- plot_fibroblast_ratio(plot_df, fibro_color_map)

  ggsave(
    filename = file.path(fig_dir, paste0("destructive_lining_fibroblast_subtype_ratio_by_", x_col, ".png")),
    plot = p_ratio,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# plot composition per sample with HA and other separated in the x-axis lower strip
sample_ratio_df <- build_fibroblast_sample_ratio_plot_data(
  meta_df = obj_fibro@meta.data,
  sample_col = "sample_id",
  condition_col = "condition",
  fill_col = label_col,
  min_cells = min_cells_per_sample_ratio_plot,
  low_cell_label = low_cell_count_label
)

write.csv(sample_ratio_df[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total", "ratio", "is_low_cell_count")], file.path(res_dir, "destructive_lining_fibroblast_subtype_ratio_by_sample.csv"), row.names = FALSE)

p_sample_ratio <- plot_fibroblast_sample_ratio(sample_ratio_df, fibro_color_map, low_cell_label = low_cell_count_label)

ggsave(
  filename = file.path(fig_dir, "destructive_lining_fibroblast_subtype_ratio_by_sample_grouped_condition.png"),
  plot = p_sample_ratio,
  width = 13,
  height = 6.5,
  dpi = 600
)

# plot the canonical marker dotplot for final fibroblast states
marker_features <- unlist(marker_genes_destructive_lining_fibroblast, use.names = FALSE)
marker_features <- unique(marker_features[marker_features %in% rownames(obj_fibro)])

p_dot <- DotPlot(
  object = obj_fibro,
  features = marker_features,
  group.by = label_col,
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
    legend.key.width = grid::unit(0.35, "cm"),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.spacing.x = grid::unit(0.10, "cm"),
    plot.margin = margin(20, 35, 35, 25)
  )

ggsave(
  filename = file.path(fig_dir, "dotplot_canonical_markers_by_destructive_lining_fibroblast_subtype.png"),
  plot = p_dot,
  width = 15,
  height = 6.8,
  dpi = 600
)

# plot iron-related FeaturePlots to inspect iron-handling states
iron_features_use <- destructive_lining_fibroblast_iron_features[
  destructive_lining_fibroblast_iron_features %in% rownames(obj_fibro)
]

p_feature <- FeaturePlot(
  object = obj_fibro,
  features = iron_features_use,
  reduction = reduction_name,
  raster = FALSE,
  order = TRUE,
  cols = c("#F3F1EC", "#7A1F2B"),
  ncol = 4
)

ggsave(
  filename = file.path(fig_dir, "featureplot_iron_related_genes_destructive_lining_fibroblast_subtype.png"),
  plot = p_feature,
  width = 15,
  height = 11,
  dpi = 600
)

# focused HMOX1 plots in the final fibroblast states
hmox1_gene <- "HMOX1"

hmox1_df <- FetchData(obj_fibro, vars = c(hmox1_gene, label_col))
colnames(hmox1_df) <- c("expression", "subtype")

hmox1_summary <- aggregate(
  expression ~ subtype,
  data = hmox1_df,
  FUN = function(x) c(mean = mean(x), median = median(x), pct_expr = mean(x > 0) * 100)
)

hmox1_summary <- data.frame(
  subtype = hmox1_summary$subtype,
  mean_expression = hmox1_summary$expression[, "mean"],
  median_expression = hmox1_summary$expression[, "median"],
  pct_expressing = hmox1_summary$expression[, "pct_expr"],
  stringsAsFactors = FALSE
)

write.csv(hmox1_summary, file.path(res_dir, "HMOX1_expression_summary_destructive_lining_fibroblast_subtypes.csv"), row.names = FALSE)
cat("HMOX1 summary:\n")
print(hmox1_summary)

p_hmox1_feature <- FeaturePlot(
  object = obj_fibro,
  features = hmox1_gene,
  reduction = reduction_name,
  raster = FALSE,
  order = TRUE,
  pt.size = 1.05,
  cols = c("#F3F1EC", "#7A1F2B")
) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 24) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 28, color = "black"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 19, color = "black"),
    plot.margin = margin(28, 28, 28, 28)
  )

ggsave(
  filename = file.path(fig_dir, "featureplot_HMOX1_destructive_lining_fibroblast_subtypes.png"),
  plot = p_hmox1_feature,
  width = 10,
  height = 8.2,
  dpi = 600
)

p_hmox1_dot <- DotPlot(
  object = obj_fibro,
  features = hmox1_gene,
  group.by = label_col,
  cols = c("#E8E2DC", "#7A1F2B"),
  dot.scale = 8,
  col.min = 0,
  col.max = 3
) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.justification = c(0.65, 0.5),
    legend.key.width = grid::unit(0.45, "cm"),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.spacing.x = grid::unit(0.10, "cm"),
    plot.margin = margin(20, 35, 25, 45)
  )

ggsave(
  filename = file.path(fig_dir, "dotplot_HMOX1_destructive_lining_fibroblast_subtypes.png"),
  plot = p_hmox1_dot,
  width = 9,
  height = 4.8,
  dpi = 600
)

cat("\n============================================================\n")
cat("Final destructive lining fibroblast annotation complete.\n")
cat("Updated object: ", output_full_object, "\n")
cat("Results dir   : ", res_dir, "\n")
cat("Figures dir   : ", fig_dir, "\n")
cat("Batch effect  : mitigated in the subset by Harmony on sample_id.\n")
cat("============================================================\n")
