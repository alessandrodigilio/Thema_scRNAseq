###########################################
### Final macrophage subtype annotation ###
###########################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")
source("src/macrophage_subclusters/utils.R")

# directories
input_subset_object <- file.path(data_dir, "integrated_object", "macrophages_subclustered.rds")
input_full_object <- file.path(data_dir, "integrated_object", "annotated_macrophage_subclusters.rds")
out_dir <- file.path(data_dir, "integrated_object")
fig_dir <- file.path(figures_dir, "macrophage_final_annotation")
res_dir <- file.path(results_dir, "macrophage_final_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

output_full_object <- file.path(out_dir, "annotated_macrophage_states.rds")

# set parameters
cluster_col <- "macrophage_subcluster"
label_col <- "macrophage_subtype"
reduction_name <- "umap.harmony.macrophage"
pt_size <- 0.14

# load macrophage subset and full object with macrophage subcluster IDs
obj_macro <- readRDS(input_subset_object)
obj_full <- readRDS(input_full_object)
cat("Macrophage cells:", ncol(obj_macro), "\n")
cat("Full object cells:", ncol(obj_full), "\n")

# normalize RNA expression used by marker plots
DefaultAssay(obj_macro) <- "RNA"
obj_macro <- JoinLayers(obj_macro, assay = "RNA")
obj_macro <- NormalizeData(obj_macro, assay = "RNA", verbose = FALSE)

# apply the final macrophage subtype labels
obj_macro <- add_macrophage_subtype_labels(
  obj = obj_macro,
  cluster_col = cluster_col,
  label_col = label_col
)

cat("Macrophage subtype labels after mapping:\n")
print(table(obj_macro@meta.data[[label_col]], useNA = "ifany"))

# copy macrophage subtype labels back into the full atlas object
obj_full@meta.data[[label_col]] <- NA_character_
obj_full@meta.data[Cells(obj_macro), label_col] <- as.character(obj_macro@meta.data[[label_col]])

# save the final subcluster annotation table
cluster_ids <- as.character(obj_macro@meta.data[[cluster_col]])
cluster_levels <- sort_cluster_levels(unique(cluster_ids))
cluster_counts <- table(cluster_ids)

# df with cluster IDs, subtype labels, and cell counts
cluster_annotation_table <- data.frame(
  macrophage_subcluster = cluster_levels,
  macrophage_subtype = unname(macrophage_subcluster_labels[cluster_levels]),
  n_cells = as.integer(cluster_counts[cluster_levels]),
  stringsAsFactors = FALSE
)
# save
write.csv(cluster_annotation_table, file.path(res_dir, "macrophage_subcluster_annotation_table.csv"), row.names = FALSE)

# save the canonical marker table (biologically meaningful)
canonical_marker_table <- make_macrophage_marker_table(obj_macro)
write.csv(canonical_marker_table, file.path(res_dir, "macrophage_canonical_marker_table.csv"), row.names = FALSE)

# save the updated full object with macrophage subtype labels
saveRDS(obj_full, output_full_object)

# set colors
macro_levels <- levels(obj_macro@meta.data[[label_col]])
macro_color_map <- get_macrophage_colors(macro_levels)

# plot macrophage UMAP
p_umap <- DimPlot(
  object = obj_macro,
  reduction = reduction_name,
  group.by = label_col,
  raster = FALSE,
  pt.size = pt_size
) +
  scale_color_manual(values = macro_color_map, drop = FALSE) +
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

ggsave(
  filename = file.path(fig_dir, "umap_by_macrophage_subtype.png"),
  plot = p_umap,
  width = 11,
  height = 6.5,
  dpi = 600
)

# plot macrophage UMAP split by condition labels
for (split_col in c("condition", "condition_all")) {
  split_values <- as.character(obj_macro@meta.data[[split_col]])
  split_values[is.na(split_values) | trimws(split_values) == ""] <- "NA"
  obj_macro@meta.data[[split_col]] <- factor(split_values, levels = unique(split_values))

  p_umap_split <- DimPlot(
    object = obj_macro,
    reduction = reduction_name,
    group.by = label_col,
    split.by = split_col,
    raster = FALSE,
    pt.size = pt_size
  ) +
    scale_color_manual(values = macro_color_map, drop = FALSE) +
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
      strip.background = element_blank(),
      strip.text = element_text(size = 13, face = "bold", color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

  ggsave(
    filename = file.path(fig_dir, paste0("umap_by_macrophage_subtype_split_", split_col, ".png")),
    plot = p_umap_split,
    width = 13,
    height = 6.5,
    dpi = 600
  )
}

# plot macrophage subtype ratios by condition labels
for (x_col in c("condition", "condition_all")) {
  plot_df <- build_ratio_plot_data(
    meta_df = obj_macro@meta.data,
    x_col = x_col,
    fill_col = label_col
  )

  p_ratio <- make_ratio_plot(plot_df, unname(macro_color_map))

  ggsave(
    filename = file.path(fig_dir, paste0("macrophage_subtype_ratio_by_", x_col, ".png")),
    plot = p_ratio,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# plot canonical markers used to annotate macrophage subtypes
marker_features <- unlist(marker_genes_macrophage, use.names = FALSE)
marker_features <- unique(marker_features[marker_features %in% rownames(obj_macro)])

dotplot_df <- build_macrophage_marker_dotplot_data(
  object = obj_macro,
  features = marker_features,
  group_col = label_col
)

p_dot <- plot_macrophage_marker_dotplot(
  dotplot_df = dotplot_df,
  low_col = "#F3F1EC",
  mid_col = "#F3F1EC",
  high_col = "#d84b1c"
)

ggsave(
  filename = file.path(fig_dir, "dotplot_canonical_markers_by_macrophage_subtype.png"),
  plot = p_dot,
  width = 15,
  height = 6.5,
  dpi = 600
)

# plot iron-related markers to inspect iron-handling states
iron_features_use <- macrophage_iron_features[macrophage_iron_features %in% rownames(obj_macro)]

p_feature_iron <- FeaturePlot(
  object = obj_macro,
  features = iron_features_use,
  reduction = reduction_name,
  raster = FALSE,
  order = TRUE,
  cols = c("#F3F1EC", "#7A1F2B"),
  ncol = 4
)

ggsave(
  filename = file.path(fig_dir, "featureplot_iron_related_genes_macrophage_subtype.png"),
  plot = p_feature_iron,
  width = 15,
  height = 11,
  dpi = 600
)

# ----------------------------------------------------------------------------------------------------------- #
# --- Among the macrophage subtypes, the red-pulp-like / iron-recycling states are of particular interest --- #
# Check of the markers to see if their profile is consistent with the expected biology of these states ------ #
# ----------------------------------------------------------------------------------------------------------- #

# plot red-pulp-like / iron-recycling markers
red_pulp_features_use <- red_pulp_like_features[red_pulp_like_features %in% rownames(obj_macro)]

# feature plot of red-pulp-like markers
p_feature_red_pulp <- FeaturePlot(
  object = obj_macro,
  features = red_pulp_features_use,
  reduction = reduction_name,
  raster = FALSE,
  order = TRUE,
  cols = c("#F4EDF8", "#6F2DBD"),
  ncol = 4
)
ggsave(
  filename = file.path(fig_dir, "featureplot_red_pulp_like_markers_macrophage_subtype.png"),
  plot = p_feature_red_pulp,
  width = 15,
  height = 10.5,
  dpi = 600
)

# bubble plot
red_pulp_dotplot_df <- build_macrophage_marker_dotplot_data(
  object = obj_macro,
  features = red_pulp_features_use,
  group_col = label_col
)
p_dot_red_pulp <- plot_macrophage_marker_dotplot(
  dotplot_df = red_pulp_dotplot_df,
  low_col = "#F4EDF8",
  mid_col = "#E4D4F1",
  high_col = "#6F2DBD"
)
ggsave(
  filename = file.path(fig_dir, "dotplot_red_pulp_like_markers_by_macrophage_subtype.png"),
  plot = p_dot_red_pulp,
  width = 12,
  height = 6.5,
  dpi = 600
)

cat("\n============================================================\n")
cat("Final macrophage annotation complete.\n")
cat("Updated object: ", output_full_object, "\n")
cat("Results dir   : ", res_dir, "\n")
cat("Figures dir   : ", fig_dir, "\n")
cat("Batch effect  : mitigated in the subset by Harmony on sample_id.\n")
cat("============================================================\n")
