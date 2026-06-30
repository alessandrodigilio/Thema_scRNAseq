###########################################
### Final endothelial subtype annotation ###
###########################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")
source("src/endothelial_subclusters/utils.R")

# directories
input_subset_object <- file.path(data_dir, "integrated_object", "activated_endothelial_subclustered.rds")
input_full_object <- file.path(data_dir, "integrated_object", "annotated_endothelial_subclusters.rds")
out_dir <- file.path(data_dir, "integrated_object")
fig_dir <- file.path(endothelial_figures_dir, "activated_endothelial_final_annotation")
res_dir <- file.path(endothelial_results_dir, "activated_endothelial_final_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

output_full_object <- file.path(out_dir, "annotated_endothelial_states.rds")

# set parameters
cluster_col <- "endothelial_subcluster"
label_col <- "endothelial_subtype"
reduction_name <- "umap.harmony.endothelial"
pt_size <- 0.14

# load endothelial subset and full object with endothelial subcluster IDs
obj_endo <- readRDS(input_subset_object)
obj_full <- readRDS(input_full_object)
cat("Activated endothelial cells:", ncol(obj_endo), "\n")
cat("Full object cells:", ncol(obj_full), "\n")

# normalize RNA expression used by marker plots
DefaultAssay(obj_endo) <- "RNA"
obj_endo <- JoinLayers(obj_endo, assay = "RNA")
obj_endo <- NormalizeData(obj_endo, assay = "RNA", verbose = FALSE)

# apply the final endothelial subtype labels
obj_endo <- add_endothelial_subtype_labels(
  obj = obj_endo,
  cluster_col = cluster_col,
  label_col = label_col
)

cat("Endothelial subtype labels after mapping:\n")
print(table(obj_endo@meta.data[[label_col]], useNA = "ifany"))

# copy endothelial subtype labels back into the full atlas object
obj_full@meta.data[[label_col]] <- NA_character_
obj_full@meta.data[Cells(obj_endo), label_col] <- as.character(obj_endo@meta.data[[label_col]])

# save annotation table
cluster_ids <- as.character(obj_endo@meta.data[[cluster_col]])
cluster_levels <- sort_cluster_levels(unique(cluster_ids))
cluster_counts <- table(cluster_ids)

cluster_annotation_table <- data.frame(
  endothelial_subcluster = cluster_levels,
  endothelial_subtype = unname(endothelial_subcluster_labels[cluster_levels]),
  n_cells = as.integer(cluster_counts[cluster_levels]),
  stringsAsFactors = FALSE
)
#save
write.csv(cluster_annotation_table, file.path(res_dir, "activated_endothelial_subcluster_annotation_table.csv"), row.names = FALSE)

# save canonical marker table
canonical_marker_table <- make_endothelial_marker_table(obj_endo)
write.csv(canonical_marker_table, file.path(res_dir, "activated_endothelial_canonical_marker_table.csv"), row.names = FALSE)

# save the updated full object with endothelial subtype labels
saveRDS(obj_full, output_full_object)

# set colors
endo_levels <- levels(obj_endo@meta.data[[label_col]])
endo_color_map <- get_endothelial_colors(endo_levels)

# plot endothelial UMAP
p_umap <- DimPlot(
  object = obj_endo,
  reduction = reduction_name,
  group.by = label_col,
  raster = FALSE,
  pt.size = pt_size
) +
  scale_color_manual(values = endo_color_map, drop = FALSE) +
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
ggsave(file.path(fig_dir, "umap_by_endothelial_subtype.png"), p_umap, width = 11, height = 6.5, dpi = 600)

# plot endothelial UMAP split by condition labels
for (split_col in c("condition", "condition_all")) {
  split_values <- as.character(obj_endo@meta.data[[split_col]])
  split_values[is.na(split_values) | trimws(split_values) == ""] <- "NA"
  obj_endo@meta.data[[split_col]] <- factor(split_values, levels = unique(split_values))

  p_umap_split <- DimPlot(
    object = obj_endo,
    reduction = reduction_name,
    group.by = label_col,
    split.by = split_col,
    raster = FALSE,
    pt.size = pt_size
  ) +
    scale_color_manual(values = endo_color_map, drop = FALSE) +
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

  ggsave(file.path(fig_dir, paste0("umap_by_endothelial_subtype_split_", split_col, ".png")), p_umap_split, width = 13, height = 6.5, dpi = 600)
}

# plot endothelial subtype ratios by condition labels
for (x_col in c("condition", "condition_all")) {
  plot_df <- build_ratio_plot_data(
    meta_df = obj_endo@meta.data,
    x_col = x_col,
    fill_col = label_col
  )

  p_ratio <- make_ratio_plot(plot_df, unname(endo_color_map))

  ggsave(file.path(fig_dir, paste0("endothelial_subtype_ratio_by_", x_col, ".png")), p_ratio, width = 12, height = 6, dpi = 600)
}

# plot canonical markers used to annotate endothelial subtypes
marker_features <- unlist(marker_genes_endothelial, use.names = FALSE)
marker_features <- unique(marker_features[marker_features %in% rownames(obj_endo)])

dotplot_df <- build_endothelial_marker_dotplot_data(
  object = obj_endo,
  features = marker_features,
  group_col = label_col
)

p_dot <- plot_endothelial_marker_dotplot(
  dotplot_df = dotplot_df,
  low_col = "#F3F1EC",
  mid_col = "#F3F1EC",
  high_col = "#d84b1c"
)
ggsave(file.path(fig_dir, "dotplot_canonical_markers_by_endothelial_subtype.png"), p_dot, width = 15, height = 6.5, dpi = 600)

# plot iron-related markers to inspect iron-handling states
iron_features_use <- endothelial_iron_features[endothelial_iron_features %in% rownames(obj_endo)]
p_feature_iron <- FeaturePlot(
  object = obj_endo,
  features = iron_features_use,
  reduction = reduction_name,
  raster = FALSE,
  order = TRUE,
  cols = c("#F3F1EC", "#7A1F2B"),
  ncol = 4
)
ggsave(file.path(fig_dir, "featureplot_iron_related_genes_endothelial_subtype.png"), p_feature_iron, width = 15, height = 10.5, dpi = 600)

cat("\n============================================================\n")
cat("Final endothelial annotation complete.\n")
cat("Updated object: ", output_full_object, "\n", sep = "")
cat("Results dir   : ", res_dir, "\n", sep = "")
cat("Figures dir   : ", fig_dir, "\n", sep = "")
cat("Batch effect  : mitigated in the subset by Harmony on sample_id.\n")
cat("============================================================\n")
