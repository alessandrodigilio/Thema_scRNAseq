###########################################
### Activated endothelial subclustering ###
###########################################

# subset the annotated object to the activated endothelial cells,
# rerun dimensional reduction and clustering on the subset, compute
# subcluster markers, write review tables and save an updated object
# with a second endothelial-specific annotation layer.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(harmony)
})

# work from the project root
setwd("~/Thema_R")
source("src/global_config.R")

# create output directories used by this step
input_object <- file.path(data_dir, "integrated_object", "annotated.rds")

out_dir <- file.path(data_dir, "integrated_object")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(figures_dir, "activated_endothelial_subclustering")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

res_dir <- file.path(results_dir, "activated_endothelial_subclustering")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

output_subset_object <- file.path(out_dir, "activated_endothelial_subclustered.rds")
output_full_object <- file.path(out_dir, "annotated_endothelial_subclusters.rds")
manual_labels_file <- file.path(res_dir, "activated_endothelial_subcluster_manual_labels.csv")

# set parameters
parent_celltype <- "Activated endothelial cells"
group_col <- "condition"
batch_col <- "sample_id"
reduction_name <- "umap.harmony.endothelial"
cluster_col <- "endothelial_subcluster"
n_pcs_endo <- 20
harmony_dims_endo <- 1:20
cluster_res_grid_endo <- seq(0.2, 0.8, by = 0.1)
selected_cluster_res_endo <- 0.3
min_pct <- 0.25
logfc_thr <- 0.25
n_save <- 5
n_review_markers <- 3
n_top_labels_review <- 3

iron_markers <- c(
  "HMOX1", "SLC40A1", "TFRC", "STEAP3", "FTH1", "FTL",
  "SLC11A2", "CP", "NCOA4", "GPX4", "ACSL4", "AIFM2",
  "NQO1", "GCLC", "GCLM", "ALOX5", "ALOX15", "SAT1"
)

set.seed(1234)

sort_cluster_levels <- function(x) {
  x <- unique(as.character(x))
  suppressWarnings(x_num <- as.integer(x))
  if (all(!is.na(x_num))) return(as.character(sort(x_num)))
  sort(x)
}

summarize_label_composition <- function(meta_df, cluster_col, label_col, prefix, n_top = 3) {
  if (!label_col %in% colnames(meta_df)) return(NULL)

  df <- meta_df %>%
    filter(!is.na(.data[[label_col]]) & trimws(.data[[label_col]]) != "") %>%
    group_by(.data[[cluster_col]], .data[[label_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(freq = n / sum(n)) %>%
    arrange(.data[[cluster_col]], desc(freq), desc(n)) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= n_top) %>%
    ungroup()

  if (nrow(df) == 0) return(NULL)

  label_wide <- data.frame(
    cluster_id = as.character(df[[cluster_col]]),
    rank = df$rank,
    value = as.character(df[[label_col]]),
    stringsAsFactors = FALSE
  ) %>%
    pivot_wider(
      names_from = rank,
      values_from = value,
      names_glue = paste0(prefix, "_top{rank}_label")
    )

  freq_wide <- data.frame(
    cluster_id = as.character(df[[cluster_col]]),
    rank = df$rank,
    value = df$freq,
    stringsAsFactors = FALSE
  ) %>%
    pivot_wider(
      names_from = rank,
      values_from = value,
      names_glue = paste0(prefix, "_top{rank}_freq")
    )

  left_join(label_wide, freq_wide, by = "cluster_id")
}

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

# subset to activated endothelial cells only
cells_use <- rownames(obj_full@meta.data)[obj_full$cell_type == parent_celltype]

if (length(cells_use) == 0) {
  stop("No cells found for the selected endothelial population")
}

obj <- subset(obj_full, cells = cells_use)
cat("Activated endothelial subset cells:", ncol(obj), "\n")
print(table(obj$cell_type))

# recompute subset embedding and clustering
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = n_pcs_endo, verbose = FALSE)

if (!batch_col %in% colnames(obj@meta.data)) {
  stop("Batch column not found in metadata: ", batch_col)
}

cat("\nRunning Harmony on activated endothelial subset...\n")
obj <- RunHarmony(
  object = obj,
  group.by.vars = batch_col,
  reduction = "pca",
  dims.use = harmony_dims_endo,
  reduction.save = "harmony",
  verbose = FALSE
)

obj <- RunUMAP(
  object = obj,
  reduction = "harmony",
  dims = harmony_dims_endo,
  reduction.name = reduction_name,
  reduction.key = "UMAPENDO_",
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  verbose = FALSE
)

obj <- FindNeighbors(
  object = obj,
  reduction = "harmony",
  dims = harmony_dims_endo,
  k.param = 20,
  verbose = FALSE
)

cluster_summary <- data.frame(
  resolution = cluster_res_grid_endo,
  n_clusters = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_along(cluster_res_grid_endo)) {
  res_here <- cluster_res_grid_endo[i]
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

write.csv(cluster_summary, file.path(res_dir, "activated_endothelial_cluster_resolution_summary.csv"), row.names = FALSE)

selected_cluster_col <- paste0("RNA_snn_res.", selected_cluster_res_endo)
if (!selected_cluster_col %in% colnames(obj@meta.data)) {
  stop("Selected clustering column not found: ", selected_cluster_col)
}

obj[[cluster_col]] <- as.character(obj@meta.data[[selected_cluster_col]])
cluster_levels <- sort_cluster_levels(obj@meta.data[[cluster_col]])
Idents(obj) <- factor(obj@meta.data[[cluster_col]], levels = cluster_levels)

# marker finding for endothelial subclusters
cat("\nRunning FindAllMarkers on activated endothelial subclusters...\n")
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

write.csv(markers, file.path(res_dir, "all_markers_activated_endothelial_subclusters.csv"), row.names = FALSE)

top_markers_save <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = n_save, with_ties = FALSE) %>%
  ungroup()

write.csv(top_markers_save, file.path(res_dir, "top_markers_per_activated_endothelial_subcluster.csv"), row.names = FALSE)

top_markers_review <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = n_review_markers, with_ties = FALSE) %>%
  summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
top_markers_review$cluster <- as.character(top_markers_review$cluster)

# review table
cluster_counts <- obj@meta.data %>%
  count(.data[[cluster_col]], name = "n_cells")
cluster_counts[[cluster_col]] <- as.character(cluster_counts[[cluster_col]])

cluster_review <- data.frame(cluster_id = cluster_levels, stringsAsFactors = FALSE) %>%
  left_join(cluster_counts, by = c("cluster_id" = cluster_col)) %>%
  left_join(top_markers_review, by = c("cluster_id" = "cluster"))

for (meta_col in c("cell_type", "sample_id", "condition_all", "condition", "sex")) {
  meta_top <- summarize_label_composition(
    meta_df = obj@meta.data,
    cluster_col = cluster_col,
    label_col = meta_col,
    prefix = meta_col,
    n_top = n_top_labels_review
  )

  if (!is.null(meta_top)) {
    cluster_review <- cluster_review %>%
      left_join(meta_top, by = "cluster_id")
  }
}

write.csv(cluster_review, file.path(res_dir, "activated_endothelial_subcluster_review_table.csv"), row.names = FALSE)

if (!file.exists(manual_labels_file)) {
  manual_template <- cluster_review %>%
    transmute(
      cluster_id = cluster_id,
      marker_hint = top_markers,
      manual_label = "",
      notes = ""
    )

  write.csv(manual_template, manual_labels_file, row.names = FALSE)
  cat("\nManual label template created:\n", manual_labels_file, "\n")
} else {
  cat("\nUsing existing manual label template:\n", manual_labels_file, "\n")
}

# save subset object
saveRDS(obj, output_subset_object)

# add subcluster layer back to the full annotated object
obj_full$endothelial_subcluster <- NA_character_
obj_full$endothelial_subcluster[Cells(obj)] <- as.character(obj$endothelial_subcluster)
obj_full$endothelial_parent_type <- NA_character_
obj_full$endothelial_parent_type[Cells(obj)] <- as.character(obj$cell_type)

saveRDS(obj_full, output_full_object)

# review plots
cat("\nPlotting activated endothelial subcluster review panels...\n")

review_plots <- list(
  DimPlot(
    object = obj,
    reduction = reduction_name,
    group.by = cluster_col,
    raster = FALSE,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Activated endothelial subclusters"),
  DimPlot(
    object = obj,
    reduction = reduction_name,
    group.by = "cell_type",
    raster = FALSE,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Parent endothelial type")
)

for (plot_col in c("sample_id", "condition")) {
  if (plot_col %in% colnames(obj@meta.data)) {
    review_plots[[length(review_plots) + 1]] <- DimPlot(
      object = obj,
      reduction = reduction_name,
      group.by = plot_col,
      raster = FALSE,
      label = TRUE,
      repel = TRUE
    ) + ggtitle(plot_col)
  }
}

p_review <- wrap_plots(review_plots, ncol = 2)

ggsave(
  file.path(fig_dir, "activated_endothelial_subclusters_review.png"),
  p_review,
  width = 16,
  height = 12,
  dpi = 600
)

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
    file.path(fig_dir, "dotplot_top_markers_activated_endothelial_subclusters.png"),
    p_marker_dot,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# dotplot focused on iron-related genes
iron_features <- iron_markers[iron_markers %in% rownames(obj)]

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
    file.path(fig_dir, "dotplot_iron_related_markers_activated_endothelial_subclusters.png"),
    p_iron_dot,
    width = 12,
    height = 6,
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Activated endothelial subclustering complete.\n")
cat("Subset object : ", output_subset_object, "\n")
cat("Updated object: ", output_full_object, "\n")
cat("Results dir   : ", res_dir, "\n")
cat("Figures dir   : ", fig_dir, "\n")
cat("============================================================\n")
