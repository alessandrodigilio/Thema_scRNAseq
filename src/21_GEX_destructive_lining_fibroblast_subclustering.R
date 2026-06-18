# ---------------------------------------
# 21_GEX_destructive_lining_fibroblast_subclustering.R
# Subclustering of MMP3+ lining fibroblasts
# ---------------------------------------

# Subset the annotated object to destructive/MMP3+ lining fibroblasts,
# rerun dimensional reduction and clustering on the subset, compute
# subcluster markers, write review tables and save an updated object
# with a fibroblast-specific subcluster layer.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(harmony)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
INPUT_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_annotated.rds")

OUT_DIR <- file.path(DATA_DIR, "integrated_object")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(FIGURES_DIR, "destructive_lining_fibroblast_subclustering")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RES_DIR <- file.path(RESULTS_DIR, "destructive_lining_fibroblast_subclustering")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_SUBSET_OBJECT <- file.path(OUT_DIR, "gex_destructive_lining_fibroblasts_subclustered.rds")
OUTPUT_FULL_OBJECT <- file.path(OUT_DIR, "gex_annotated_destructive_lining_fibroblast_subclusters.rds")
MANUAL_LABELS_FILE <- file.path(RES_DIR, "destructive_lining_fibroblast_subcluster_manual_labels.csv")

# Set parameters
PARENT_CELLTYPE <- "Destructive lining fibroblasts (MMP3+)"
GROUP_COL <- "condition"
BATCH_COL <- "sample_id"
REDUCTION_NAME <- "umap.harmony.destructive.lining.fibroblast"
CLUSTER_COL <- "destructive_lining_fibroblast_subcluster"
N_PCS_FIBRO <- 20
HARMONY_DIMS_FIBRO <- 1:20
CLUSTER_RES_GRID_FIBRO <- seq(0.2, 0.8, by = 0.1)
SELECTED_CLUSTER_RES_FIBRO <- 0.3
MIN_PCT <- 0.25
LOGFC_THR <- 0.25
N_SAVE <- 5
N_REVIEW_MARKERS <- 3
N_TOP_LABELS_REVIEW <- 3

FIBROBLAST_MARKERS <- c(
  "MMP3", "MMP1", "CXCL1", "CXCL6", "IL6", "PTGS2",
  "PRG4", "CLIC5", "DEFB1", "PDPN", "THY1",
  "COL1A1", "COL1A2", "COL3A1", "COL5A3", "COL6A1",
  "ADAMTS4", "ADAMTS9", "TIMP1", "TIMP2", "HAS1", "VCAM1", "VEGFA"
)

IRON_MARKERS <- c(
  "HMOX1", "SLC40A1", "TFRC", "STEAP3", "FTH1", "FTL",
  "SLC11A2", "CP", "NCOA4", "GPX4", "ACSL4", "AIFM2",
  "NQO1", "GCLC", "GCLM", "ALOX5", "ALOX15", "SAT1",
  "SLC25A37", "SLC25A28", "SLC39A14", "SLC39A8"
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
    dplyr::filter(!is.na(.data[[label_col]]) & trimws(.data[[label_col]]) != "") %>%
    dplyr::group_by(.data[[cluster_col]], .data[[label_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::mutate(freq = n / sum(n)) %>%
    dplyr::arrange(.data[[cluster_col]], dplyr::desc(freq), dplyr::desc(n)) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::filter(rank <= n_top) %>%
    dplyr::ungroup()

  if (nrow(df) == 0) return(NULL)

  label_wide <- data.frame(
    cluster_id = as.character(df[[cluster_col]]),
    rank = df$rank,
    value = as.character(df[[label_col]]),
    stringsAsFactors = FALSE
  ) %>%
    tidyr::pivot_wider(
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
    tidyr::pivot_wider(
      names_from = rank,
      values_from = value,
      names_glue = paste0(prefix, "_top{rank}_freq")
    )

  dplyr::left_join(label_wide, freq_wide, by = "cluster_id")
}

# Load annotated object
if (!file.exists(INPUT_OBJECT)) {
  stop("Missing annotated object: ", INPUT_OBJECT)
}

cat("Loading annotated object...\n")
obj_full <- readRDS(INPUT_OBJECT)
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

# Subset to destructive/MMP3+ lining fibroblasts only
cells_use <- rownames(obj_full@meta.data)[obj_full$cell_type == PARENT_CELLTYPE]

if (length(cells_use) == 0) {
  stop("No cells found for the selected destructive lining fibroblast population")
}

obj <- subset(obj_full, cells = cells_use)
cat("Destructive lining fibroblast subset cells:", ncol(obj), "\n")
print(table(obj$cell_type))

# Recompute subset embedding and clustering
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = N_PCS_FIBRO, verbose = FALSE)

if (!BATCH_COL %in% colnames(obj@meta.data)) {
  stop("Batch column not found in metadata: ", BATCH_COL)
}

cat("\nRunning Harmony on destructive lining fibroblast subset...\n")
obj <- harmony::RunHarmony(
  object = obj,
  group.by.vars = BATCH_COL,
  reduction = "pca",
  dims.use = HARMONY_DIMS_FIBRO,
  reduction.save = "harmony",
  verbose = FALSE
)

obj <- RunUMAP(
  object = obj,
  reduction = "harmony",
  dims = HARMONY_DIMS_FIBRO,
  reduction.name = REDUCTION_NAME,
  reduction.key = "UMAPDLF_",
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  verbose = FALSE
)

obj <- FindNeighbors(
  object = obj,
  reduction = "harmony",
  dims = HARMONY_DIMS_FIBRO,
  k.param = 20,
  verbose = FALSE
)

cluster_summary <- data.frame(
  resolution = CLUSTER_RES_GRID_FIBRO,
  n_clusters = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_along(CLUSTER_RES_GRID_FIBRO)) {
  res_here <- CLUSTER_RES_GRID_FIBRO[i]
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

write.csv(
  cluster_summary,
  file.path(RES_DIR, "destructive_lining_fibroblast_cluster_resolution_summary.csv"),
  row.names = FALSE
)

selected_cluster_col <- paste0("RNA_snn_res.", SELECTED_CLUSTER_RES_FIBRO)
if (!selected_cluster_col %in% colnames(obj@meta.data)) {
  stop("Selected clustering column not found: ", selected_cluster_col)
}

obj[[CLUSTER_COL]] <- as.character(obj@meta.data[[selected_cluster_col]])
cluster_levels <- sort_cluster_levels(obj@meta.data[[CLUSTER_COL]])
Idents(obj) <- factor(obj@meta.data[[CLUSTER_COL]], levels = cluster_levels)

# Marker finding for destructive lining fibroblast subclusters
cat("\nRunning FindAllMarkers on destructive lining fibroblast subclusters...\n")
markers <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = MIN_PCT,
  logfc.threshold = LOGFC_THR,
  test.use = "wilcox",
  verbose = FALSE
)

cat("Total marker genes found:", nrow(markers), "\n")

write.csv(
  markers,
  file.path(RES_DIR, "all_markers_destructive_lining_fibroblast_subclusters.csv"),
  row.names = FALSE
)

top_markers_save <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = N_SAVE, with_ties = FALSE) %>%
  dplyr::ungroup()

write.csv(
  top_markers_save,
  file.path(RES_DIR, "top_markers_per_destructive_lining_fibroblast_subcluster.csv"),
  row.names = FALSE
)

top_markers_review <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = N_REVIEW_MARKERS, with_ties = FALSE) %>%
  dplyr::summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
top_markers_review$cluster <- as.character(top_markers_review$cluster)

# Review table
cluster_counts <- obj@meta.data %>%
  dplyr::count(.data[[CLUSTER_COL]], name = "n_cells")
cluster_counts[[CLUSTER_COL]] <- as.character(cluster_counts[[CLUSTER_COL]])

cluster_review <- data.frame(cluster_id = cluster_levels, stringsAsFactors = FALSE) %>%
  dplyr::left_join(cluster_counts, by = c("cluster_id" = CLUSTER_COL)) %>%
  dplyr::left_join(top_markers_review, by = c("cluster_id" = "cluster"))

for (meta_col in c("cell_type", "sample_id", "condition_all", "condition", "sex")) {
  meta_top <- summarize_label_composition(
    meta_df = obj@meta.data,
    cluster_col = CLUSTER_COL,
    label_col = meta_col,
    prefix = meta_col,
    n_top = N_TOP_LABELS_REVIEW
  )

  if (!is.null(meta_top)) {
    cluster_review <- cluster_review %>%
      dplyr::left_join(meta_top, by = "cluster_id")
  }
}

write.csv(
  cluster_review,
  file.path(RES_DIR, "destructive_lining_fibroblast_subcluster_review_table.csv"),
  row.names = FALSE
)

if (!file.exists(MANUAL_LABELS_FILE)) {
  manual_template <- cluster_review %>%
    dplyr::transmute(
      cluster_id = cluster_id,
      marker_hint = top_markers,
      manual_label = "",
      notes = ""
    )

  write.csv(manual_template, MANUAL_LABELS_FILE, row.names = FALSE)
  cat("\nManual label template created:\n", MANUAL_LABELS_FILE, "\n")
} else {
  cat("\nUsing existing manual label template:\n", MANUAL_LABELS_FILE, "\n")
}

# Save subset object
saveRDS(obj, OUTPUT_SUBSET_OBJECT)

# Add subcluster layer back to the full annotated object
obj_full$destructive_lining_fibroblast_subcluster <- NA_character_
obj_full$destructive_lining_fibroblast_subcluster[Cells(obj)] <- as.character(obj[[CLUSTER_COL]][, 1])
obj_full$destructive_lining_fibroblast_parent_type <- NA_character_
obj_full$destructive_lining_fibroblast_parent_type[Cells(obj)] <- as.character(obj$cell_type)

saveRDS(obj_full, OUTPUT_FULL_OBJECT)

# Review plots
cat("\nPlotting destructive lining fibroblast subcluster review panels...\n")

review_plots <- list(
  DimPlot(
    object = obj,
    reduction = REDUCTION_NAME,
    group.by = CLUSTER_COL,
    raster = FALSE,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("MMP3+ lining fibroblast subclusters"),
  DimPlot(
    object = obj,
    reduction = REDUCTION_NAME,
    group.by = "cell_type",
    raster = FALSE,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Parent fibroblast type")
)

for (plot_col in c("sample_id", "condition")) {
  if (plot_col %in% colnames(obj@meta.data)) {
    review_plots[[length(review_plots) + 1]] <- DimPlot(
      object = obj,
      reduction = REDUCTION_NAME,
      group.by = plot_col,
      raster = FALSE,
      label = TRUE,
      repel = TRUE
    ) + ggtitle(plot_col)
  }
}

p_review <- wrap_plots(review_plots, ncol = 2)

ggsave(
  file.path(FIG_DIR, "destructive_lining_fibroblast_subclusters_review.png"),
  p_review,
  width = 16,
  height = 12,
  dpi = 600
)

# Dotplot of top subcluster markers
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
    file.path(FIG_DIR, "dotplot_top_markers_destructive_lining_fibroblast_subclusters.png"),
    p_marker_dot,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# Dotplot focused on lining/remodeling fibroblast markers
fibro_features <- FIBROBLAST_MARKERS[FIBROBLAST_MARKERS %in% rownames(obj)]

if (length(fibro_features) > 0) {
  p_fibro_dot <- DotPlot(
    object = obj,
    features = fibro_features,
    cols = c("white", "#C65A5A"),
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
    file.path(FIG_DIR, "dotplot_lining_remodeling_markers_destructive_lining_fibroblast_subclusters.png"),
    p_fibro_dot,
    width = 14,
    height = 6,
    dpi = 600
  )
}

# Dotplot focused on iron-related genes
iron_features <- IRON_MARKERS[IRON_MARKERS %in% rownames(obj)]

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
    file.path(FIG_DIR, "dotplot_iron_related_markers_destructive_lining_fibroblast_subclusters.png"),
    p_iron_dot,
    width = 14,
    height = 6,
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Destructive lining fibroblast subclustering complete.\n")
cat("Subset object : ", OUTPUT_SUBSET_OBJECT, "\n")
cat("Updated object: ", OUTPUT_FULL_OBJECT, "\n")
cat("Results dir   : ", RES_DIR, "\n")
cat("Figures dir   : ", FIG_DIR, "\n")
cat("============================================================\n")
