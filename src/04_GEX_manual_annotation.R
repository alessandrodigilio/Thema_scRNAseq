# ---------------------------------------
# 04_GEX_manual_annotation.R
# Manual annotation of integrated GEX
# ---------------------------------------

# Load the integrated object, compute cluster markers, write a cluster
# review table and a manual label template for downstream annotation.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
INPUT_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_integrated.rds")

FIG_DIR <- file.path(FIGURES_DIR, "cluster_naming")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RES_DIR <- file.path(RESULTS_DIR, "cluster_naming")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

MANUAL_LABELS_FILE <- file.path(RES_DIR, "cluster_manual_labels.csv")

# Set parameters
CLUSTER_COL <- "gex_cluster"
REDUCTION_NAME <- "umap.harmony.rna"
MIN_PCT <- 0.25
LOGFC_THR <- 0.25
N_SAVE <- 5
N_REVIEW_MARKERS <- 3
N_TOP_LABELS_REVIEW <- 3
MARKERS_TO_INSPECT <- character(0)

set.seed(1234)

# Sort cluster identifiers numerically when possible
sort_cluster_levels <- function(x) {
  x <- unique(as.character(x))
  suppressWarnings(x_num <- as.integer(x))
  if (all(!is.na(x_num))) return(as.character(sort(x_num)))
  sort(x)
}

# Summarize the top labels per cluster in a wide review table.
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

# Load the integrated object
if (!file.exists(INPUT_OBJECT)) {
  stop("Missing integrated object: ", INPUT_OBJECT)
}

cat("Loading integrated object...\n")
obj <- readRDS(INPUT_OBJECT)
cat("Cells:", ncol(obj), "\n")

if (!CLUSTER_COL %in% colnames(obj@meta.data)) {
  stop("Cluster column not found in metadata: ", CLUSTER_COL)
}

cluster_levels <- sort_cluster_levels(obj@meta.data[[CLUSTER_COL]])
cat("Clusters:", length(cluster_levels), "\n")

Idents(obj) <- factor(obj@meta.data[[CLUSTER_COL]], levels = cluster_levels)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")

# Normalize the RNA assay only if the data layer is still empty
rna_data_layer <- tryCatch(
  LayerData(obj[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  cat("Normalizing RNA assay...\n")
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

# Run FindAllMarkers once to support manual review
cat("\nRunning FindAllMarkers...\n")
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
  file.path(RES_DIR, "all_markers.csv"),
  row.names = FALSE
)

top_markers_save <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = N_SAVE, with_ties = FALSE) %>%
  dplyr::ungroup()

write.csv(
  top_markers_save,
  file.path(RES_DIR, "top_markers_per_cluster.csv"),
  row.names = FALSE
)

top_markers_review <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = N_REVIEW_MARKERS, with_ties = FALSE) %>%
  dplyr::summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop")
top_markers_review$cluster <- as.character(top_markers_review$cluster)

# Summarize metadata information per cluster for manual review
cluster_counts <- obj@meta.data %>%
  dplyr::count(.data[[CLUSTER_COL]], name = "n_cells")
cluster_counts[[CLUSTER_COL]] <- as.character(cluster_counts[[CLUSTER_COL]])

cluster_review <- data.frame(cluster_id = cluster_levels, stringsAsFactors = FALSE) %>%
  dplyr::left_join(cluster_counts, by = c("cluster_id" = CLUSTER_COL)) %>%
  dplyr::left_join(top_markers_review, by = c("cluster_id" = "cluster"))

for (meta_col in c("sample_id", "condition_all", "condition", "sex")) {
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
  file.path(RES_DIR, "cluster_review_table.csv"),
  row.names = FALSE
)

# Write a template for manual labels if it does not exist yet
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

# Plot a compact review panel with cluster IDs and sample IDs
if (REDUCTION_NAME %in% names(obj@reductions)) {
  review_plots <- list(
    DimPlot(
      object = obj,
      reduction = REDUCTION_NAME,
      group.by = CLUSTER_COL,
      raster = FALSE,
      label = TRUE,
      repel = TRUE
    ) + ggtitle("RNA clusters")
  )

  if ("sample_id" %in% colnames(obj@meta.data)) {
    review_plots[[length(review_plots) + 1]] <- DimPlot(
      object = obj,
      reduction = REDUCTION_NAME,
      group.by = "sample_id",
      raster = FALSE,
      label = TRUE,
      repel = TRUE
    ) + ggtitle("Samples")
  }

  p_review <- wrap_plots(review_plots, ncol = length(review_plots))

  ggsave(
    file.path(FIG_DIR, "RNA_clusters_review.png"),
    p_review,
    width = 8 * length(review_plots),
    height = 8,
    dpi = 600
  )
}

# Dotplot of the top marker genes found per cluster
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
    ) +
    guides(
      color = guide_colorbar(
        title = "Average Expression",
        title.position = "top",
        barwidth = 5,
        barheight = 0.5
      ),
      size = guide_legend(
        title = "Percent Expressed",
        title.position = "top",
        nrow = 1,
        byrow = TRUE
      )
    )

  ggsave(
    file.path(FIG_DIR, "dotplot_top_markers.png"),
    p_marker_dot,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# Plot selected marker UMAPs only when markers were explicitly requested
marker_features <- unique(trimws(MARKERS_TO_INSPECT))
marker_features <- marker_features[marker_features != ""]
marker_features <- marker_features[marker_features %in% rownames(obj)]

if (length(marker_features) > 0 && REDUCTION_NAME %in% names(obj@reductions)) {
  cat("\nPlotting marker UMAPs for manual review...\n")

  for (feature_name in marker_features) {
    p_marker_umap <- FeaturePlot(
      object = obj,
      features = feature_name,
      reduction = REDUCTION_NAME,
      raster = FALSE,
      order = TRUE
    ) +
      ggtitle(paste0("UMAP - ", feature_name))

    ggsave(
      filename = file.path(FIG_DIR, paste0("marker_umap_", feature_name, ".png")),
      plot = p_marker_umap,
      width = 8,
      height = 7,
      dpi = 600
    )
  }
}

cat("\n============================================================\n")
cat("Cluster review preparation complete.\n")
cat("Files ready for manual annotation:\n")
cat(" - ", file.path(RES_DIR, "cluster_review_table.csv"), "\n", sep = "")
cat(" - ", file.path(RES_DIR, "top_markers_per_cluster.csv"), "\n", sep = "")
cat(" - ", MANUAL_LABELS_FILE, "\n", sep = "")
cat("============================================================\n")
