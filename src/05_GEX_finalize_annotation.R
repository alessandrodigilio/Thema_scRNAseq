# ---------------------------------------
# 05_GEX_finalize_annotation.R
# Final GEX annotation and plotting
# ---------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
INPUT_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_integrated.rds")

OUT_DIR <- file.path(DATA_DIR, "integrated_object")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(FIGURES_DIR, "final_annotation")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RES_DIR <- file.path(RESULTS_DIR, "final_annotation")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_OBJECT <- file.path(OUT_DIR, "gex_annotated.rds")

# Set parameters
CLUSTER_COL <- "gex_cluster"
REDUCTION_NAME <- "umap.harmony.rna"
PT_SIZE <- 0.01
FEATURE_MARKERS <- c(
  "PRG4", "CLIC5", "HBEGF",
  "MFAP5", "CXCL12", "SFRP4",
  "MMP3", "MMP1", "INHBA",
  "IL1B", "C1QC", "FOLR2",
  "SELE", "ACKR1", "EMCN",
  "RGS5", "ACTA2", "MYH11",
  "CD1C", "FCER1A", "CLEC10A",
  "HMOX1", "SLC40A1"
)
FEATURE_NCOL <- 4

# Build a color vector for the labels present in the object.
build_celltype_colors <- function(labels) {
  labels <- unique(as.character(labels))
  cols <- cluster_name_colors[labels]
  missing_labels <- labels[is.na(cols)]

  if (length(missing_labels) > 0) {
    fallback <- grDevices::hcl.colors(length(missing_labels), palette = "Set 3")
    names(fallback) <- missing_labels
    cols[names(fallback)] <- fallback
  }

  cols[labels]
}

build_ratio_plot_data <- function(meta_df, x_col, fill_col) {
  valid_rows <- !is.na(meta_df[[x_col]]) & trimws(as.character(meta_df[[x_col]])) != ""
  valid_rows <- valid_rows & !is.na(meta_df[[fill_col]]) & trimws(as.character(meta_df[[fill_col]])) != ""

  plot_df <- meta_df[valid_rows, c(x_col, fill_col), drop = FALSE]
  if (nrow(plot_df) == 0) return(NULL)

  plot_df[[x_col]] <- as.character(plot_df[[x_col]])
  plot_df[[fill_col]] <- as.character(plot_df[[fill_col]])

  count_df <- as.data.frame(table(plot_df[[x_col]], plot_df[[fill_col]]), stringsAsFactors = FALSE)
  colnames(count_df) <- c("group", "cell_type", "n_cells")
  count_df <- count_df[count_df$n_cells > 0, , drop = FALSE]
  if (nrow(count_df) == 0) return(NULL)

  totals <- aggregate(n_cells ~ group, data = count_df, FUN = sum)
  count_df <- merge(count_df, totals, by = "group", suffixes = c("", "_total"), sort = FALSE)
  count_df$ratio <- count_df$n_cells / count_df$n_cells_total
  count_df$group <- factor(count_df$group, levels = unique(count_df$group))
  fill_levels <- unique(as.character(meta_df[[fill_col]]))
  count_df$cell_type <- factor(count_df$cell_type, levels = fill_levels)
  count_df
}

make_ratio_plot <- function(plot_df, fill_colors) {
  ggplot(plot_df, aes(x = group, y = ratio, fill = cell_type)) +
    geom_col(width = 0.92, color = NA) +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Ratio", fill = NULL) +
    theme_classic(base_size = 24) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 24, color = "black"),
      axis.text.y = element_text(size = 24, color = "black"),
      axis.title.y = element_text(size = 24, color = "black"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length = grid::unit(0.22, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.position = "right",
      legend.key.size = grid::unit(0.45, "cm"),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Load annotated input object
if (!file.exists(INPUT_OBJECT)) {
  stop("Missing annotated object: ", INPUT_OBJECT)
}

if (!exists("GEX_CLUSTER_CELLTYPE")) {
  stop("GEX_CLUSTER_CELLTYPE not found in global_config.R")
}

if (!exists("cluster_name_colors")) {
  stop("cluster_name_colors not found in global_config.R")
}

if (!exists("marker_genes")) {
  stop("marker_genes not found in global_config.R")
}

cat("Loading annotated object...\n")
obj <- readRDS(INPUT_OBJECT)
cat("Cells:", ncol(obj), "\n")

if (!CLUSTER_COL %in% colnames(obj@meta.data)) {
  stop("Cluster column not found in metadata: ", CLUSTER_COL)
}

if (!REDUCTION_NAME %in% names(obj@reductions)) {
  stop("Reduction not found in object: ", REDUCTION_NAME)
}

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")

# Normalize the RNA assay if the data layer is still missing or empty
rna_data_layer <- tryCatch(
  LayerData(obj[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  cat("Normalizing RNA assay...\n")
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

# Apply final cluster-to-cell-type annotation from global config
cluster_ids <- as.character(obj@meta.data[[CLUSTER_COL]])
obj$cell_type <- unname(GEX_CLUSTER_CELLTYPE[cluster_ids])
obj$cell_type[is.na(obj$cell_type)] <- "Unknown"

# Save a simple cluster annotation table
cluster_levels <- unique(cluster_ids)
suppressWarnings(cluster_levels_num <- as.integer(cluster_levels))
if (all(!is.na(cluster_levels_num))) {
  cluster_levels <- as.character(sort(cluster_levels_num))
} else {
  cluster_levels <- sort(cluster_levels)
}

cluster_annotation_table <- data.frame(
  cluster_id = cluster_levels,
  cell_type = unname(GEX_CLUSTER_CELLTYPE[cluster_levels]),
  n_cells = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(cluster_annotation_table))) {
  cluster_annotation_table$n_cells[i] <- sum(cluster_ids == cluster_annotation_table$cluster_id[i])
}

write.csv(
  cluster_annotation_table,
  file.path(RES_DIR, "cluster_annotation_table.csv"),
  row.names = FALSE
)

# Save a table of canonical markers used for the final dotplot
marker_rows <- list()
marker_index <- 1
marker_names <- names(marker_genes)

for (ct in marker_names) {
  genes_here <- unique(marker_genes[[ct]])
  genes_here <- genes_here[genes_here %in% rownames(obj)]
  if (length(genes_here) == 0) next

  marker_rows[[marker_index]] <- data.frame(
    cell_type = rep(ct, length(genes_here)),
    gene = genes_here,
    stringsAsFactors = FALSE
  )
  marker_index <- marker_index + 1
}

if (length(marker_rows) == 0) {
  stop("No marker genes from global_config.R were found in the object")
}

canonical_marker_table <- do.call(rbind, marker_rows)
write.csv(
  canonical_marker_table,
  file.path(RES_DIR, "canonical_marker_table.csv"),
  row.names = FALSE
)

# Plot and save UMAP colored by final cell type
cat("\nPlotting final UMAP...\n")
final_colors <- build_celltype_colors(obj$cell_type)

p_umap <- DimPlot(
  object = obj,
  reduction = REDUCTION_NAME,
  group.by = "cell_type",
  cols = final_colors,
  raster = FALSE,
  pt.size = PT_SIZE
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

print(p_umap)

ggsave(
  filename = file.path(FIG_DIR, "umap_by_celltype.png"),
  plot = p_umap,
  width = 10,
  height = 6,
  dpi = 600
)

for (split_col in c("condition", "condition_all")) {
  if (split_col %in% colnames(obj@meta.data)) {
    split_values <- as.character(obj@meta.data[[split_col]])
    split_values[is.na(split_values) | trimws(split_values) == ""] <- "NA"
    obj@meta.data[[split_col]] <- factor(split_values, levels = unique(split_values))

    cat("\nPlotting", split_col, "-specific UMAP panel...\n")

    p_umap_split <- DimPlot(
      object = obj,
      reduction = REDUCTION_NAME,
      group.by = "cell_type",
      split.by = split_col,
      cols = final_colors,
      raster = FALSE,
      pt.size = PT_SIZE
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

    print(p_umap_split)

    ggsave(
      filename = file.path(FIG_DIR, paste0("umap_by_celltype_split_", split_col, ".png")),
      plot = p_umap_split,
      width = 7 * length(levels(obj@meta.data[[split_col]])),
      height = 6,
      dpi = 600
    )
  } else {
    warning(split_col, " not found in metadata; skipping split UMAP")
  }
}

cat("\nPlotting cell type ratio barplots...\n")
ratio_colors <- build_celltype_colors(obj$cell_type)
meta_df <- obj@meta.data

for (group_col in c("sample_id", "condition", "condition_all")) {
  if (!group_col %in% colnames(meta_df)) {
    warning(group_col, " not found in metadata; skipping ratio plot")
    next
  }

  group_vals <- as.character(meta_df[[group_col]])
  group_vals[is.na(group_vals) | trimws(group_vals) == ""] <- "NA"
  meta_df[[group_col]] <- group_vals

  ratio_df <- build_ratio_plot_data(meta_df, group_col, "cell_type")
  if (is.null(ratio_df)) next

  p_ratio <- make_ratio_plot(ratio_df, ratio_colors)
  print(p_ratio)

  plot_width <- if (group_col == "sample_id") {
    max(8, 1.2 * length(unique(ratio_df$group)))
  } else {
    max(12, 3.2 * length(unique(ratio_df$group)))
  }

  ggsave(
    filename = file.path(FIG_DIR, paste0("celltype_ratio_by_", group_col, ".png")),
    plot = p_ratio,
    width = plot_width,
    height = 10,
    dpi = 600
  )
}

# Plot and save dotplot using the canonical markers from global config
cat("\nPlotting canonical marker dotplot...\n")
celltype_order <- unique(cluster_annotation_table$cell_type)
celltype_order <- celltype_order[!is.na(celltype_order)]
celltype_order <- celltype_order[celltype_order %in% names(marker_genes)]
Idents(obj) <- factor(obj$cell_type, levels = rev(celltype_order))

dot_features <- unique(as.character(canonical_marker_table$gene))
dot_features <- dot_features[dot_features %in% rownames(obj)]

if (length(dot_features) > 0) {
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

  print(p_dot)

  ggsave(
    filename = file.path(FIG_DIR, "dotplot_canonical_markers_by_celltype.png"),
    plot = p_dot,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# Plot and save feature plots for selected marker genes
FEATURE_MARKERS <- unique(trimws(FEATURE_MARKERS))
FEATURE_MARKERS <- FEATURE_MARKERS[FEATURE_MARKERS != ""]
FEATURE_MARKERS <- FEATURE_MARKERS[FEATURE_MARKERS %in% rownames(obj)]

if (length(FEATURE_MARKERS) > 0) {
  cat("\nPlotting feature plots for selected markers...\n")

  p_feature_panel <- FeaturePlot(
    object = obj,
    features = FEATURE_MARKERS,
    reduction = REDUCTION_NAME,
    raster = FALSE,
    order = TRUE,
    ncol = FEATURE_NCOL
  )

  print(p_feature_panel)

  ggsave(
    filename = file.path(FIG_DIR, "featureplot_selected_markers.png"),
    plot = p_feature_panel,
    width = 4 * min(FEATURE_NCOL, length(FEATURE_MARKERS)),
    height = 3.8 * ceiling(length(FEATURE_MARKERS) / FEATURE_NCOL),
    dpi = 600
  )
}

# Save final annotated object
cat("\nSaving final annotated object...\n")
saveRDS(obj, OUTPUT_OBJECT)

cat("\n============================================================\n")
cat("Final annotation plotting complete.\n")
cat("Object saved : ", OUTPUT_OBJECT, "\n")
cat("Figures dir  : ", FIG_DIR, "\n")
cat("Results dir  : ", RES_DIR, "\n")
cat("============================================================\n")
