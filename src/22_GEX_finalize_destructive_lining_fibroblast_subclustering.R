# ---------------------------------------
# 22_GEX_finalize_destructive_lining_fibroblast_subclustering.R
# Final MMP3+ lining fibroblast subcluster annotation and plotting
# ---------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
INPUT_SUBSET_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_destructive_lining_fibroblasts_subclustered.rds")
INPUT_FULL_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_annotated_destructive_lining_fibroblast_subclusters.rds")

OUT_DIR <- file.path(DATA_DIR, "integrated_object")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(FIGURES_DIR, "destructive_lining_fibroblast_final_annotation")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RES_DIR <- file.path(RESULTS_DIR, "destructive_lining_fibroblast_final_annotation")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_FULL_OBJECT <- file.path(OUT_DIR, "gex_annotated_destructive_lining_fibroblast_states.rds")

# Set parameters
CLUSTER_COL <- "destructive_lining_fibroblast_subcluster"
LABEL_COL <- "destructive_lining_fibroblast_subtype"
REDUCTION_NAME <- "umap.harmony.destructive.lining.fibroblast"
PT_SIZE <- 0.18
MIN_CELLS_PER_SAMPLE_RATIO_PLOT <- 5
LOW_CELL_COUNT_LABEL <- "Low cell count (<5 cells)"

DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS <- c(
  "0" = "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)",
  "1" = "Activated MMP3+ lining fibroblast cells (ID1+)",
  "2" = "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)",
  "3" = "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)",
  "4" = "MMP3+ lining fibroblast cells (FAM184A+)",
  "5" = "HA-enriched SFRP2+ matrix fibroblast-like cells"
)

destructive_lining_fibroblast_subtype_colors <- c(
  "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)" = "#5F7EA6",
  "Activated MMP3+ lining fibroblast cells (ID1+)" = "#D98F5C",
  "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)" = "#C65A5A",
  "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)" = "#4C9F8A",
  "MMP3+ lining fibroblast cells (FAM184A+)" = "#9B7AAE",
  "HA-enriched SFRP2+ matrix fibroblast-like cells" = "#C7B24A",
  "Low cell count (<5 cells)" = "white"
)

condition_colors <- c(
  "other" = "#B65A5A",
  "HA" = "#5B8DB8"
)

marker_genes_fibroblast <- list(
  "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)" = c("HLA-DRA", "CD74", "CSN1S1", "AMTN"),
  "Activated MMP3+ lining fibroblast cells (ID1+)" = c("ID1", "FOS", "SERTAD1", "NFE2L3"),
  "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)" = c("CCL7", "CCRL2", "CCL20", "CXCL1", "BIRC3"),
  "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)" = c("ITGB8", "ADAMTSL1", "SEMA3A", "TXNIP", "ZNF385B"),
  "MMP3+ lining fibroblast cells (FAM184A+)" = c("FAM184A", "C5orf64", "ARHGAP15", "SERPINE3"),
  "HA-enriched SFRP2+ matrix fibroblast-like cells" = c("SFRP2", "SFRP1", "COMP", "PODN", "IGF1", "MFAP5")
)

IRON_FEATURES <- c(
  "HMOX1", "SLC40A1", "TFRC", "STEAP3", "FTH1", "FTL",
  "SLC11A2", "CP", "NCOA4", "GPX4", "ACSL4", "AIFM2",
  "NQO1", "GCLC", "GCLM", "ALOX5", "ALOX15", "SAT1",
  "SLC25A37", "SLC25A28", "SLC39A14", "SLC39A8"
)

build_label_colors <- function(labels, color_map) {
  labels <- unique(as.character(labels))
  cols <- color_map[labels]
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
  fill_levels <- unique(as.character(plot_df[[fill_col]]))
  count_df$cell_type <- factor(count_df$cell_type, levels = fill_levels)
  count_df
}

build_sample_ratio_plot_data <- function(meta_df, sample_col, condition_col, fill_col) {
  valid_rows <- !is.na(meta_df[[sample_col]]) & trimws(as.character(meta_df[[sample_col]])) != ""
  valid_rows <- valid_rows & !is.na(meta_df[[condition_col]]) & trimws(as.character(meta_df[[condition_col]])) != ""
  valid_rows <- valid_rows & !is.na(meta_df[[fill_col]]) & trimws(as.character(meta_df[[fill_col]])) != ""

  plot_df <- meta_df[valid_rows, c(sample_col, condition_col, fill_col), drop = FALSE]
  if (nrow(plot_df) == 0) return(NULL)

  colnames(plot_df) <- c("sample_id", "condition", "cell_type")
  plot_df$sample_id <- as.character(plot_df$sample_id)
  plot_df$condition <- ifelse(as.character(plot_df$condition) == "HA", "HA", "other")
  plot_df$condition <- factor(plot_df$condition, levels = c("HA", "other"))
  plot_df$cell_type <- as.character(plot_df$cell_type)

  count_df <- as.data.frame(table(plot_df$sample_id, plot_df$condition, plot_df$cell_type), stringsAsFactors = FALSE)
  colnames(count_df) <- c("sample_id", "condition", "cell_type", "n_cells")
  count_df <- count_df[count_df$n_cells > 0, , drop = FALSE]
  if (nrow(count_df) == 0) return(NULL)

  totals <- aggregate(n_cells ~ sample_id, data = count_df, FUN = sum)
  colnames(totals)[2] <- "sample_total"
  count_df <- merge(count_df, totals, by = "sample_id", sort = FALSE)
  count_df$ratio <- count_df$n_cells / count_df$sample_total

  sample_info <- unique(plot_df[, c("sample_id", "condition"), drop = FALSE])
  sample_info <- sample_info[order(sample_info$condition, sample_info$sample_id), , drop = FALSE]

  low_count_samples <- totals$sample_id[totals$sample_total < MIN_CELLS_PER_SAMPLE_RATIO_PLOT]
  if (length(low_count_samples) > 0) {
    low_count_rows <- merge(
      data.frame(sample_id = low_count_samples, stringsAsFactors = FALSE),
      unique(count_df[, c("sample_id", "condition", "sample_total"), drop = FALSE]),
      by = "sample_id",
      sort = FALSE
    )
    low_count_rows$cell_type <- LOW_CELL_COUNT_LABEL
    low_count_rows$n_cells <- low_count_rows$sample_total
    low_count_rows$ratio <- 1

    count_df <- count_df[!count_df$sample_id %in% low_count_samples, , drop = FALSE]
    count_df <- rbind(
      count_df[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total", "ratio"), drop = FALSE],
      low_count_rows[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total", "ratio"), drop = FALSE]
    )
  }

  count_df$is_low_cell_count <- count_df$sample_id %in% low_count_samples
  count_df$sample_id <- factor(count_df$sample_id, levels = sample_info$sample_id)
  count_df$condition <- factor(as.character(count_df$condition), levels = c("HA", "other"))
  count_df$cell_type <- factor(count_df$cell_type, levels = c(levels(meta_df[[fill_col]]), LOW_CELL_COUNT_LABEL))
  count_df
}

make_ratio_plot <- function(plot_df, fill_colors) {
  ggplot(plot_df, aes(x = group, y = ratio, fill = cell_type)) +
    geom_col(width = 0.92, color = NA) +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Ratio", fill = NULL) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.title.y = element_text(size = 18, color = "black"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length = grid::unit(0.22, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.4, "cm"),
      legend.box = "vertical",
      plot.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))
}

make_sample_ratio_plot <- function(plot_df, fill_colors) {
  sample_fill_colors <- fill_colors
  sample_fill_colors[LOW_CELL_COUNT_LABEL] <- "white"

  ggplot(plot_df, aes(x = sample_id, y = ratio, fill = cell_type)) +
    geom_col(width = 0.92, color = "black", linewidth = 0.25) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_manual(
      values = sample_fill_colors,
      breaks = setdiff(names(sample_fill_colors), LOW_CELL_COUNT_LABEL),
      drop = FALSE
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Ratio", fill = NULL) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.title.y = element_text(size = 18, color = "black"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length = grid::unit(0.22, "cm"),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 17, face = "bold", color = "black", margin = margin(t = 8)),
      panel.spacing.x = grid::unit(0.35, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.4, "cm"),
      legend.box = "vertical",
      plot.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 30, 20)
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))
}

# Load subset and updated full object
if (!file.exists(INPUT_SUBSET_OBJECT)) {
  stop("Missing destructive lining fibroblast subset object: ", INPUT_SUBSET_OBJECT)
}

if (!file.exists(INPUT_FULL_OBJECT)) {
  stop("Missing annotated full object with destructive lining fibroblast subclusters: ", INPUT_FULL_OBJECT)
}

cat("Loading destructive lining fibroblast subset object...\n")
obj_fibro <- readRDS(INPUT_SUBSET_OBJECT)
cat("Destructive lining fibroblast cells:", ncol(obj_fibro), "\n")

cat("Loading full annotated object...\n")
obj_full <- readRDS(INPUT_FULL_OBJECT)
cat("Full object cells:", ncol(obj_full), "\n")

if (!CLUSTER_COL %in% colnames(obj_fibro@meta.data)) {
  stop("Destructive lining fibroblast subcluster column not found in subset object: ", CLUSTER_COL)
}

if (!REDUCTION_NAME %in% names(obj_fibro@reductions)) {
  stop("Reduction not found in destructive lining fibroblast subset object: ", REDUCTION_NAME)
}

DefaultAssay(obj_fibro) <- "RNA"
obj_fibro <- JoinLayers(obj_fibro, assay = "RNA")

rna_data_layer <- tryCatch(
  LayerData(obj_fibro[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  cat("Normalizing RNA assay...\n")
  obj_fibro <- NormalizeData(obj_fibro, assay = "RNA", verbose = FALSE)
}

# Apply the final subcluster-to-state annotation
cluster_ids <- as.character(obj_fibro@meta.data[[CLUSTER_COL]])
mapped_labels <- unname(DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS[cluster_ids])
mapped_labels[is.na(mapped_labels)] <- "Unknown"

obj_fibro@meta.data[[LABEL_COL]] <- factor(
  mapped_labels,
  levels = unname(DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS)
)
obj_fibro@meta.data[[LABEL_COL]] <- droplevels(obj_fibro@meta.data[[LABEL_COL]])

cat("Destructive lining fibroblast subtype labels after mapping:\n")
print(table(obj_fibro@meta.data[[LABEL_COL]], useNA = "ifany"))

if (all(is.na(obj_fibro@meta.data[[LABEL_COL]]))) {
  stop("All destructive lining fibroblast subtype labels are NA after mapping. Check ", CLUSTER_COL, " and label mapping.")
}

obj_full@meta.data[[LABEL_COL]] <- NA_character_
obj_full@meta.data[Cells(obj_fibro), LABEL_COL] <- as.character(obj_fibro@meta.data[[LABEL_COL]])

# Save the final annotation table
cluster_levels <- unique(cluster_ids)
suppressWarnings(cluster_levels_num <- as.integer(cluster_levels))
if (all(!is.na(cluster_levels_num))) {
  cluster_levels <- as.character(sort(cluster_levels_num))
} else {
  cluster_levels <- sort(cluster_levels)
}

cluster_annotation_table <- data.frame(
  destructive_lining_fibroblast_subcluster = cluster_levels,
  destructive_lining_fibroblast_subtype = unname(DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS[cluster_levels]),
  n_cells = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(cluster_annotation_table))) {
  cluster_annotation_table$n_cells[i] <- sum(cluster_ids == cluster_annotation_table$destructive_lining_fibroblast_subcluster[i])
}

write.csv(
  cluster_annotation_table,
  file.path(RES_DIR, "destructive_lining_fibroblast_subcluster_annotation_table.csv"),
  row.names = FALSE
)

# Save the canonical markers used for the dotplot
marker_rows <- list()
marker_index <- 1

for (ct in names(marker_genes_fibroblast)) {
  genes_here <- unique(marker_genes_fibroblast[[ct]])
  genes_here <- genes_here[genes_here %in% rownames(obj_fibro)]
  if (length(genes_here) == 0) next

  marker_rows[[marker_index]] <- data.frame(
    destructive_lining_fibroblast_subtype = rep(ct, length(genes_here)),
    gene = genes_here,
    stringsAsFactors = FALSE
  )
  marker_index <- marker_index + 1
}

canonical_marker_table <- do.call(rbind, marker_rows)
write.csv(
  canonical_marker_table,
  file.path(RES_DIR, "destructive_lining_fibroblast_canonical_marker_table.csv"),
  row.names = FALSE
)

# Save the updated full object with final fibroblast state labels
saveRDS(obj_full, OUTPUT_FULL_OBJECT)

cat("\nPlotting destructive lining fibroblast UMAPs...\n")
fibro_levels <- levels(obj_fibro@meta.data[[LABEL_COL]])
fibro_color_map <- build_label_colors(fibro_levels, destructive_lining_fibroblast_subtype_colors)

p_umap <- DimPlot(
  object = obj_fibro,
  reduction = REDUCTION_NAME,
  group.by = LABEL_COL,
  raster = FALSE,
  pt.size = PT_SIZE
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
  filename = file.path(FIG_DIR, "umap_by_destructive_lining_fibroblast_subtype.png"),
  plot = p_umap,
  width = 12,
  height = 6.5,
  dpi = 600
)

# UMAP colored directly by HA vs other condition
if ("condition" %in% colnames(obj_fibro@meta.data)) {
  condition_values <- ifelse(as.character(obj_fibro@meta.data$condition) == "HA", "HA", "other")
  obj_fibro@meta.data$condition <- factor(condition_values, levels = c("HA", "other"))
  condition_color_map <- condition_colors[levels(obj_fibro@meta.data$condition)]

  p_umap_condition <- DimPlot(
    object = obj_fibro,
    reduction = REDUCTION_NAME,
    group.by = "condition",
    raster = FALSE,
    pt.size = PT_SIZE
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
    filename = file.path(FIG_DIR, "umap_by_condition.png"),
    plot = p_umap_condition,
    width = 9,
    height = 6.5,
    dpi = 600
  )
}

# Plot split UMAPs by condition and condition_all
for (split_col in c("condition", "condition_all")) {
  if (split_col %in% colnames(obj_fibro@meta.data)) {
    split_values <- as.character(obj_fibro@meta.data[[split_col]])
    split_values[is.na(split_values) | trimws(split_values) == ""] <- "NA"
    obj_fibro@meta.data[[split_col]] <- factor(split_values, levels = unique(split_values))
    obj_fibro@meta.data[[split_col]] <- droplevels(obj_fibro@meta.data[[split_col]])

    p_umap_split <- DimPlot(
      object = obj_fibro,
      reduction = REDUCTION_NAME,
      group.by = LABEL_COL,
      split.by = split_col,
      raster = FALSE,
      pt.size = PT_SIZE
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
      filename = file.path(FIG_DIR, paste0("umap_by_destructive_lining_fibroblast_subtype_split_", split_col, ".png")),
      plot = p_umap_split,
      width = 13,
      height = 6.5,
      dpi = 600
    )
  }
}

# Plot the two usual composition barplots by condition and condition_all
for (x_col in c("condition", "condition_all")) {
  if (!x_col %in% colnames(obj_fibro@meta.data)) next

  plot_df <- build_ratio_plot_data(
    meta_df = obj_fibro@meta.data,
    x_col = x_col,
    fill_col = LABEL_COL
  )

  if (is.null(plot_df)) next

  p_ratio <- make_ratio_plot(plot_df, fibro_color_map)

  ggsave(
    filename = file.path(FIG_DIR, paste0("destructive_lining_fibroblast_subtype_ratio_by_", x_col, ".png")),
    plot = p_ratio,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# Plot composition per sample with HA and other separated in the x-axis lower strip
if (all(c("sample_id", "condition") %in% colnames(obj_fibro@meta.data))) {
  sample_ratio_df <- build_sample_ratio_plot_data(
    meta_df = obj_fibro@meta.data,
    sample_col = "sample_id",
    condition_col = "condition",
    fill_col = LABEL_COL
  )

  if (!is.null(sample_ratio_df)) {
    write.csv(
      sample_ratio_df[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total", "ratio", "is_low_cell_count")],
      file.path(RES_DIR, "destructive_lining_fibroblast_subtype_ratio_by_sample.csv"),
      row.names = FALSE
    )

    p_sample_ratio <- make_sample_ratio_plot(sample_ratio_df, fibro_color_map)

    ggsave(
      filename = file.path(FIG_DIR, "destructive_lining_fibroblast_subtype_ratio_by_sample_grouped_condition.png"),
      plot = p_sample_ratio,
      width = 13,
      height = 6.5,
      dpi = 600
    )
  }
}

# Plot the canonical marker dotplot for final fibroblast states
marker_features <- unlist(marker_genes_fibroblast, use.names = FALSE)
marker_features <- unique(marker_features[marker_features %in% rownames(obj_fibro)])

if (length(marker_features) > 0) {
  p_dot <- DotPlot(
    object = obj_fibro,
    features = marker_features,
    group.by = LABEL_COL,
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
    filename = file.path(FIG_DIR, "dotplot_canonical_markers_by_destructive_lining_fibroblast_subtype.png"),
    plot = p_dot,
    width = 15,
    height = 6.8,
    dpi = 600
  )
}

# Plot iron-related FeaturePlots to inspect iron-handling states
iron_features_use <- IRON_FEATURES[IRON_FEATURES %in% rownames(obj_fibro)]

if (length(iron_features_use) > 0) {
  p_feature <- FeaturePlot(
    object = obj_fibro,
    features = iron_features_use,
    reduction = REDUCTION_NAME,
    raster = FALSE,
    order = TRUE,
    cols = c("#F3F1EC", "#7A1F2B"),
    ncol = 4
  )

  ggsave(
    filename = file.path(FIG_DIR, "featureplot_iron_related_genes_destructive_lining_fibroblast_subtype.png"),
    plot = p_feature,
    width = 15,
    height = 11,
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Final destructive lining fibroblast annotation complete.\n")
cat("Updated object: ", OUTPUT_FULL_OBJECT, "\n")
cat("Results dir   : ", RES_DIR, "\n")
cat("Figures dir   : ", FIG_DIR, "\n")
cat("Batch effect  : mitigated in the subset by Harmony on sample_id.\n")
cat("============================================================\n")
