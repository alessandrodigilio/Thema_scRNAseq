# ---------------------------------------
# 17_GEX_finalize_activated_endothelial_subclustering.R
# Final activated endothelial subcluster annotation and plotting
# ---------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
INPUT_SUBSET_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_activated_endothelial_subclustered.rds")
INPUT_FULL_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_annotated_endothelial_subclusters.rds")

OUT_DIR <- file.path(DATA_DIR, "integrated_object")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(FIGURES_DIR, "activated_endothelial_final_annotation")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RES_DIR <- file.path(RESULTS_DIR, "activated_endothelial_final_annotation")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_FULL_OBJECT <- file.path(OUT_DIR, "gex_annotated_endothelial_states.rds")

# Set parameters
CLUSTER_COL <- "endothelial_subcluster"
LABEL_COL <- "endothelial_subtype"
REDUCTION_NAME <- "umap.harmony.endothelial"
PT_SIZE <- 0.14

ENDOTHELIAL_SUBCLUSTER_LABELS <- c(
  "0" = "Endothelial cells (PLXNA4+)",
  "1" = "Stress-response endothelial cells (HSPA6+)",
  "2" = "Activated endothelial cells (IL6+)",
  "3" = "Endothelial cells (ZNF385B+)",
  "4" = "Endothelial cells (EDNRB+)",
  "5" = "Arterial-like endothelial cells (GJA5+)",
  "6" = "Endothelial cells (SLC2A14+)",
  "7" = "Mixed stromal-like cells",
  "8" = "Mural-like cells"
)

endothelial_subtype_colors <- c(
  "Endothelial cells (PLXNA4+)" = "#D46A6A",
  "Stress-response endothelial cells (HSPA6+)" = "#C7B24A",
  "Activated endothelial cells (IL6+)" = "#84B547",
  "Endothelial cells (ZNF385B+)" = "#7C7FD4",
  "Endothelial cells (EDNRB+)" = "#49B6A5",
  "Arterial-like endothelial cells (GJA5+)" = "#4E88C7",
  "Endothelial cells (SLC2A14+)" = "#57C6D9",
  "Mixed stromal-like cells" = "#C58A8A",
  "Mural-like cells" = "#C86DD7"
)

marker_genes_endothelial <- list(
  "Endothelial cells (PLXNA4+)" = c("PLXNA4", "NR2F2-AS1", "PLCXD3"),
  "Stress-response endothelial cells (HSPA6+)" = c("HSPA6", "G0S2", "MMP1"),
  "Activated endothelial cells (IL6+)" = c("IL6", "RGS16", "HES1"),
  "Endothelial cells (ZNF385B+)" = c("ZNF385B", "LINC01411", "KCNQ5"),
  "Endothelial cells (EDNRB+)" = c("EDNRB", "BTNL9", "RCAN2"),
  "Arterial-like endothelial cells (GJA5+)" = c("GJA5", "TMEM178A", "LINC00639"),
  "Endothelial cells (SLC2A14+)" = c("SLC2A14", "SNX31", "LRRC23"),
  "Mixed stromal-like cells" = c("COMP", "IGF1", "COL1A1"),
  "Mural-like cells" = c("RGS5", "AVPR1A", "ABCC9")
)

IRON_FEATURES <- c(
  "HMOX1", "SLC40A1", "TFRC", "STEAP3", "FTH1", "FTL",
  "SLC11A2", "CP", "NCOA4", "GPX4", "ACSL4", "AIFM2",
  "NQO1", "GCLC", "GCLM", "ALOX5", "ALOX15", "SAT1"
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

build_marker_dotplot_data <- function(object, features, group_col) {
  if (length(features) == 0) return(NULL)

  plot_df <- FetchData(object, vars = c(group_col, features))
  colnames(plot_df)[1] <- "group"
  plot_df$group <- as.character(plot_df$group)
  plot_df <- plot_df[!is.na(plot_df$group) & plot_df$group != "", , drop = FALSE]
  if (nrow(plot_df) == 0) return(NULL)

  group_levels <- levels(object@meta.data[[group_col]])
  if (is.null(group_levels)) {
    group_levels <- unique(plot_df$group)
  }

  res_list <- vector("list", length(features))

  for (i in seq_along(features)) {
    g <- features[i]
    df_g <- data.frame(
      group = plot_df$group,
      expr = plot_df[[g]],
      stringsAsFactors = FALSE
    )

    avg_expr <- tapply(df_g$expr, df_g$group, mean, na.rm = TRUE)
    pct_expr <- tapply(df_g$expr > 0, df_g$group, mean, na.rm = TRUE) * 100

    avg_expr <- avg_expr[group_levels]
    pct_expr <- pct_expr[group_levels]

    avg_expr[is.na(avg_expr)] <- 0
    pct_expr[is.na(pct_expr)] <- 0

    if (stats::sd(avg_expr) > 0) {
      avg_scaled <- as.numeric(scale(avg_expr))
    } else {
      avg_scaled <- rep(0, length(avg_expr))
    }

    res_list[[i]] <- data.frame(
      group = group_levels,
      gene = g,
      avg_scaled = avg_scaled,
      pct_expr = as.numeric(pct_expr),
      stringsAsFactors = FALSE
    )
  }

  res_df <- do.call(rbind, res_list)
  res_df <- res_df[res_df$pct_expr > 0, , drop = FALSE]
  if (nrow(res_df) == 0) return(NULL)

  res_df$group <- factor(res_df$group, levels = group_levels)
  res_df$gene <- factor(res_df$gene, levels = features)
  res_df
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

# Load subset and updated full object
if (!file.exists(INPUT_SUBSET_OBJECT)) {
  stop("Missing endothelial subset object: ", INPUT_SUBSET_OBJECT)
}

if (!file.exists(INPUT_FULL_OBJECT)) {
  stop("Missing annotated full object with endothelial subclusters: ", INPUT_FULL_OBJECT)
}

cat("Loading activated endothelial subset object...\n")
obj_endo <- readRDS(INPUT_SUBSET_OBJECT)
cat("Activated endothelial cells:", ncol(obj_endo), "\n")

cat("Loading full annotated object...\n")
obj_full <- readRDS(INPUT_FULL_OBJECT)
cat("Full object cells:", ncol(obj_full), "\n")

if (!CLUSTER_COL %in% colnames(obj_endo@meta.data)) {
  stop("Endothelial subcluster column not found in subset object: ", CLUSTER_COL)
}

if (!REDUCTION_NAME %in% names(obj_endo@reductions)) {
  stop("Reduction not found in endothelial subset object: ", REDUCTION_NAME)
}

DefaultAssay(obj_endo) <- "RNA"
obj_endo <- JoinLayers(obj_endo, assay = "RNA")

rna_data_layer <- tryCatch(
  LayerData(obj_endo[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  cat("Normalizing RNA assay...\n")
  obj_endo <- NormalizeData(obj_endo, assay = "RNA", verbose = FALSE)
}

# Apply final subcluster-to-state annotation
cluster_ids <- as.character(obj_endo@meta.data[[CLUSTER_COL]])
mapped_labels <- unname(ENDOTHELIAL_SUBCLUSTER_LABELS[cluster_ids])
mapped_labels[is.na(mapped_labels)] <- "Unknown"

obj_endo@meta.data[[LABEL_COL]] <- factor(
  mapped_labels,
  levels = unname(ENDOTHELIAL_SUBCLUSTER_LABELS)
)
obj_endo@meta.data[[LABEL_COL]] <- droplevels(obj_endo@meta.data[[LABEL_COL]])

cat("Endothelial subtype labels after mapping:\n")
print(table(obj_endo@meta.data[[LABEL_COL]], useNA = "ifany"))

if (all(is.na(obj_endo@meta.data[[LABEL_COL]]))) {
  stop("All endothelial subtype labels are NA after mapping. Check ", CLUSTER_COL, " and label mapping.")
}

obj_full@meta.data[[LABEL_COL]] <- NA_character_
obj_full@meta.data[Cells(obj_endo), LABEL_COL] <- as.character(obj_endo@meta.data[[LABEL_COL]])

# Save the final endothelial annotation table
cluster_levels <- unique(cluster_ids)
suppressWarnings(cluster_levels_num <- as.integer(cluster_levels))
if (all(!is.na(cluster_levels_num))) {
  cluster_levels <- as.character(sort(cluster_levels_num))
} else {
  cluster_levels <- sort(cluster_levels)
}

cluster_annotation_table <- data.frame(
  endothelial_subcluster = cluster_levels,
  endothelial_subtype = unname(ENDOTHELIAL_SUBCLUSTER_LABELS[cluster_levels]),
  n_cells = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(cluster_annotation_table))) {
  cluster_annotation_table$n_cells[i] <- sum(cluster_ids == cluster_annotation_table$endothelial_subcluster[i])
}

write.csv(
  cluster_annotation_table,
  file.path(RES_DIR, "activated_endothelial_subcluster_annotation_table.csv"),
  row.names = FALSE
)

# Save the canonical markers used for the endothelial dotplot
marker_rows <- list()
marker_index <- 1

for (ct in names(marker_genes_endothelial)) {
  genes_here <- unique(marker_genes_endothelial[[ct]])
  genes_here <- genes_here[genes_here %in% rownames(obj_endo)]
  if (length(genes_here) == 0) next

  marker_rows[[marker_index]] <- data.frame(
    endothelial_subtype = rep(ct, length(genes_here)),
    gene = genes_here,
    stringsAsFactors = FALSE
  )
  marker_index <- marker_index + 1
}

canonical_marker_table <- do.call(rbind, marker_rows)
write.csv(
  canonical_marker_table,
  file.path(RES_DIR, "activated_endothelial_canonical_marker_table.csv"),
  row.names = FALSE
)

# Save the updated full object with the final endothelial state labels
saveRDS(obj_full, OUTPUT_FULL_OBJECT)

# Plot the endothelial UMAP colored by final endothelial state
cat("\nPlotting activated endothelial UMAP...\n")
endo_levels <- levels(obj_endo@meta.data[[LABEL_COL]])
endo_color_map <- build_label_colors(endo_levels, endothelial_subtype_colors)

p_umap <- DimPlot(
  object = obj_endo,
  reduction = REDUCTION_NAME,
  group.by = LABEL_COL,
  raster = FALSE,
  pt.size = PT_SIZE
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

ggsave(
  filename = file.path(FIG_DIR, "umap_by_endothelial_subtype.png"),
  plot = p_umap,
  width = 11,
  height = 6.5,
  dpi = 600
)

# Plot split UMAPs by condition and condition_all
for (split_col in c("condition", "condition_all")) {
  if (split_col %in% colnames(obj_endo@meta.data)) {
    split_values <- as.character(obj_endo@meta.data[[split_col]])
    split_values[is.na(split_values) | trimws(split_values) == ""] <- "NA"
    obj_endo@meta.data[[split_col]] <- factor(split_values, levels = unique(split_values))
    obj_endo@meta.data[[split_col]] <- droplevels(obj_endo@meta.data[[split_col]])

    p_umap_split <- DimPlot(
      object = obj_endo,
      reduction = REDUCTION_NAME,
      group.by = LABEL_COL,
      split.by = split_col,
      raster = FALSE,
      pt.size = PT_SIZE
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

    ggsave(
      filename = file.path(FIG_DIR, paste0("umap_by_endothelial_subtype_split_", split_col, ".png")),
      plot = p_umap_split,
      width = 13,
      height = 6.5,
      dpi = 600
    )
  }
}

# Plot the two usual composition barplots
for (x_col in c("condition", "condition_all")) {
  if (!x_col %in% colnames(obj_endo@meta.data)) next

  plot_df <- build_ratio_plot_data(
    meta_df = obj_endo@meta.data,
    x_col = x_col,
    fill_col = LABEL_COL
  )

  if (is.null(plot_df)) next

  p_ratio <- make_ratio_plot(plot_df, unname(endo_color_map))

  ggsave(
    filename = file.path(FIG_DIR, paste0("endothelial_subtype_ratio_by_", x_col, ".png")),
    plot = p_ratio,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# Plot the canonical marker bubble plot for the final endothelial states
marker_features <- unlist(marker_genes_endothelial, use.names = FALSE)
marker_features <- unique(marker_features[marker_features %in% rownames(obj_endo)])
dotplot_df <- build_marker_dotplot_data(
  object = obj_endo,
  features = marker_features,
  group_col = LABEL_COL
)

if (!is.null(dotplot_df) && nrow(dotplot_df) > 0) {
  p_dot <- ggplot(dotplot_df, aes(x = gene, y = group, size = pct_expr, color = avg_scaled)) +
    geom_point() +
    scale_size_continuous(name = "Percent Expressed", range = c(0.6, 6)) +
    scale_color_gradient2(
      name = "Average Expression",
      low = "#F3F1EC",
      mid = "#F3F1EC",
      high = "#d84b1c",
      midpoint = 0,
      limits = c(min(dotplot_df$avg_scaled), max(dotplot_df$avg_scaled))
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
      legend.key.width = grid::unit(0.45, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      plot.margin = margin(20, 35, 35, 25)
    ) +
    guides(
      color = guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(2.4, "cm"),
        barheight = grid::unit(0.35, "cm")
      ),
      size = guide_legend(title.position = "top")
    )

  ggsave(
    filename = file.path(FIG_DIR, "dotplot_canonical_markers_by_endothelial_subtype.png"),
    plot = p_dot,
    width = 15,
    height = 6.5,
    dpi = 600
  )
}

# Plot iron-related FeaturePlots to inspect iron-handling states
iron_features_use <- IRON_FEATURES[IRON_FEATURES %in% rownames(obj_endo)]

if (length(iron_features_use) > 0) {
  p_feature_iron <- FeaturePlot(
    object = obj_endo,
    features = iron_features_use,
    reduction = REDUCTION_NAME,
    raster = FALSE,
    order = TRUE,
    cols = c("#F3F1EC", "#7A1F2B"),
    ncol = 4
  )

  ggsave(
    filename = file.path(FIG_DIR, "featureplot_iron_related_genes_endothelial_subtype.png"),
    plot = p_feature_iron,
    width = 15,
    height = 10.5,
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Final activated endothelial annotation complete.\n")
cat("Updated object: ", OUTPUT_FULL_OBJECT, "\n")
cat("Results dir   : ", RES_DIR, "\n")
cat("Figures dir   : ", FIG_DIR, "\n")
cat("Batch effect  : mitigated in the subset by Harmony on sample_id.\n")
cat("============================================================\n")
