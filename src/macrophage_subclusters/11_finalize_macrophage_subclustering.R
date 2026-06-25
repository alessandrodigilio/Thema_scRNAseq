###########################################
### Final macrophage subtype annotation ###
###########################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# work from the project root
setwd("~/Thema_R")
source("src/global_config.R")

# create output directories used by this step
input_subset_object <- file.path(data_dir, "integrated_object", "macrophages_subclustered.rds")
input_full_object <- file.path(data_dir, "integrated_object", "annotated_macrophage_subclusters.rds")

out_dir <- file.path(data_dir, "integrated_object")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(figures_dir, "macrophage_final_annotation")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

res_dir <- file.path(results_dir, "macrophage_final_annotation")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

output_full_object <- file.path(out_dir, "annotated_macrophage_states.rds")

# set parameters
cluster_col <- "macrophage_subcluster"
label_col <- "macrophage_subtype"
reduction_name <- "umap.harmony.macrophage"
pt_size <- 0.14

macrophage_subcluster_labels <- c(
  "0" = "Inflammatory macrophages (KANK1+)",
  "1" = "Inflammatory macrophages (THBS1+)",
  "2" = "Macrophage-like state (AMTN+)",
  "3" = "Resident macrophages (HSPA6+)",
  "4" = "Red-pulp-like resident macrophages (MERTK+)",
  "5" = "Mixed macrophage-like cells (RNASE1+)",
  "6" = "Plasma-like contaminants",
  "7" = "Low-confidence cells",
  "8" = "Proliferating macrophages"
)

macrophage_subtype_colors <- c(
  "Inflammatory macrophages (KANK1+)" = "#C65A5A",
  "Inflammatory macrophages (THBS1+)" = "#b0e17b",
  "Macrophage-like state (AMTN+)" = "#D98F5C",
  "Resident macrophages (HSPA6+)" = "#8c674b",
  "Red-pulp-like resident macrophages (MERTK+)" = "#4C9F8A",
  "Mixed macrophage-like cells (RNASE1+)" = "#5f90b3",
  "Plasma-like contaminants" = "#B58ACF",
  "Low-confidence cells" = "#9FA4A9",
  "Proliferating macrophages" = "#D95FA7"
)

marker_genes_macrophage <- list(
  "Inflammatory macrophages (KANK1+)" = c("KANK1", "KHDRBS3", "BAALC"),
  "Inflammatory macrophages (THBS1+)" = c("THBS1", "IL2RA", "TNIP3"),
  "Macrophage-like state (AMTN+)" = c("AMTN", "SULF1", "ITGBL1"),
  "Resident macrophages (HSPA6+)" = c("HSPA6", "HSPE1-MOB4", "ENSG00000293472"),
  "Red-pulp-like resident macrophages (MERTK+)" = c("MERTK", "CD163", "FCGR3A"),
  "Mixed macrophage-like cells (RNASE1+)" = c("RNASE1", "CRIP1", "HCST"),
  "Plasma-like contaminants" = c("IGHA1", "MZB1", "IGKC"),
  "Low-confidence cells" = c("CCDC200", "RNU5B-1", "ENSG00000285646"),
  "Proliferating macrophages" = c("MKI67", "UBE2C", "CDCA3")
)

iron_features <- c(
  "SPIC", "HMOX1", "SLC40A1", "TFRC", "STEAP3", "FTH1", "FTL",
  "SLC11A2", "CP", "NCOA4", "GPX4", "ACSL4", "AIFM2",
  "NQO1", "GCLC", "GCLM", "ALOX5", "ALOX15", "SAT1"
)

red_pulp_like_features <- c(
  "CD163",    # paper isolation strategy: CD163 high red pulp macrophages
  "CD14",     # paper isolation strategy: CD14 low compared with monocytes
  "FCGR2A",   # fc gamma receptor IIA, reported as part of the RPM Fc receptor signature
  "FCGR3A",   # fc gamma receptor IIIA, reported as part of the RPM Fc receptor signature
  "FCGR2B",   # fc gamma receptor IIB, useful as low/absent comparator from the paper
  "CYBB",     # gp91PHOX, reported as highly expressed in human RPM
  "SLC48A1",  # heme transporter HRG-1 from phagolysosome to cytosol
  "HMOX1",    # heme degradation enzyme HO-1 and iron recycling
  "SLC40A1"   # ferroportin, iron export
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

# load macrophage subset and updated full object
if (!file.exists(input_subset_object)) {
  stop("Missing macrophage subset object: ", input_subset_object)
}

if (!file.exists(input_full_object)) {
  stop("Missing annotated full object with macrophage subclusters: ", input_full_object)
}

cat("Loading macrophage subset object...\n")
obj_macro <- readRDS(input_subset_object)
cat("Macrophage cells:", ncol(obj_macro), "\n")

cat("Loading full annotated object...\n")
obj_full <- readRDS(input_full_object)
cat("Full object cells:", ncol(obj_full), "\n")

if (!cluster_col %in% colnames(obj_macro@meta.data)) {
  stop("Macrophage subcluster column not found in subset object: ", cluster_col)
}

if (!reduction_name %in% names(obj_macro@reductions)) {
  stop("Reduction not found in macrophage subset object: ", reduction_name)
}

DefaultAssay(obj_macro) <- "RNA"
obj_macro <- JoinLayers(obj_macro, assay = "RNA")

rna_data_layer <- tryCatch(
  LayerData(obj_macro[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  cat("Normalizing RNA assay...\n")
  obj_macro <- NormalizeData(obj_macro, assay = "RNA", verbose = FALSE)
}

# apply the final subcluster-to-state annotation
cluster_ids <- as.character(obj_macro@meta.data[[cluster_col]])
mapped_labels <- unname(macrophage_subcluster_labels[cluster_ids])
mapped_labels[is.na(mapped_labels)] <- "Unknown"

# store the final labels directly in metadata before plotting
obj_macro@meta.data[[label_col]] <- factor(
  mapped_labels,
  levels = unname(macrophage_subcluster_labels)
)
obj_macro@meta.data[[label_col]] <- droplevels(obj_macro@meta.data[[label_col]])

cat("Macrophage subtype labels after mapping:\n")
print(table(obj_macro@meta.data[[label_col]], useNA = "ifany"))

if (all(is.na(obj_macro@meta.data[[label_col]]))) {
  stop("All macrophage subtype labels are NA after mapping. Check ", cluster_col, " and label mapping.")
}

obj_full@meta.data[[label_col]] <- NA_character_
obj_full@meta.data[Cells(obj_macro), label_col] <- as.character(obj_macro@meta.data[[label_col]])

# save the final macrophage annotation table
cluster_levels <- unique(cluster_ids)
suppressWarnings(cluster_levels_num <- as.integer(cluster_levels))
if (all(!is.na(cluster_levels_num))) {
  cluster_levels <- as.character(sort(cluster_levels_num))
} else {
  cluster_levels <- sort(cluster_levels)
}

cluster_annotation_table <- data.frame(
  macrophage_subcluster = cluster_levels,
  macrophage_subtype = unname(macrophage_subcluster_labels[cluster_levels]),
  n_cells = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(cluster_annotation_table))) {
  cluster_annotation_table$n_cells[i] <- sum(cluster_ids == cluster_annotation_table$macrophage_subcluster[i])
}

write.csv(cluster_annotation_table, file.path(res_dir, "macrophage_subcluster_annotation_table.csv"), row.names = FALSE)

# save the canonical markers used for the macrophage dotplot
marker_rows <- list()
marker_index <- 1

for (ct in names(marker_genes_macrophage)) {
  genes_here <- unique(marker_genes_macrophage[[ct]])
  genes_here <- genes_here[genes_here %in% rownames(obj_macro)]
  if (length(genes_here) == 0) next

  marker_rows[[marker_index]] <- data.frame(
    macrophage_subtype = rep(ct, length(genes_here)),
    gene = genes_here,
    stringsAsFactors = FALSE
  )
  marker_index <- marker_index + 1
}

canonical_marker_table <- do.call(rbind, marker_rows)
write.csv(canonical_marker_table, file.path(res_dir, "macrophage_canonical_marker_table.csv"), row.names = FALSE)

# save the updated full object with the final macrophage state labels
saveRDS(obj_full, output_full_object)

# plot the macrophage UMAP colored by final macrophage state
cat("\nPlotting macrophage UMAP...\n")
macro_levels <- levels(obj_macro@meta.data[[label_col]])
macro_color_map <- build_label_colors(macro_levels, macrophage_subtype_colors)

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

# plot split UMAPs by condition and condition_all
for (split_col in c("condition", "condition_all")) {
  if (split_col %in% colnames(obj_macro@meta.data)) {
    split_values <- as.character(obj_macro@meta.data[[split_col]])
    split_values[is.na(split_values) | trimws(split_values) == ""] <- "NA"
    obj_macro@meta.data[[split_col]] <- factor(split_values, levels = unique(split_values))
    obj_macro@meta.data[[split_col]] <- droplevels(obj_macro@meta.data[[split_col]])

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
}

# plot the two usual composition barplots
for (x_col in c("condition", "condition_all")) {
  if (!x_col %in% colnames(obj_macro@meta.data)) next

  plot_df <- build_ratio_plot_data(
    meta_df = obj_macro@meta.data,
    x_col = x_col,
    fill_col = label_col
  )

  if (is.null(plot_df)) next

  p_ratio <- make_ratio_plot(plot_df, unname(macro_color_map))

  ggsave(
    filename = file.path(fig_dir, paste0("macrophage_subtype_ratio_by_", x_col, ".png")),
    plot = p_ratio,
    width = 12,
    height = 6,
    dpi = 600
  )
}

# plot the canonical marker bubble plot for the final macrophage states
marker_features <- unlist(marker_genes_macrophage, use.names = FALSE)
marker_features <- unique(marker_features[marker_features %in% rownames(obj_macro)])
dotplot_df <- build_marker_dotplot_data(
  object = obj_macro,
  features = marker_features,
  group_col = label_col
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
    filename = file.path(fig_dir, "dotplot_canonical_markers_by_macrophage_subtype.png"),
    plot = p_dot,
    width = 15,
    height = 6.5,
    dpi = 600
  )
}

# plot iron-related FeaturePlots to inspect iron-handling states
iron_features_use <- iron_features[iron_features %in% rownames(obj_macro)]

if (length(iron_features_use) > 0) {
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
}

# plot red-pulp-like / iron-recycling markers separately
red_pulp_features_use <- red_pulp_like_features[red_pulp_like_features %in% rownames(obj_macro)]

if (length(red_pulp_features_use) > 0) {
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
}

# plot a dedicated red-pulp-like marker bubble plot across macrophage subtypes
red_pulp_dotplot_df <- build_marker_dotplot_data(
  object = obj_macro,
  features = red_pulp_features_use,
  group_col = label_col
)

if (!is.null(red_pulp_dotplot_df) && nrow(red_pulp_dotplot_df) > 0) {
  p_dot_red_pulp <- ggplot(
    red_pulp_dotplot_df,
    aes(x = gene, y = group, size = pct_expr, color = avg_scaled)
  ) +
    geom_point() +
    scale_size_continuous(name = "Percent Expressed", range = c(0.6, 6)) +
    scale_color_gradient2(
      name = "Average Expression",
      low = "#F4EDF8",
      mid = "#E4D4F1",
      high = "#6F2DBD",
      midpoint = 0,
      limits = c(min(red_pulp_dotplot_df$avg_scaled), max(red_pulp_dotplot_df$avg_scaled))
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
    filename = file.path(fig_dir, "dotplot_red_pulp_like_markers_by_macrophage_subtype.png"),
    plot = p_dot_red_pulp,
    width = 12,
    height = 6.5,
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Final macrophage annotation complete.\n")
cat("Updated object: ", output_full_object, "\n")
cat("Results dir   : ", res_dir, "\n")
cat("Figures dir   : ", fig_dir, "\n")
cat("Batch effect  : mitigated in the subset by Harmony on sample_id.\n")
cat("============================================================\n")
