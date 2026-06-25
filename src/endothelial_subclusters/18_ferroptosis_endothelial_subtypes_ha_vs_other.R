############################################################
### Endothelial subtype ferroptosis scoring: HA vs other ###
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(readxl)
})

# work from the project root
setwd("~/Thema_R")
source("src/global_config.R")

# create output directories used by this step
input_object <- file.path(data_dir, "integrated_object", "activated_endothelial_subclustered.rds")
geneset_file <- ferroptosis_geneset_file

fig_dir <- file.path(figures_dir, "endothelial_ferroptosis")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

res_dir <- file.path(results_dir, "endothelial_ferroptosis")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
cluster_col <- "endothelial_subcluster"
label_col <- "endothelial_subtype"
reduction_name <- "umap.harmony.endothelial"
group_col <- "condition"
ha_label <- "HA"
pt_size <- 0.14
ferroptosis_low_color <- "#F5E6E8"
ferroptosis_high_color <- "#7A1F2B"

endothelial_subcluster_labels <- c(
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

wrap_each_word <- function(x) {
  x <- gsub(" / ", "/", as.character(x), fixed = TRUE)
  vapply(strsplit(x, " "), function(words) {
    paste(words, collapse = "\n")
  }, character(1))
}

if (!file.exists(input_object)) {
  stop("Missing endothelial subset object: ", input_object)
}

if (!file.exists(geneset_file)) {
  stop("Missing ferroptosis gene set: ", geneset_file)
}

cat("Loading activated endothelial subset object...\n")
obj <- readRDS(input_object)
cat("Cells:", ncol(obj), "\n")

if (!cluster_col %in% colnames(obj@meta.data)) {
  stop("Endothelial subcluster column not found in metadata: ", cluster_col)
}

if (!reduction_name %in% names(obj@reductions)) {
  stop("Reduction not found in object: ", reduction_name)
}

if (!group_col %in% colnames(obj@meta.data)) {
  stop("Group column not found in metadata: ", group_col)
}

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")

rna_data_layer <- tryCatch(
  LayerData(obj[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  cat("Normalizing RNA assay...\n")
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

# add endothelial subtype labels from the finalized annotation
cluster_ids <- as.character(obj@meta.data[[cluster_col]])
mapped_labels <- unname(endothelial_subcluster_labels[cluster_ids])
mapped_labels[is.na(mapped_labels)] <- "Unknown"
obj@meta.data[[label_col]] <- factor(
  mapped_labels,
  levels = unname(endothelial_subcluster_labels)
)
obj@meta.data[[label_col]] <- droplevels(obj@meta.data[[label_col]])

cat("Endothelial subtype labels:\n")
print(table(obj@meta.data[[label_col]], useNA = "ifany"))

genes_df <- as.data.frame(read_excel(geneset_file))
genes_vec <- unique(trimws(as.character(genes_df[[1]])))
genes_vec <- genes_vec[!is.na(genes_vec) & genes_vec != ""]
genes_vec[genes_vec == "TRFC"] <- "TFRC"
genes_use <- intersect(genes_vec, rownames(obj))

if (length(genes_use) == 0) {
  stop("No ferroptosis genes from Excel were found in the endothelial subset")
}

cat("Ferroptosis genes in Excel:", length(genes_vec), "\n")
cat("Ferroptosis genes found in endothelial subset:", length(genes_use), "\n")

obj <- AddModuleScore(
  object = obj,
  features = list(genes_use),
  assay = "RNA",
  name = "Ferroptosis_Score",
  ctrl = 100,
  seed = 1234
)

obj$ferroptosis_score <- obj$Ferroptosis_Score1
obj$HA_vs_other <- ifelse(as.character(obj@meta.data[[group_col]]) == ha_label, "HA", "other")
obj$HA_vs_other <- factor(obj$HA_vs_other, levels = c("HA", "other"))
score_max <- max(obj$ferroptosis_score, na.rm = TRUE)

summary_global <- aggregate(
  ferroptosis_score ~ HA_vs_other,
  data = obj@meta.data,
  FUN = function(x) c(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    n = length(x)
  )
)

summary_global_out <- data.frame(
  group = summary_global$HA_vs_other,
  mean_score = summary_global$ferroptosis_score[, "mean"],
  median_score = summary_global$ferroptosis_score[, "median"],
  sd_score = summary_global$ferroptosis_score[, "sd"],
  n_cells = summary_global$ferroptosis_score[, "n"],
  stringsAsFactors = FALSE
)

write.csv(summary_global_out, file.path(res_dir, "ferroptosis_score_global_HA_vs_other.csv"), row.names = FALSE)

summary_subtype <- aggregate(
  ferroptosis_score ~ endothelial_subtype + HA_vs_other,
  data = obj@meta.data,
  FUN = function(x) c(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    n = length(x)
  )
)

summary_subtype_out <- data.frame(
  endothelial_subtype = summary_subtype$endothelial_subtype,
  group = summary_subtype$HA_vs_other,
  mean_score = summary_subtype$ferroptosis_score[, "mean"],
  median_score = summary_subtype$ferroptosis_score[, "median"],
  sd_score = summary_subtype$ferroptosis_score[, "sd"],
  n_cells = summary_subtype$ferroptosis_score[, "n"],
  stringsAsFactors = FALSE
)

write.csv(summary_subtype_out, file.path(res_dir, "ferroptosis_score_by_endothelial_subtype_HA_vs_other.csv"), row.names = FALSE)

if ("sample_id" %in% colnames(obj@meta.data)) {
  summary_sample <- aggregate(
    ferroptosis_score ~ sample_id + endothelial_subtype + HA_vs_other,
    data = obj@meta.data,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      n = length(x)
    )
  )

  summary_sample_out <- data.frame(
    sample_id = summary_sample$sample_id,
    endothelial_subtype = summary_sample$endothelial_subtype,
    group = summary_sample$HA_vs_other,
    mean_score = summary_sample$ferroptosis_score[, "mean"],
    median_score = summary_sample$ferroptosis_score[, "median"],
    n_cells = summary_sample$ferroptosis_score[, "n"],
    stringsAsFactors = FALSE
  )

  write.csv(summary_sample_out, file.path(res_dir, "ferroptosis_score_by_sample_endothelial_subtype_HA_vs_other.csv"), row.names = FALSE)
}

# ferroptosis score UMAP in endothelial subset
cat("\nPlotting ferroptosis score UMAP...\n")
p_score <- FeaturePlot(
  object = obj,
  features = "ferroptosis_score",
  reduction = reduction_name,
  raster = FALSE,
  order = TRUE,
  cols = c(ferroptosis_low_color, ferroptosis_high_color),
  min.cutoff = 0,
  max.cutoff = score_max
)

ggsave(
  filename = file.path(fig_dir, "umap_ferroptosis_score_endothelial_subtypes.png"),
  plot = p_score,
  width = 8,
  height = 7,
  dpi = 600
)

# split UMAP HA vs other
cat("Plotting split ferroptosis score UMAP...\n")
p_score_split <- FeaturePlot(
  object = obj,
  features = "ferroptosis_score",
  reduction = reduction_name,
  split.by = "HA_vs_other",
  raster = FALSE,
  order = TRUE,
  cols = c(ferroptosis_low_color, ferroptosis_high_color),
  min.cutoff = 0,
  max.cutoff = score_max
)

ggsave(
  filename = file.path(fig_dir, "umap_ferroptosis_score_endothelial_subtypes_split_HA_vs_other.png"),
  plot = p_score_split,
  width = 14,
  height = 7,
  dpi = 600
)

# violin HA vs other
cat("Plotting ferroptosis score violin plots...\n")
vln_df_group <- obj@meta.data[, c("HA_vs_other", "ferroptosis_score"), drop = FALSE]

p_vln_group <- ggplot(vln_df_group, aes(x = HA_vs_other, y = ferroptosis_score, fill = HA_vs_other)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("HA" = "#7A1F2B", "other" = "#D7B7BC")) +
  labs(x = NULL, y = "Ferroptosis_Score") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "none"
  )

ggsave(
  filename = file.path(fig_dir, "violin_ferroptosis_score_endothelial_subtypes_HA_vs_other.png"),
  plot = p_vln_group,
  width = 6,
  height = 6,
  dpi = 600
)

# violin by endothelial subtype split HA vs other
subtype_levels <- levels(obj@meta.data[[label_col]])
obj@meta.data[[label_col]] <- factor(obj@meta.data[[label_col]], levels = subtype_levels)

vln_df_subtype <- obj@meta.data[, c(label_col, "HA_vs_other", "ferroptosis_score"), drop = FALSE]

p_vln_subtype <- ggplot(vln_df_subtype, aes(x = endothelial_subtype, y = ferroptosis_score, fill = HA_vs_other)) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    position = position_dodge(width = 0.8),
    color = NA
  ) +
  scale_fill_manual(values = c("HA" = "#7A1F2B", "other" = "#D7B7BC")) +
  labs(x = NULL, y = "Ferroptosis_Score", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "top"
  )

ggsave(
  filename = file.path(fig_dir, "violin_ferroptosis_score_by_endothelial_subtype_HA_vs_other.png"),
  plot = p_vln_subtype,
  width = 14,
  height = 7,
  dpi = 600
)

# heatmap of ferroptosis genes by endothelial subtype and HA status
dot_features <- genes_use[genes_use %in% rownames(obj)]
if (length(dot_features) > 0) {
  heatmap_subtypes <- levels(obj@meta.data[[label_col]])

  if (length(heatmap_subtypes) > 0) {
    obj$subtype_group <- paste(obj@meta.data[[label_col]], obj$HA_vs_other, sep = " | ")

    group_levels <- as.vector(rbind(
      paste(heatmap_subtypes, "HA", sep = " | "),
      paste(heatmap_subtypes, "other", sep = " | ")
    ))
    obj$subtype_group <- factor(obj$subtype_group, levels = group_levels)

    avg_expr <- AverageExpression(
      object = obj,
      assays = "RNA",
      features = dot_features,
      group.by = "subtype_group",
      slot = "data",
      verbose = FALSE
    )$RNA

    avg_expr <- avg_expr[, group_levels[group_levels %in% colnames(avg_expr)], drop = FALSE]

    if (ncol(avg_expr) > 0) {
      scaled_expr <- t(scale(t(as.matrix(avg_expr))))
      scaled_expr[is.na(scaled_expr)] <- 0

      heatmap_df <- as.data.frame(as.table(scaled_expr), stringsAsFactors = FALSE)
      colnames(heatmap_df) <- c("gene", "group", "z_score")
      heatmap_df$gene <- factor(heatmap_df$gene, levels = rev(dot_features))
      heatmap_df$group <- factor(heatmap_df$group, levels = colnames(scaled_expr))
      heatmap_groups <- strsplit(as.character(heatmap_df$group), " | ", fixed = TRUE)
      heatmap_df$endothelial_subtype <- factor(vapply(heatmap_groups, `[`, character(1), 1), levels = heatmap_subtypes)
      heatmap_df$HA_vs_other <- factor(vapply(heatmap_groups, `[`, character(1), 2), levels = c("HA", "other"))

      p_heatmap <- ggplot(heatmap_df, aes(x = HA_vs_other, y = gene, fill = z_score)) +
        geom_tile(color = "white", linewidth = 0.25) +
        facet_grid(
          cols = vars(endothelial_subtype),
          scales = "free_x",
          space = "free_x",
          labeller = labeller(endothelial_subtype = wrap_each_word)
        ) +
        scale_fill_gradient2(
          low = "#4C78A8",
          mid = "white",
          high = "#7A1F2B",
          midpoint = 0,
          name = "Scaled\nexpression"
        ) +
        scale_x_discrete(position = "bottom") +
        scale_y_discrete(position = "right") +
        theme_classic(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          panel.spacing.x = grid::unit(0.22, "cm"),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 8.5, color = "black", lineheight = 0.95),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 11),
          legend.position = "top",
          legend.key.width = grid::unit(0.5, "cm"),
          legend.key.height = grid::unit(0.35, "cm"),
          legend.spacing.x = grid::unit(0.08, "cm"),
          plot.margin = margin(20, 35, 35, 25)
        )

      ggsave(
        filename = file.path(fig_dir, "heatmap_ferroptosis_genes_by_endothelial_subtype_HA_vs_other.png"),
        plot = p_heatmap,
        width = max(12, 1.15 * length(unique(heatmap_df$endothelial_subtype))),
        height = max(10, 0.16 * nrow(scaled_expr)),
        dpi = 600
      )
    }
  }
}

cat("\n============================================================\n")
cat("Endothelial ferroptosis analysis complete.\n")
cat("Figures dir  : ", fig_dir, "\n")
cat("Results dir  : ", res_dir, "\n")
cat("============================================================\n")
