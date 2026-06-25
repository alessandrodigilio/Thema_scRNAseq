#####################################################
### HMOX1 in destructive lining fibroblast states ###
#####################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# work from the project root
setwd("~/Thema_R")
source("src/global_config.R")

# input and output
input_subset_object <- file.path(data_dir, "integrated_object", "destructive_lining_fibroblasts_subclustered.rds")
fig_dir <- file.path(figures_dir, "destructive_lining_fibroblast_final_annotation")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

res_dir <- file.path(results_dir, "destructive_lining_fibroblast_final_annotation")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
gene_to_plot <- "HMOX1"
cluster_col <- "destructive_lining_fibroblast_subcluster"
label_col <- "destructive_lining_fibroblast_subtype"
reduction_name <- "umap.harmony.destructive.lining.fibroblast"
pt_size <- 1.05

destructive_lining_fibroblast_subcluster_labels <- c(
  "0" = "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)",
  "1" = "Activated MMP3+ lining fibroblast cells (ID1+)",
  "2" = "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)",
  "3" = "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)",
  "4" = "MMP3+ lining fibroblast cells (FAM184A+)",
  "5" = "HA-enriched SFRP2+ matrix fibroblast-like cells"
)

# load subset object
if (!file.exists(input_subset_object)) {
  stop("Missing destructive lining fibroblast subset object: ", input_subset_object)
}

cat("Loading destructive lining fibroblast subset object...\n")
obj <- readRDS(input_subset_object)
cat("Cells:", ncol(obj), "\n")

if (!cluster_col %in% colnames(obj@meta.data)) {
  stop("Missing cluster column: ", cluster_col)
}

if (!reduction_name %in% names(obj@reductions)) {
  stop("Missing reduction: ", reduction_name)
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

if (!gene_to_plot %in% rownames(obj)) {
  stop("Gene not found in object: ", gene_to_plot)
}

# apply the same final labels used in script 22
cluster_ids <- as.character(obj@meta.data[[cluster_col]])
mapped_labels <- unname(destructive_lining_fibroblast_subcluster_labels[cluster_ids])
mapped_labels[is.na(mapped_labels)] <- "Unknown"

obj@meta.data[[label_col]] <- factor(
  mapped_labels,
  levels = unname(destructive_lining_fibroblast_subcluster_labels)
)
obj@meta.data[[label_col]] <- droplevels(obj@meta.data[[label_col]])

cat("Subtypes:\n")
print(table(obj@meta.data[[label_col]], useNA = "ifany"))

# save a compact HMOX1 expression summary by subtype
plot_df <- FetchData(obj, vars = c(gene_to_plot, label_col))
colnames(plot_df) <- c("expression", "subtype")
summary_df <- aggregate(
  expression ~ subtype,
  data = plot_df,
  FUN = function(x) c(mean = mean(x), median = median(x), pct_expr = mean(x > 0) * 100)
)
summary_out <- data.frame(
  subtype = summary_df$subtype,
  mean_expression = summary_df$expression[, "mean"],
  median_expression = summary_df$expression[, "median"],
  pct_expressing = summary_df$expression[, "pct_expr"],
  stringsAsFactors = FALSE
)

write.csv(summary_out, file.path(res_dir, "HMOX1_expression_summary_destructive_lining_fibroblast_subtypes.csv"), row.names = FALSE)

cat("HMOX1 summary:\n")
print(summary_out)

# featurePlot for HMOX1
p_feature <- FeaturePlot(
  object = obj,
  features = gene_to_plot,
  reduction = reduction_name,
  raster = FALSE,
  order = TRUE,
  pt.size = pt_size,
  cols = c("#F3F1EC", "#7A1F2B")
) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 24) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 28, color = "black"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 19, color = "black"),
    plot.margin = margin(28, 28, 28, 28)
  )

ggsave(
  filename = file.path(fig_dir, "featureplot_HMOX1_destructive_lining_fibroblast_subtypes.png"),
  plot = p_feature,
  width = 10,
  height = 8.2,
  dpi = 600
)

# dotplot for HMOX1 across final subtypes
p_dot <- DotPlot(
  object = obj,
  features = gene_to_plot,
  group.by = label_col,
  cols = c("#E8E2DC", "#7A1F2B"),
  dot.scale = 8,
  col.min = 0,
  col.max = 3
) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.justification = c(0.65, 0.5),
    legend.key.width = grid::unit(0.45, "cm"),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.spacing.x = grid::unit(0.10, "cm"),
    plot.margin = margin(20, 35, 25, 45)
  )

ggsave(
  filename = file.path(fig_dir, "dotplot_HMOX1_destructive_lining_fibroblast_subtypes.png"),
  plot = p_dot,
  width = 9,
  height = 4.8,
  dpi = 600
)

cat("\n============================================================\n")
cat("HMOX1 destructive lining fibroblast plots complete.\n")
cat("FeaturePlot: ", file.path(fig_dir, "featureplot_HMOX1_destructive_lining_fibroblast_subtypes.png"), "\n")
cat("DotPlot    : ", file.path(fig_dir, "dotplot_HMOX1_destructive_lining_fibroblast_subtypes.png"), "\n")
cat("Summary    : ", file.path(res_dir, "HMOX1_expression_summary_destructive_lining_fibroblast_subtypes.csv"), "\n")
cat("============================================================\n")
