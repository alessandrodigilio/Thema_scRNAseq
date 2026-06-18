# ---------------------------------------
# 22b_GEX_HMOX1_destructive_lining_fibroblast_plots.R
# HMOX1 plots in MMP3+ lining fibroblast subclusters
# ---------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Input and output
INPUT_SUBSET_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_destructive_lining_fibroblasts_subclustered.rds")
FIG_DIR <- file.path(FIGURES_DIR, "destructive_lining_fibroblast_final_annotation")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RES_DIR <- file.path(RESULTS_DIR, "destructive_lining_fibroblast_final_annotation")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# Set parameters
GENE_TO_PLOT <- "HMOX1"
CLUSTER_COL <- "destructive_lining_fibroblast_subcluster"
LABEL_COL <- "destructive_lining_fibroblast_subtype"
REDUCTION_NAME <- "umap.harmony.destructive.lining.fibroblast"
PT_SIZE <- 1.05

DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS <- c(
  "0" = "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)",
  "1" = "Activated MMP3+ lining fibroblast cells (ID1+)",
  "2" = "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)",
  "3" = "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)",
  "4" = "MMP3+ lining fibroblast cells (FAM184A+)",
  "5" = "HA-enriched SFRP2+ matrix fibroblast-like cells"
)

# Load subset object
if (!file.exists(INPUT_SUBSET_OBJECT)) {
  stop("Missing destructive lining fibroblast subset object: ", INPUT_SUBSET_OBJECT)
}

cat("Loading destructive lining fibroblast subset object...\n")
obj <- readRDS(INPUT_SUBSET_OBJECT)
cat("Cells:", ncol(obj), "\n")

if (!CLUSTER_COL %in% colnames(obj@meta.data)) {
  stop("Missing cluster column: ", CLUSTER_COL)
}

if (!REDUCTION_NAME %in% names(obj@reductions)) {
  stop("Missing reduction: ", REDUCTION_NAME)
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

if (!GENE_TO_PLOT %in% rownames(obj)) {
  stop("Gene not found in object: ", GENE_TO_PLOT)
}

# Apply the same final labels used in script 22
cluster_ids <- as.character(obj@meta.data[[CLUSTER_COL]])
mapped_labels <- unname(DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS[cluster_ids])
mapped_labels[is.na(mapped_labels)] <- "Unknown"

obj@meta.data[[LABEL_COL]] <- factor(
  mapped_labels,
  levels = unname(DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS)
)
obj@meta.data[[LABEL_COL]] <- droplevels(obj@meta.data[[LABEL_COL]])

cat("Subtypes:\n")
print(table(obj@meta.data[[LABEL_COL]], useNA = "ifany"))

# Save a compact HMOX1 expression summary by subtype
plot_df <- FetchData(obj, vars = c(GENE_TO_PLOT, LABEL_COL))
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

write.csv(
  summary_out,
  file.path(RES_DIR, "HMOX1_expression_summary_destructive_lining_fibroblast_subtypes.csv"),
  row.names = FALSE
)

cat("HMOX1 summary:\n")
print(summary_out)

# FeaturePlot for HMOX1
p_feature <- FeaturePlot(
  object = obj,
  features = GENE_TO_PLOT,
  reduction = REDUCTION_NAME,
  raster = FALSE,
  order = TRUE,
  pt.size = PT_SIZE,
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
  filename = file.path(FIG_DIR, "featureplot_HMOX1_destructive_lining_fibroblast_subtypes.png"),
  plot = p_feature,
  width = 10,
  height = 8.2,
  dpi = 600
)

# DotPlot for HMOX1 across final subtypes
p_dot <- DotPlot(
  object = obj,
  features = GENE_TO_PLOT,
  group.by = LABEL_COL,
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
  filename = file.path(FIG_DIR, "dotplot_HMOX1_destructive_lining_fibroblast_subtypes.png"),
  plot = p_dot,
  width = 9,
  height = 4.8,
  dpi = 600
)

cat("\n============================================================\n")
cat("HMOX1 destructive lining fibroblast plots complete.\n")
cat("FeaturePlot: ", file.path(FIG_DIR, "featureplot_HMOX1_destructive_lining_fibroblast_subtypes.png"), "\n")
cat("DotPlot    : ", file.path(FIG_DIR, "dotplot_HMOX1_destructive_lining_fibroblast_subtypes.png"), "\n")
cat("Summary    : ", file.path(RES_DIR, "HMOX1_expression_summary_destructive_lining_fibroblast_subtypes.csv"), "\n")
cat("============================================================\n")
