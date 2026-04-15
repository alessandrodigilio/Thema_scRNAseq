# ===============================================================
#   20b_GEX_boxplot_HMOX1_stress_response_endothelial_HA_vs_other.R
# ===============================================================

# Simple sample-level boxplot for one gene in the stress-response
# endothelial subtype comparing HA vs other.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Input and output
INPUT_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_annotated_endothelial_states.rds")
FIG_DIR <- file.path(FIGURES_DIR, "endothelial_gene_violin")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Set parameters
GENE_TO_PLOT <- "HMOX1"
GROUP_COL <- "condition"
SAMPLE_COL <- "sample_id"
SUBTYPE_COL <- "endothelial_subtype"
TARGET_SUBTYPE <- "Stress-response endothelial cells (HSPA6+)"
GROUP_LEVELS <- c("other", "HA")
GROUP_COLORS <- c("other" = "#B65A5A", "HA" = "#5B8DB8")
BOX_ALPHA <- 0.9
OUTPUT_FILE <- file.path(FIG_DIR, "boxplot_HMOX1_stress_response_endothelial_cells_HA_vs_other.png")

# Load object
if (!file.exists(INPUT_OBJECT)) {
  stop("Missing endothelial-annotated object: ", INPUT_OBJECT)
}

cat("Loading endothelial-annotated object...\n")
obj <- readRDS(INPUT_OBJECT)
cat("Cells:", ncol(obj), "\n")

# Check metadata and gene
required_cols <- c(GROUP_COL, SAMPLE_COL, SUBTYPE_COL)
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) obj)

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

# Keep only the target endothelial subtype
cells_use <- rownames(obj@meta.data)[obj@meta.data[[SUBTYPE_COL]] == TARGET_SUBTYPE]
if (length(cells_use) == 0) {
  stop("No cells found for subtype: ", TARGET_SUBTYPE)
}

obj_sub <- subset(obj, cells = cells_use)
cat("Cells in target subtype:", ncol(obj_sub), "\n")

# Prepare sample-level plotting data
plot_df <- FetchData(obj_sub, vars = c(GENE_TO_PLOT, GROUP_COL, SAMPLE_COL))
colnames(plot_df) <- c("expression", "condition", "sample_id")
plot_df$condition <- ifelse(as.character(plot_df$condition) == "HA", "HA", "other")
plot_df <- plot_df[
  !is.na(plot_df$condition) &
    !is.na(plot_df$sample_id) &
    as.character(plot_df$sample_id) != "",
  ,
  drop = FALSE
]

sample_df <- aggregate(
  expression ~ sample_id + condition,
  data = plot_df,
  FUN = mean
)
sample_df$condition <- factor(sample_df$condition, levels = GROUP_LEVELS)

cat("Samples by condition:\n")
print(table(sample_df$condition))
cat("Sample-level mean expression by condition:\n")
print(tapply(sample_df$expression, sample_df$condition, median, na.rm = TRUE))
cat("Sample-level values:\n")
print(sample_df[order(sample_df$condition, sample_df$sample_id), , drop = FALSE])

# Build boxplot
p_box <- ggplot(sample_df, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(
    width = 0.44,
    alpha = BOX_ALPHA,
    color = "black",
    linewidth = 0.45,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(
    title = paste0(GENE_TO_PLOT, " in ", TARGET_SUBTYPE),
    x = NULL,
    y = "Mean norm. expression per sample"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 17, hjust = 0.5, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1.1, color = "black"),
    axis.ticks = element_line(linewidth = 1.1, color = "black"),
    axis.ticks.length = grid::unit(0.2, "cm"),
    plot.margin = margin(18, 18, 18, 18)
  )

print(p_box)
ggsave(
  filename = OUTPUT_FILE,
  plot = p_box,
  width = 8,
  height = 7,
  dpi = 600
)

cat("\n============================================================\n")
cat("Boxplot complete.\n")
cat("Subtype : ", TARGET_SUBTYPE, "\n")
cat("Gene    : ", GENE_TO_PLOT, "\n")
cat("Output  : ", OUTPUT_FILE, "\n")
cat("============================================================\n")
