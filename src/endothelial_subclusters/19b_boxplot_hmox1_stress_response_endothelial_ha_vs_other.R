##################################################
### HMOX1 in stress-response endothelial cells ###
##################################################

# simple sample-level boxplot for one gene in the stress-response
# endothelial subtype comparing HA vs other.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

setwd("~/Thema_R")
source("src/global_config.R")

# input and output
input_object <- file.path(data_dir, "integrated_object", "annotated_endothelial_states.rds")
fig_dir <- file.path(figures_dir, "endothelial_gene_violin")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
target_gene <- "HMOX1"
group_col <- "condition"
sample_col <- "sample_id"
subtype_col <- "endothelial_subtype"
stress_response_subtype <- "Stress-response endothelial cells (HSPA6+)"
group_levels <- c("other", "HA")
group_colors <- c("other" = "#B65A5A", "HA" = "#5B8DB8")
box_alpha <- 0.9
output_file <- file.path(fig_dir, "boxplot_HMOX1_stress_response_endothelial_cells_HA_vs_other.png")

# load object
if (!file.exists(input_object)) {
  stop("Missing endothelial-annotated object: ", input_object)
}

cat("Loading endothelial-annotated object...\n")
obj <- readRDS(input_object)
cat("Cells:", ncol(obj), "\n")

# check metadata and gene
required_cols <- c(group_col, sample_col, subtype_col)
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

if (!target_gene %in% rownames(obj)) {
  stop("Gene not found in object: ", target_gene)
}

# keep only the target endothelial subtype
cells_use <- rownames(obj@meta.data)[obj@meta.data[[subtype_col]] == stress_response_subtype]
if (length(cells_use) == 0) {
  stop("No cells found for subtype: ", stress_response_subtype)
}

obj_sub <- subset(obj, cells = cells_use)
cat("Cells in target subtype:", ncol(obj_sub), "\n")

# prepare sample-level plotting data
plot_df <- FetchData(obj_sub, vars = c(target_gene, group_col, sample_col))
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
sample_df$condition <- factor(sample_df$condition, levels = group_levels)

cat("Samples by condition:\n")
print(table(sample_df$condition))
cat("Sample-level mean expression by condition:\n")
print(tapply(sample_df$expression, sample_df$condition, median, na.rm = TRUE))
cat("Sample-level values:\n")
print(sample_df[order(sample_df$condition, sample_df$sample_id), , drop = FALSE])

# build boxplot
p_box <- ggplot(sample_df, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(
    width = 0.44,
    alpha = box_alpha,
    color = "black",
    linewidth = 0.45,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = paste0(target_gene, " in ", stress_response_subtype),
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
  filename = output_file,
  plot = p_box,
  width = 8,
  height = 7,
  dpi = 600
)

cat("\n============================================================\n")
cat("Boxplot complete.\n")
cat("Subtype : ", stress_response_subtype, "\n")
cat("Gene    : ", target_gene, "\n")
cat("Output  : ", output_file, "\n")
cat("============================================================\n")
