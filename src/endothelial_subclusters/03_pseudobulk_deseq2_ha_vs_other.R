#############################################################
### Endothelial subtype pseudobulk expression: HA vs other ###
#############################################################

# run sample-level pseudobulk DESeq2 for each endothelial subtype

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
  library(openxlsx)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")
source("src/endothelial_subclusters/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "annotated_endothelial_states.rds")
res_dir <- file.path(endothelial_results_dir, "pseudobulk_deseq2_endothelial_subtypes")
fig_dir <- file.path(endothelial_figures_dir, "pseudobulk_deseq2_endothelial_subtypes")
volcano_dir <- file.path(fig_dir, "volcano_plot")
gene_fig_dir <- file.path(endothelial_figures_dir, "endothelial_gene_violin")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_fig_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
condition_col <- "condition"
sample_col <- "sample_id"
subtype_col <- "endothelial_subtype"
group_levels <- c("other", "HA")
group_colors <- c("other" = "#B65A5A", "HA" = "#5B8DB8")

# parameters for pseudobulk DESeq2
min_cells_per_sample <- 10
min_pseudobulk_count <- 10
min_detected_cell_fraction <- 0.05
padj_thr <- 0.05
top_n_labels <- 15
log2fc_plot_limit <- 5
target_gene <- "HMOX1"
stress_response_subtype <- "Stress-response endothelial cells (HSPA6+)"
iron_gene_files <- iron_related_geneset_files

# load final endothelial-annotated object
cat("Loading endothelial-annotated object...\n")
obj <- readRDS(input_object)
cat("Cells:", ncol(obj), "\n")

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
rna_counts <- LayerData(obj[["RNA"]], layer = "counts")

# metadata used for pseudobulk
meta_cols <- c(condition_col, sample_col, subtype_col, "nCount_RNA", "nFeature_RNA", "percent.mt")
meta_df <- obj@meta.data[, meta_cols, drop = FALSE]
colnames(meta_df)[colnames(meta_df) == condition_col] <- "condition"
colnames(meta_df)[colnames(meta_df) == sample_col] <- "sample_id"
colnames(meta_df)[colnames(meta_df) == subtype_col] <- "cell_type"
meta_df$condition <- factor(ifelse(as.character(meta_df$condition) == "HA", "HA", "other"), levels = group_levels)
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$cell_type <- as.character(meta_df$cell_type)
meta_df <- meta_df[!is.na(meta_df$sample_id) & meta_df$sample_id != "" & !is.na(meta_df$cell_type) & meta_df$cell_type != "", , drop = FALSE]

endothelial_subtypes <- names(endothelial_subtype_colors)
endothelial_subtypes <- endothelial_subtypes[endothelial_subtypes %in% unique(meta_df$cell_type)]

# run pseudobulk DESeq2 subtype by subtype
cat("Running endothelial subtype pseudobulk DESeq2...\n")
result_files <- character()
summary_rows <- list()

for (subtype in endothelial_subtypes) {
  res <- run_pseudobulk_deseq2_celltype(
    cell_type = subtype,
    rna_counts = rna_counts,
    meta_df = meta_df,
    res_dir = res_dir,
    group_levels = group_levels,
    group_colors = group_colors,
    min_cells_per_sample = min_cells_per_sample,
    min_pseudobulk_count = min_pseudobulk_count,
    min_detected_cell_fraction = min_detected_cell_fraction,
    padj_thr = padj_thr
  )

  summary_rows[[subtype]] <- res$summary
  if (!is.na(res$result_file)) result_files <- c(result_files, res$result_file)
}

summary_df <- do.call(rbind, summary_rows)
summary_df$endothelial_subtype <- summary_df$cell_type
write.csv(summary_df, file.path(res_dir, "deseq2_pseudobulk_summary_by_endothelial_subtype.csv"), row.names = FALSE)

# volcano plots
save_pseudobulk_volcanoes(result_files, volcano_dir, group_colors, padj_thr, top_n_labels, log2fc_plot_limit)
save_cumulative_pseudobulk_volcano(result_files, summary_df, volcano_dir, endothelial_subtype_colors, padj_thr, top_n_labels, log2fc_plot_limit)
file.rename(
  file.path(volcano_dir, "volcano_cumulative_HA_vs_other.png"),
  file.path(volcano_dir, "volcano_cumulative_HA_vs_other_endothelial_subtypes.png")
)

# significant DEG summary
deg_summary_df <- summarize_pseudobulk_significant_genes(summary_df, res_dir)
colnames(deg_summary_df)[colnames(deg_summary_df) == "cell_type"] <- "endothelial_subtype"
write.csv(deg_summary_df, file.path(res_dir, "deseq2_significant_genes_summary_endothelial_subtypes.csv"), row.names = FALSE)

# iron-related DEG bubble heatmap
iron_genes <- load_excel_gene_list(iron_gene_files)
iron_results <- scan_iron_genes_in_pseudobulk(summary_df, res_dir, iron_genes, padj_thr)
p_iron_bubble <- make_iron_gene_bubble_plot(iron_results$significant_hits, summary_df$cell_type)

if (!is.null(p_iron_bubble)) {
  ggsave(
    file.path(fig_dir, "bubble_heatmap_iron_related_significant_genes_endothelial_subtypes.png"),
    p_iron_bubble,
    width = 12,
    height = max(6, 0.28 * length(unique(iron_results$significant_hits$gene)) + 2),
    dpi = 600
  )
}

# ------------------------------------------------------------------------ #
# --- HMOX1 expression in stress-response endothelial cells (HSPA6+) --- #
# ------------------------------------------------------------------------ #

# HMOX1 is a key gene in the stress-response/iron-related pathways

# plot HMOX1 in stress-response endothelial cells
cat("Plotting HMOX1 in stress-response endothelial cells...\n")
obj <- prepare_rna_assay_for_scoring(obj)
hmox1_plot <- plot_gene_by_sample_in_endothelial_subtype(
  obj = obj,
  gene = target_gene,
  endothelial_subtype = stress_response_subtype,
  group_col = condition_col,
  sample_col = sample_col,
  subtype_col = subtype_col,
  group_levels = group_levels,
  group_colors = group_colors
)
# save
write.csv(hmox1_plot$data, file.path(res_dir, "HMOX1_stress_response_endothelial_cells_sample_expression.csv"), row.names = FALSE)
ggsave(file.path(gene_fig_dir, "boxplot_HMOX1_stress_response_endothelial_cells_HA_vs_other.png"), hmox1_plot$plot, width = 8, height = 7, dpi = 600)

cat("\n============================================================\n")
cat("Endothelial pseudobulk analysis complete.\n")
cat("Results dir: ", res_dir, "\n", sep = "")
cat("Figures dir: ", fig_dir, "\n", sep = "")
cat("============================================================\n")
