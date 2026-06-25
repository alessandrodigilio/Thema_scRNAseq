#############################################################
### Atlas pseudobulk differential expression: HA vs other ###
#############################################################

# run sample-level pseudobulk DESeq2 for each annotated atlas cell type

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "annotated.rds")
res_dir <- file.path(results_dir, "pseudobulk_deseq2")
fig_dir <- file.path(figures_dir, "pseudobulk_deseq2")
volcano_dir <- file.path(fig_dir, "volcano_plot")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
condition_col <- "condition"
sample_col <- "sample_id"
celltype_col <- "cell_type"
group_levels <- c("other", "HA")
group_colors <- c("other" = "#B65A5A", "HA" = "#5B8DB8")

min_cells_per_sample <- 20
min_pseudobulk_count <- 10
min_detected_cell_fraction <- 0.05
padj_thr <- 0.05
top_n_labels <- 15
log2fc_plot_limit <- 5 # avoid extreme log2FC values for the volcanoplot
remove_mt_genes <- FALSE

# load final annotated object
cat("Loading final annotated object...\n")
obj <- readRDS(input_object)
cat("Cells:", ncol(obj), "\n")

DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
rna_counts <- LayerData(obj[["RNA"]], layer = "counts")

# keep metadata needed for sample-level pseudobulk
meta_df <- obj@meta.data[, c(condition_col, sample_col, celltype_col, "nCount_RNA", "nFeature_RNA", "percent.mt"), drop = FALSE]
colnames(meta_df) <- c("condition", "sample_id", "cell_type", "nCount_RNA", "nFeature_RNA", "percent.mt")

meta_df$condition <- ifelse(as.character(meta_df$condition) == "HA", "HA", "other")
meta_df$condition <- factor(meta_df$condition, levels = group_levels)
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$cell_type <- as.character(meta_df$cell_type)

cell_types <- names(cluster_name_colors)
cell_types <- cell_types[cell_types %in% unique(meta_df$cell_type)]

# run DESeq2 separately for each atlas cell type
summary_rows <- list()
result_files <- character()

# loop over cell types to run DESeq2 pseudobulk
for (cell_type in cell_types) {
  cat("\nCell type:", cell_type, "\n")
  
  # DESeq2 for the current cell type
  out <- run_pseudobulk_deseq2_celltype( # helper in utils.R
    cell_type = cell_type,
    rna_counts = rna_counts,
    meta_df = meta_df,
    res_dir = res_dir,
    group_levels = group_levels,
    group_colors = group_colors,
    min_cells_per_sample = min_cells_per_sample,
    min_pseudobulk_count = min_pseudobulk_count,
    min_detected_cell_fraction = min_detected_cell_fraction,
    padj_thr = padj_thr,
    remove_mt_genes = remove_mt_genes
  )

  summary_rows[[cell_type]] <- out$summary
  if (!is.na(out$result_file)) result_files <- c(result_files, out$result_file)
}

summary_df <- do.call(rbind, summary_rows)
# save
write.csv(summary_df, file.path(res_dir, "deseq2_pseudobulk_summary_by_celltype.csv"), row.names = FALSE)

# volcano plots from saved DESeq2 result tables
save_pseudobulk_volcanoes( # helper utils.R
  result_files = result_files,
  volcano_dir = volcano_dir,
  group_colors = group_colors,
  padj_thr = padj_thr,
  top_n_labels = top_n_labels,
  log2fc_plot_limit = log2fc_plot_limit
)

# cumulative volcano plot to see all cell types in the same plot
save_cumulative_pseudobulk_volcano( # helper utils.R
  result_files = result_files,
  summary_df = summary_df,
  volcano_dir = volcano_dir,
  cluster_name_colors = cluster_name_colors,
  padj_thr = padj_thr,
  top_n_labels = top_n_labels,
  log2fc_plot_limit = log2fc_plot_limit
)

# summarize significant genes per cell type
deg_summary_df <- summarize_pseudobulk_significant_genes(summary_df, res_dir)
write.csv(deg_summary_df, file.path(res_dir, "deseq2_significant_genes_summary.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("Pseudobulk DESeq2 analysis complete.\n")
cat("Results dir : ", res_dir, "\n", sep = "")
cat("Figures dir : ", fig_dir, "\n", sep = "")
cat("Summary     : ", file.path(res_dir, "deseq2_pseudobulk_summary_by_celltype.csv"), "\n", sep = "")
cat("DEG summary : ", file.path(res_dir, "deseq2_significant_genes_summary.csv"), "\n", sep = "")
cat("============================================================\n")
