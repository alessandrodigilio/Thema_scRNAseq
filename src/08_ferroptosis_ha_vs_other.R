######################################################
### Atlas iron and ferroptosis analysis: HA vs other ###
######################################################

# compute a ferroptosis score in the annotated atlas and summarize
# iron-related genes found in the atlas pseudobulk DESeq2 results

suppressPackageStartupMessages({
  library(Seurat)
  library(openxlsx)
  library(ggplot2)
})

# work from the project root
setwd("~/Thema_R")
source("src/global_config.R")
source("src/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "annotated.rds")
pseudobulk_res_dir <- file.path(results_dir, "pseudobulk_deseq2")
res_dir <- file.path(results_dir, "ferroptosis")
fig_dir <- file.path(figures_dir, "ferroptosis")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# parameters
geneset_file <- ferroptosis_geneset_file # in global_config.R
iron_gene_files <- c(
  ferroptosis_geneset_file,
  file.path(metadata_dir, "iron_genes", "iron_uptake_transport_genes.xlsx")
)

summary_file <- file.path(pseudobulk_res_dir, "deseq2_pseudobulk_summary_by_celltype.csv")
reduction_name <- "umap.harmony.rna"
group_col <- "condition"
ha_label <- "HA"
padj_thr <- 0.05
score_col <- "ferroptosis_score"
score_colors <- c("#F5E6E8", "#7A1F2B")
score_group_colors <- c("HA" = "#7A1F2B", "other" = "#D7B7BC")

# cell types to include in the heatmap (specific focus)
heatmap_celltypes <- c(
  "Sublining fibroblasts (SFRP2+)",
  "Lining fibroblasts (PRG4+)",
  "Inflammatory fibroblasts (ADAM12+)",
  "Remodeling fibroblasts (HTRA1+)",
  "Activated endothelial cells",
  "Inflammatory macrophages (IL1B+)",
  "Resident macrophages (C1QC+)",
  "Destructive lining fibroblasts (MMP3+)",
  "Pericytes / vascular smooth muscle cells"
)

# output files
global_score_file <- file.path(res_dir, "ferroptosis_score_global_HA_vs_other.csv")
celltype_score_file <- file.path(res_dir, "ferroptosis_score_by_celltype_HA_vs_other.csv")
sample_score_file <- file.path(res_dir, "ferroptosis_score_by_sample_celltype_HA_vs_other.csv")
iron_csv_file <- file.path(res_dir, "iron_related_genes_in_pseudobulk_DEGs_HA_vs_other.csv")
iron_xlsx_file <- file.path(res_dir, "iron_related_genes_in_pseudobulk_DEGs_HA_vs_other.xlsx")
iron_bubble_file <- file.path(fig_dir, "iron_related_genes_in_pseudobulk_DEGs_HA_vs_other_bubble_heatmap.png")

# load object and prepare RNA assay
cat("Loading annotated object...\n")
obj <- readRDS(input_object)
obj <- prepare_rna_assay_for_scoring(obj) # helper in utils.R

# load curated ferroptosis genes and compute the score
ferroptosis_genes <- load_excel_gene_list(geneset_file) # helper in utils.R
genes_use <- intersect(ferroptosis_genes, rownames(obj))

cat("Ferroptosis genes in metadata:", length(ferroptosis_genes), "\n")
cat("Ferroptosis genes found in object:", length(genes_use), "\n")

if (length(genes_use) == 0) stop("No ferroptosis genes were found in the annotated object.")

obj <- add_ferroptosis_score(obj, genes_use, score_col) # helper in utils.R
obj$HA_vs_other <- ifelse(as.character(obj@meta.data[[group_col]]) == ha_label, "HA", "other")
obj$HA_vs_other <- factor(obj$HA_vs_other, levels = c("HA", "other"))
score_max <- max(obj[[score_col]][, 1], na.rm = TRUE)

# summarize score globally, by cell type and by sample
global_summary <- summarize_score_by(obj@meta.data, score_col, "HA_vs_other") # helper in utils.R
celltype_summary <- summarize_score_by(obj@meta.data, score_col, c("cell_type", "HA_vs_other"))
sample_summary <- summarize_score_by(obj@meta.data, score_col, c("sample_id", "cell_type", "HA_vs_other"))
# save
write.csv(global_summary, global_score_file, row.names = FALSE)
write.csv(celltype_summary, celltype_score_file, row.names = FALSE)
write.csv(sample_summary, sample_score_file, row.names = FALSE)

# plot ferroptosis score on atlas UMAP
cat("Plotting ferroptosis score UMAPs...\n")
p_score <- make_score_umap(obj, score_col, reduction_name, score_colors, score_max) # helper in utils.R
# split by HA vs other to see the differenes
p_score_split <- make_score_umap(obj, score_col, reduction_name, score_colors, score_max, split_by = "HA_vs_other")
# save plots
ggsave(file.path(fig_dir, "umap_ferroptosis_score.png"), p_score, width = 8, height = 7, dpi = 600)
ggsave(file.path(fig_dir, "umap_ferroptosis_score_split_HA_vs_other.png"), p_score_split, width = 14, height = 7, dpi = 600)

# plot ferroptosis score distributions
cat("Plotting ferroptosis score distributions...\n")
celltype_levels <- names(cluster_name_colors)[names(cluster_name_colors) %in% unique(as.character(obj$cell_type))]
obj$cell_type <- factor(obj$cell_type, levels = celltype_levels)
p_vln_group <- make_score_violin_group(obj@meta.data, score_col, "HA_vs_other", score_group_colors, "Ferroptosis_Score") # helper in utils.R
p_vln_celltype <- make_score_violin_celltype(obj@meta.data, score_col, "cell_type", "HA_vs_other", score_group_colors, "Ferroptosis_Score")
# save plots
ggsave(file.path(fig_dir, "violin_ferroptosis_score_HA_vs_other.png"), p_vln_group, width = 6, height = 6, dpi = 600)
ggsave(file.path(fig_dir, "violin_ferroptosis_score_by_celltype_HA_vs_other.png"), p_vln_celltype, width = 14, height = 7, dpi = 600)

# plot average ferroptosis genes by cell type and condition
cat("Plotting ferroptosis gene heatmap...\n")
p_heatmap <- make_ferroptosis_gene_heatmap(obj, genes_use, heatmap_celltypes) # helper in utils.R

if (!is.null(p_heatmap)) {
  ggsave(
    file.path(fig_dir, "heatmap_ferroptosis_genes_by_celltype_HA_vs_other.png"),
    p_heatmap$plot,
    width = max(12, 1.15 * p_heatmap$n_celltypes),
    height = max(10, 0.16 * p_heatmap$n_genes),
    dpi = 600
  )
}

# find iron-related genes inside atlas pseudobulk DESeq2 results
cat("Scanning pseudobulk DEGs for iron-related genes...\n")
iron_genes <- load_excel_gene_list(iron_gene_files)
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
summary_df <- summary_df[summary_df$status == "tested", , drop = FALSE]

# iron genes found in the pseudobulk DEGs
iron_results <- scan_iron_genes_in_pseudobulk(summary_df, pseudobulk_res_dir, iron_genes, padj_thr) # helper in utils.R
# save
write.csv(iron_results$all_hits, iron_csv_file, row.names = FALSE)
save_iron_gene_workbook(iron_results, iron_xlsx_file) # helper in utils.R

# plot iron genes bubble heatmap
cat("Plotting iron-related gene bubble heatmap...\n")
p_iron_bubble <- make_iron_gene_bubble_plot(iron_results$significant_hits, summary_df$cell_type) # helper in utils.R

if (!is.null(p_iron_bubble)) {
  ggsave(
    iron_bubble_file,
    p_iron_bubble,
    width = 12,
    height = max(6, 0.28 * length(unique(iron_results$significant_hits$gene)) + 2),
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Atlas iron and ferroptosis analysis complete.\n")
cat("Ferroptosis genes found in object: ", length(genes_use), "\n", sep = "")
cat("Iron-related DEGs found: ", nrow(iron_results$all_hits), "\n", sep = "")
cat("Significant iron-related DEGs: ", nrow(iron_results$significant_hits), "\n", sep = "")
cat("Figures dir: ", fig_dir, "\n", sep = "")
cat("Results dir: ", res_dir, "\n", sep = "")
cat("============================================================\n")
