###########################################################
### Macrophage subtype ferroptosis scoring: HA vs other ###
###########################################################

# compute a curated ferroptosis score inside the macrophage subclustered object

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(openxlsx)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")
source("src/macrophage_subclusters/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "macrophages_subclustered.rds")
fig_dir <- file.path(figures_dir, "macrophage_ferroptosis")
res_dir <- file.path(results_dir, "macrophage_ferroptosis")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
geneset_file <- ferroptosis_geneset_file
cluster_col <- "macrophage_subcluster"
label_col <- "macrophage_subtype"
reduction_name <- "umap.harmony.macrophage"
group_col <- "condition"
ha_label <- "HA"
score_col <- "ferroptosis_score"
score_colors <- c("#F5E6E8", "#7A1F2B")
score_group_colors <- c("HA" = "#7A1F2B", "other" = "#D7B7BC")

# load macrophage object and add final subtype labels
cat("Loading macrophage subset object...\n")
obj <- readRDS(input_object)
obj <- add_macrophage_subtype_labels(obj, cluster_col, label_col)
obj <- prepare_rna_assay_for_scoring(obj)

# load curated ferroptosis genes and compute the score
ferroptosis_genes <- load_excel_gene_list(geneset_file, include_column_names = FALSE)
genes_use <- intersect(ferroptosis_genes, rownames(obj))

# check
cat("Ferroptosis genes in metadata:", length(ferroptosis_genes), "\n")
cat("Ferroptosis genes found in macrophage subset:", length(genes_use), "\n")
if (length(genes_use) == 0) stop("No ferroptosis genes were found in the macrophage subset.")

# compute ferroptosis score and add HA vs other condition label
obj <- add_ferroptosis_score(obj, genes_use, score_col)
obj$HA_vs_other <- ifelse(as.character(obj@meta.data[[group_col]]) == ha_label, "HA", "other")
obj$HA_vs_other <- factor(obj$HA_vs_other, levels = c("HA", "other"))
score_max <- max(obj[[score_col]][, 1], na.rm = TRUE)

# summarize score globally, by subtype and by sample
global_summary <- summarize_score_by(obj@meta.data, score_col, "HA_vs_other")
subtype_summary <- summarize_score_by(obj@meta.data, score_col, c(label_col, "HA_vs_other"))
sample_summary <- summarize_score_by(obj@meta.data, score_col, c("sample_id", label_col, "HA_vs_other"))
# save
write.csv(global_summary, file.path(res_dir, "ferroptosis_score_global_HA_vs_other.csv"), row.names = FALSE)
write.csv(subtype_summary, file.path(res_dir, "ferroptosis_score_by_macrophage_subtype_HA_vs_other.csv"), row.names = FALSE)
write.csv(sample_summary, file.path(res_dir, "ferroptosis_score_by_sample_macrophage_subtype_HA_vs_other.csv"), row.names = FALSE)

# plot ferroptosis score on macrophage UMAP
cat("Plotting ferroptosis score UMAPs...\n")
p_score <- make_score_umap(obj, score_col, reduction_name, score_colors, score_max)
p_score_split <- make_score_umap(obj, score_col, reduction_name, score_colors, score_max, split_by = "HA_vs_other")
# save plots
ggsave(file.path(fig_dir, "umap_ferroptosis_score_macrophage_subtypes.png"), p_score, width = 8, height = 7, dpi = 600)
ggsave(file.path(fig_dir, "umap_ferroptosis_score_macrophage_subtypes_split_HA_vs_other.png"), p_score_split, width = 14, height = 7, dpi = 600)

# plot score distributions
cat("Plotting ferroptosis score distributions...\n")
p_vln_group <- make_score_violin_group(obj@meta.data, score_col, "HA_vs_other", score_group_colors, "Ferroptosis_Score")
p_vln_subtype <- make_score_violin_celltype(obj@meta.data, score_col, label_col, "HA_vs_other", score_group_colors, "Ferroptosis_Score")

ggsave(file.path(fig_dir, "violin_ferroptosis_score_macrophage_subtypes_HA_vs_other.png"), p_vln_group, width = 6, height = 6, dpi = 600)
ggsave(file.path(fig_dir, "violin_ferroptosis_score_by_macrophage_subtype_HA_vs_other.png"), p_vln_subtype, width = 14, height = 7, dpi = 600)

# plot average ferroptosis genes by macrophage subtype and condition
cat("Plotting ferroptosis gene heatmap...\n")
p_heatmap <- make_macrophage_ferroptosis_heatmap(obj, genes_use, label_col)

ggsave(
  file.path(fig_dir, "heatmap_ferroptosis_genes_by_macrophage_subtype_HA_vs_other.png"),
  p_heatmap$plot,
  width = max(12, 1.15 * p_heatmap$n_subtypes),
  height = max(10, 0.16 * p_heatmap$n_genes),
  dpi = 600
)

cat("\n============================================================\n")
cat("Macrophage ferroptosis analysis complete.\n")
cat("Figures dir: ", fig_dir, "\n", sep = "")
cat("Results dir: ", res_dir, "\n", sep = "")
cat("============================================================\n")
