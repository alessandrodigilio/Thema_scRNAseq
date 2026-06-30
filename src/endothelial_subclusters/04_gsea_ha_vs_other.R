#######################################
### Endothelial subtype ranked GSEA ###
#######################################

# run ranked fgsea on endothelial-subtype pseudobulk DESeq2 results

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")
source("src/endothelial_subclusters/utils.R")

# directories
pseudobulk_res_dir <- file.path(endothelial_results_dir, "pseudobulk_deseq2_endothelial_subtypes")
fgsea_res_dir <- file.path(endothelial_results_dir, "fgsea_pseudobulk_endothelial_subtypes_ha_vs_other")
fgsea_table_dir <- file.path(fgsea_res_dir, "significant_tables")
fgsea_iron_table_dir <- file.path(fgsea_res_dir, "iron_related_significant_tables")
fgsea_curve_dir <- file.path(endothelial_figures_dir, "fgsea_pseudobulk_endothelial_subtypes_ha_vs_other", "gsea_curves")

dir.create(fgsea_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_iron_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_curve_dir, recursive = TRUE, showWarnings = FALSE)

# load tested endothelial subtype pseudobulk results
summary_file <- file.path(pseudobulk_res_dir, "deseq2_pseudobulk_summary_by_endothelial_subtype.csv")
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
summary_df <- summary_df[summary_df$status == "tested", , drop = FALSE]
if (!"cell_type" %in% colnames(summary_df)) summary_df$cell_type <- summary_df$endothelial_subtype

# run fgsea
cat("Loading MSigDB pathways...\n")
pathway_data <- load_msigdb_pathways("Homo sapiens")

cat("Running endothelial subtype fgsea...\n")
run_fgsea_pseudobulk(
  summary_df = summary_df,
  pseudobulk_res_dir = pseudobulk_res_dir,
  fgsea_table_dir = fgsea_table_dir,
  fgsea_iron_table_dir = fgsea_iron_table_dir,
  fgsea_curve_dir = fgsea_curve_dir,
  pathways = pathway_data$pathways,
  pathway_info = pathway_data$pathway_info,
  iron_related_patterns = iron_related_patterns,
  padj_thr = 0.05,
  min_size = 5,
  max_size = 500,
  fgsea_seed = 1234
)

cat("\n============================================================\n")
cat("Endothelial subtype fgsea complete.\n")
cat("Results dir: ", fgsea_res_dir, "\n", sep = "")
cat("Curves dir : ", fgsea_curve_dir, "\n", sep = "")
cat("============================================================\n")
