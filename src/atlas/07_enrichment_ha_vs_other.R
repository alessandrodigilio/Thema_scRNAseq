######################################
### Atlas enrichment: HA vs other ###
######################################

# run ranked fgsea and Enrichr over-representation analysis on atlas pseudobulk DESeq2 results

suppressPackageStartupMessages({
  library(fgsea)
  library(enrichR)
  library(ggplot2)
  library(patchwork)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")

# directories
pseudobulk_res_dir <- file.path(results_dir, "pseudobulk_deseq2")
res_dir <- file.path(results_dir, "enrichment_ha_vs_other")
fig_dir <- file.path(figures_dir, "enrichment_ha_vs_other")
fgsea_table_dir <- file.path(res_dir, "fgsea", "significant_tables")
fgsea_iron_table_dir <- file.path(res_dir, "fgsea", "iron_related_significant_tables")
fgsea_curve_dir <- file.path(fig_dir, "fgsea", "gsea_curves")
enrichr_table_dir <- file.path(res_dir, "enrichr", "tables")
enrichr_selected_table_dir <- file.path(res_dir, "enrichr", "selected_pathway_gene_tables")
enrichr_selected_fig_dir <- file.path(fig_dir, "enrichr", "selected_pathway_gene_panels")

for (dir_here in c( # loop for dir creation
  fgsea_table_dir,
  fgsea_iron_table_dir,
  fgsea_curve_dir,
  enrichr_table_dir,
  enrichr_selected_table_dir,
  enrichr_selected_fig_dir
)) {
  dir.create(dir_here, recursive = TRUE, showWarnings = FALSE)
}

# set fgsea parameters
padj_thr <- 0.05
min_size <- 5
max_size <- 500
fgsea_seed <- 1234
msigdb_species <- "Homo sapiens"

# set Enrichr parameters
min_degs_per_celltype <- 30 # we consider an enrichment significant when we have more than 30 degs
selected_celltypes <- c(
  # focus on fibroblasts
  "Lining fibroblasts (PRG4+)",
  "Destructive lining fibroblasts (MMP3+)"
)

# use GOBP, Reactome and KEGG onthologies
dbs <- c(
  "GO_Biological_Process_2025",
  "Reactome_2022",
  "KEGG_2021_Human"
)

db_labels <- c(
  "GO_Biological_Process_2025" = "GO BP 2025",
  "Reactome_2022" = "Reactome 2022",
  "KEGG_2021_Human" = "KEGG 2021"
)

# load tested pseudobulk cell types
summary_file <- file.path(pseudobulk_res_dir, "deseq2_pseudobulk_summary_by_celltype.csv")
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
summary_df <- summary_df[summary_df$status == "tested", , drop = FALSE]

# run ranked fgsea using all genes from each DESeq2 table
cat("Loading MSigDB pathways...\n")
pathway_data <- load_msigdb_pathways(msigdb_species)

cat("Running fgsea...\n")
run_fgsea_pseudobulk(
  summary_df = summary_df,
  pseudobulk_res_dir = pseudobulk_res_dir,
  fgsea_table_dir = fgsea_table_dir,
  fgsea_iron_table_dir = fgsea_iron_table_dir,
  fgsea_curve_dir = fgsea_curve_dir,
  pathways = pathway_data$pathways,
  pathway_info = pathway_data$pathway_info,
  iron_related_patterns = iron_related_patterns, # in global_config.R
  padj_thr = padj_thr,
  min_size = min_size,
  max_size = max_size,
  fgsea_seed = fgsea_seed
)

# run Enrichr using significant genes from each DESeq2 table
cat("Running Enrichr...\n")
run_enrichr_pseudobulk(
  summary_df = summary_df,
  pseudobulk_res_dir = pseudobulk_res_dir,
  enrichr_table_dir = enrichr_table_dir,
  selected_table_dir = enrichr_selected_table_dir,
  selected_fig_dir = enrichr_selected_fig_dir,
  dbs = dbs,
  db_labels = db_labels,
  min_degs_per_celltype = min_degs_per_celltype,
  padj_cutoff = padj_thr,
  selected_celltypes = selected_celltypes,
  iron_related_patterns = iron_related_patterns, # in global_config.R
  tgf_beta_patterns = tgf_beta_patterns # in global_config.R
)

cat("\n============================================================\n")
cat("Atlas enrichment analysis complete.\n")
cat("Results dir: ", res_dir, "\n", sep = "")
cat("Figures dir: ", fig_dir, "\n", sep = "")
cat("============================================================\n")
