#########################################################################
### Functional analysis of HA-enriched destructive lining fibroblasts ###
#########################################################################

# run functional analysis on genes differentially expressed in the two
# HA-enriched destructive/MMP3+ lining fibroblast states

suppressPackageStartupMessages({
  library(Seurat)
  library(enrichR)
  library(fgsea)
  library(ggplot2)
  library(patchwork)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/fibroblast_subcluster/utils.R")

# directories
de_res_dir <- file.path(fibroblast_results_dir, "destructive_lining_fibroblast_HA_enriched_state_DE")
de_file <- file.path(de_res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states.csv")
iron_de_file <- file.path(de_res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_iron_related.csv")
input_fibroblast_object <- file.path(data_dir, "integrated_object", "destructive_lining_fibroblasts_subclustered.rds")
func_res_dir <- file.path(fibroblast_results_dir, "functional_analysis_HA_enriched_destructive_lining_fibroblast_states")
fgsea_table_dir <- file.path(func_res_dir, "fgsea_significant_tables")
fgsea_iron_table_dir <- file.path(func_res_dir, "fgsea_iron_related_significant_tables")
enrichr_table_dir <- file.path(func_res_dir, "enrichr_ora_tables")
enrichr_selected_table_dir <- file.path(func_res_dir, "enrichr_iron_related_gene_tables")
func_fig_dir <- file.path(fibroblast_figures_dir, "functional_analysis_HA_enriched_destructive_lining_fibroblast_states")
bubble_fig_dir <- file.path(func_fig_dir, "bubble_heatmap")
fgsea_curve_dir <- file.path(func_fig_dir, "gsea_curves")
enrichr_fig_dir <- file.path(func_fig_dir, "enrichr_ora")
tgf_beta_fig_dir <- file.path(func_fig_dir, "TGFbeta_featureplots")

# create dir
dir.create(func_res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_iron_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichr_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichr_selected_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bubble_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_curve_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichr_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tgf_beta_fig_dir, recursive = TRUE, showWarnings = FALSE)

# settings
padj_thr <- 0.05
rank_method <- "avg_log2FC"
ha_enriched_group <- "HA-enriched DLF states"
other_state_group <- "Other DLF states"
iron_pathway_patterns <- c(iron_related_patterns, "oxidative stress")

# focus on selected pathways (inflammation, coagulation, wound healing, cilium, Golgi vesicle transport)
selected_gsea_curves <- c(
  "GOBP_INFLAMMATORY_RESPONSE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "GOBP_HEMOSTASIS",
  "GOBP_REGULATION_OF_COAGULATION",
  "REACTOME_DISSOLUTION_OF_FIBRIN_CLOT",
  "GOBP_RESPONSE_TO_WOUNDING",
  "GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT",
  "GOBP_CILIUM_ORGANIZATION",
  "GOBP_CELL_PROJECTION_ASSEMBLY",
  "GOBP_GOLGI_VESICLE_TRANSPORT"
)

# ontology databases for Enrichr ORA
enrichr_dbs <- c(
  "GO_Biological_Process_2025",
  "Reactome_2022",
  "KEGG_2021_Human"
)
enrichr_db_labels <- c(
  "GO_Biological_Process_2025" = "GO BP 2025",
  "Reactome_2022" = "Reactome 2022",
  "KEGG_2021_Human" = "KEGG 2021"
)

# load DE results
cat("Loading targeted FindMarkers table:\n")
cat(de_file, "\n")
de_df <- read.csv(de_file, stringsAsFactors = FALSE)

# run ranked GSEA on the full FindMarkers table
fgsea_out <- run_fibroblast_fgsea_analysis(
  de_df = de_df,
  func_res_dir = func_res_dir,
  fgsea_table_dir = fgsea_table_dir,
  fgsea_iron_table_dir = fgsea_iron_table_dir,
  fgsea_curve_dir = fgsea_curve_dir,
  selected_gsea_curves = selected_gsea_curves,
  iron_patterns = iron_pathway_patterns,
  rank_method = rank_method,
  padj_thr = padj_thr,
  min_size = 5,
  max_size = 500,
  fgsea_seed = 1234,
  use_multilevel = TRUE,
  n_permutations = 10000,
  msigdb_species = "Homo sapiens"
)

# plot significant iron-related DE genes from script 02
plot_fibroblast_iron_bubble_heatmap(
  iron_de_file = iron_de_file,
  output_table = file.path(func_res_dir, "significant_iron_related_genes_bubble_heatmap_input_HA_enriched_DLF_states_vs_other_states.csv"),
  output_file = file.path(bubble_fig_dir, "bubble_heatmap_significant_iron_related_genes_HA_enriched_DLF_states_vs_other_states.png"),
  ha_enriched_group = ha_enriched_group,
  other_state_group = other_state_group
)

# run Enrichr ORA on significant genes and keep iron-related terms
run_fibroblast_enrichr_ora(
  de_df = de_df,
  enrichr_table_dir = enrichr_table_dir,
  enrichr_selected_table_dir = enrichr_selected_table_dir,
  enrichr_fig_dir = enrichr_fig_dir,
  enrichr_dbs = enrichr_dbs,
  enrichr_db_labels = enrichr_db_labels,
  iron_patterns = iron_pathway_patterns,
  padj_thr = padj_thr,
  top_n_terms = 10,
  max_tries = 4,
  retry_wait_sec = 20
)

# plot TGF-beta-related genes on the destructive lining fibroblast UMAP
cat("\nPlotting TGF-beta-related genes in destructive lining fibroblasts...\n")
plot_tgf_beta_featureplots(
  input_object = input_fibroblast_object,
  output_file = file.path(tgf_beta_fig_dir, "featureplot_TGFbeta_destructive_lining_fibroblasts.png"),
  reduction_name = "umap.harmony.destructive.lining.fibroblast",
  features = tgf_beta_fibroblast_features
)

# save compact summary table
summary_df <- data.frame(
  comparison = "HA-enriched DLF states vs other DLF states",
  rank_method = rank_method,
  n_ranked_genes = length(fgsea_out$ranks),
  n_significant_pathways = nrow(fgsea_out$sig_df),
  n_iron_related_significant_pathways = nrow(fgsea_out$iron_sig_df),
  stringsAsFactors = FALSE
)
# save
write.csv(summary_df, file.path(func_res_dir, "fgsea_summary_HA_enriched_DLF_states_vs_other_DLF_states.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("Functional analysis for HA-enriched destructive lining fibroblast states complete.\n")
cat("All results             : ", fgsea_out$all_file, "\n", sep = "")
cat("Significant table       : ", fgsea_out$sig_file, "\n", sep = "")
cat("Iron-related significant: ", fgsea_out$iron_file, "\n", sep = "")
cat("Bubble heatmap figures  : ", bubble_fig_dir, "\n", sep = "")
cat("GSEA curves             : ", fgsea_curve_dir, "\n", sep = "")
cat("Enrichr ORA panels      : ", enrichr_fig_dir, "\n", sep = "")
cat("TGF-beta FeaturePlots   : ", tgf_beta_fig_dir, "\n", sep = "")
cat("Interpretation          : positive NES = enriched toward HA-enriched DLF states; negative NES = enriched toward other DLF states.\n")
cat("============================================================\n")
