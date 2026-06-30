##############################################################
### DE of HA-enriched destructive lining fibroblast states ###
##############################################################

# compare the two HA-enriched destructive/MMP3+ lining fibroblast
# subclusters against the remaining destructive lining fibroblast states
# with cell-level DE

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/fibroblast_subcluster/utils.R")

# directories
input_subset_object <- file.path(data_dir, "integrated_object", "destructive_lining_fibroblasts_subclustered.rds")
res_dir <- file.path(results_dir, "destructive_lining_fibroblast_HA_enriched_state_DE")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(figures_dir, "destructive_lining_fibroblast_HA_enriched_state_DE")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
cluster_col <- "destructive_lining_fibroblast_subcluster"
label_col <- "destructive_lining_fibroblast_subtype"
state_group_col <- "destructive_lining_fibroblast_state_group"
condition_col <- "condition"
sample_col <- "sample_id"

ha_enriched_clusters <- c("2", "5") # cluster labels for the two HA-enriched destructive lining fibroblast subclusters
other_state_clusters <- c("0", "1", "3", "4") # cluster labels for the remaining destructive lining fibroblast subclusters
ha_enriched_group <- "HA-enriched DLF states" 
other_state_group <- "Other DLF states"
group_levels <- c(other_state_group, ha_enriched_group)

padj_thr <- 0.05
group_colors <- c(
  "Other DLF states" = "#8FA0A8",
  "HA-enriched DLF states" = "#C65A5A"
)

# genes related to iron metabolism, WNT signaling, extracellular matrix, chemokines, and inflammatory mediators
selected_genes <- c(
  "HMOX1", "NQO1", "FTL", "FTH1", "CP", "SLC40A1", "TFRC",
  "SFRP2", "SFRP1", "COMP", "PODN", "IGF1", "MFAP5",
  "CCL7", "CXCL1", "CCL20", "CCRL2", "BIRC3",
  "MMP3", "MMP1", "IL6", "PTGS2"
)

# load destructive lining fibroblast subset object
if (!file.exists(input_subset_object)) {
  stop("Missing destructive lining fibroblast subset object: ", input_subset_object)
}
cat("Loading destructive lining fibroblast subset object...\n")
obj <- readRDS(input_subset_object)
cat("Cells:", ncol(obj), "\n")

# check that metadata columns are present
required_cols <- c(cluster_col, condition_col, sample_col)
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
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

# apply final labels and define the targeted state groups
cluster_ids <- as.character(obj@meta.data[[cluster_col]])
obj <- add_destructive_lining_fibroblast_subtype_labels(obj, cluster_col = cluster_col, label_col = label_col)

state_group <- rep(NA_character_, length(cluster_ids))
state_group[cluster_ids %in% ha_enriched_clusters] <- ha_enriched_group
state_group[cluster_ids %in% other_state_clusters] <- other_state_group

obj@meta.data[[state_group_col]] <- factor(state_group, levels = group_levels)
obj <- subset(obj, cells = rownames(obj@meta.data)[!is.na(obj@meta.data[[state_group_col]])])
obj@meta.data[[state_group_col]] <- droplevels(obj@meta.data[[state_group_col]])

# check 
cat("Cells by final subtype:\n")
print(table(obj@meta.data[[label_col]], useNA = "ifany"))
cat("Cells by targeted state group:\n")
print(table(obj@meta.data[[state_group_col]], useNA = "ifany"))
cat("Cells by sample and targeted state group:\n")
print(table(obj@meta.data[[sample_col]], obj@meta.data[[state_group_col]]))

# ---------------------------------------- #
# -------- Cell-level DE analysis -------- #
# ---------------------------------------- #

# In this case the DE is cell type due to the limited number of samples and the fact that 
# the two HA-enriched destructive lining fibroblast subclusters are only present in HA samples
# --> DE between the two broad state groups (HA-enriched DLF states vs other DLF states) 

# run cell-level differential expression
cat("\nRunning FindMarkers: HA-enriched DLF states vs other DLF states...\n")
Idents(obj) <- obj@meta.data[[state_group_col]]
de_df <- FindMarkers(
  object = obj,
  assay = "RNA",
  ident.1 = ha_enriched_group,
  ident.2 = other_state_group,
  test.use = "wilcox",
  min.pct = 0.10, # require at least 10% of cells in either group to express the gene
  logfc.threshold = 0.10, 
  only.pos = FALSE,
  verbose = FALSE
)

# DE results table
de_df$gene <- rownames(de_df)
de_df <- de_df[, c("gene", setdiff(colnames(de_df), "gene"))]
de_df <- de_df[order(de_df$p_val_adj, de_df$p_val), , drop = FALSE]
de_df$direction <- ifelse(de_df$avg_log2FC > 0, "Higher in HA-enriched DLF states", "Higher in other DLF states")
de_df$is_significant <- !is.na(de_df$p_val_adj) & de_df$p_val_adj < padj_thr
# save
write.csv(de_df, file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states.csv"), row.names = FALSE)
sig_df <- de_df[de_df$is_significant, , drop = FALSE]
write.csv(sig_df, file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_significant.csv"), row.names = FALSE)

# subset DE results for iron-related genes
iron_genes <- load_gene_set_files(iron_related_geneset_files)
iron_de_df <- de_df[toupper(de_df$gene) %in% iron_genes, , drop = FALSE]
write.csv(iron_de_df, file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_iron_related.csv"), row.names = FALSE)
selected_de_df <- de_df[de_df$gene %in% selected_genes, , drop = FALSE]
write.csv(selected_de_df, file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_selected_genes.csv"), row.names = FALSE)

# SUMMARY
cat("Top DE genes:\n")
print(head(de_df, 20))
cat("Significant genes at adjusted p < ", padj_thr, ": ", nrow(sig_df), "\n", sep = "")
cat("Iron-related genes found in DE table:", nrow(iron_de_df), "\n")
cat("Significant iron-related genes:\n")
print(iron_de_df[iron_de_df$is_significant, c("gene", "avg_log2FC", "p_val_adj", "direction"), drop = FALSE])

# dotplot of selected genes across the two broad state groups
selected_genes_use <- selected_genes[selected_genes %in% rownames(obj)]

if (length(selected_genes_use) > 0) {
  p_dot <- DotPlot(
    object = obj,
    features = selected_genes_use,
    group.by = state_group_col,
    cols = c("#E8E2DC", "#7A1F2B"),
    dot.scale = 8,
    col.min = 0,
    col.max = 3
  ) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right") +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"),
      axis.text.y = element_text(size = 15, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      legend.title = element_text(size = 13, color = "black"),
      legend.text = element_text(size = 12, color = "black"),
      legend.position = "top",
      plot.margin = margin(20, 35, 35, 25)
    )

  ggsave(
    filename = file.path(fig_dir, "dotplot_selected_genes_HA_enriched_DLF_states_vs_other_states.png"),
    plot = p_dot,
    width = 13,
    height = 4.8,
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Targeted HA-enriched destructive lining fibroblast state DE complete.\n")
cat("Comparison  : ", ha_enriched_group, " vs ", other_state_group, "\n", sep = "")
cat("Results dir : ", res_dir, "\n")
cat("Figures dir : ", fig_dir, "\n")
cat("Main table  : ", file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states.csv"), "\n")
cat("Caution     : FindMarkers is exploratory cell-level DE.\n")
cat("============================================================\n")
