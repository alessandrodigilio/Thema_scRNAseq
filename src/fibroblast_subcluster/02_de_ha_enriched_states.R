##############################################################
### DE of HA-enriched destructive lining fibroblast states ###
##############################################################

# compare the two HA-enriched destructive/MMP3+ lining fibroblast
# subclusters against the remaining destructive lining fibroblast states.
# this is exploratory cell-level DE, followed by sample-level summaries
# to check whether key signals are driven by only one patient.

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

ha_enriched_clusters <- c("2", "5")
other_state_clusters <- c("0", "1", "3", "4")
ha_enriched_group <- "HA-enriched DLF states"
other_state_group <- "Other DLF states"
group_levels <- c(other_state_group, ha_enriched_group)

padj_thr <- 0.05
top_n_labels <- 20
min_cells_per_sample_state <- 5
group_colors <- c(
  "Other DLF states" = "#8FA0A8",
  "HA-enriched DLF states" = "#C65A5A"
)

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

cat("Cells by final subtype:\n")
print(table(obj@meta.data[[label_col]], useNA = "ifany"))
cat("Cells by targeted state group:\n")
print(table(obj@meta.data[[state_group_col]], useNA = "ifany"))
cat("Cells by sample and targeted state group:\n")
print(table(obj@meta.data[[sample_col]], obj@meta.data[[state_group_col]]))

# save sample-level composition of HA-enriched versus other DLF states
meta_df <- obj@meta.data[, c(sample_col, condition_col, cluster_col, label_col, state_group_col), drop = FALSE]
colnames(meta_df) <- c("sample_id", "condition", "subcluster", "subtype", "state_group")
meta_df$condition <- ifelse(as.character(meta_df$condition) == "HA", "HA", "other")
meta_df$condition <- factor(meta_df$condition, levels = c("other", "HA"))
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$state_group <- factor(as.character(meta_df$state_group), levels = group_levels)

composition_df <- as.data.frame(table(meta_df$sample_id, meta_df$condition, meta_df$state_group), stringsAsFactors = FALSE)
colnames(composition_df) <- c("sample_id", "condition", "state_group", "n_cells")
composition_df <- composition_df[composition_df$n_cells > 0, , drop = FALSE]
sample_totals <- aggregate(n_cells ~ sample_id, data = composition_df, FUN = sum)
colnames(sample_totals)[2] <- "sample_total_cells"
composition_df <- merge(composition_df, sample_totals, by = "sample_id", sort = FALSE)
composition_df$fraction <- composition_df$n_cells / composition_df$sample_total_cells
composition_df$condition <- factor(as.character(composition_df$condition), levels = c("other", "HA"))
composition_df$state_group <- factor(as.character(composition_df$state_group), levels = group_levels)

write.csv(composition_df, file.path(res_dir, "sample_composition_HA_enriched_DLF_states_vs_other_states.csv"), row.names = FALSE)

sample_order <- unique(composition_df[order(composition_df$condition, composition_df$sample_id), "sample_id"])
composition_df$sample_id <- factor(composition_df$sample_id, levels = sample_order)

p_composition <- ggplot(composition_df, aes(x = sample_id, y = fraction, fill = state_group)) +
  geom_col(width = 0.88, color = "black", linewidth = 0.25) +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction", fill = NULL) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 17, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.line = element_line(linewidth = 1.1, color = "black"),
    axis.ticks = element_line(linewidth = 1.1, color = "black"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 17, face = "bold", color = "black", margin = margin(t = 8)),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    panel.spacing.x = grid::unit(0.35, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 30, 20)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  filename = file.path(fig_dir, "sample_fraction_HA_enriched_DLF_states_vs_other_states.png"),
  plot = p_composition,
  width = 11,
  height = 5.8,
  dpi = 600
)

# run exploratory cell-level differential expression
cat("\nRunning FindMarkers: HA-enriched DLF states vs other DLF states...\n")
Idents(obj) <- obj@meta.data[[state_group_col]]

de_df <- FindMarkers(
  object = obj,
  assay = "RNA",
  ident.1 = ha_enriched_group,
  ident.2 = other_state_group,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE,
  verbose = FALSE
)

de_df$gene <- rownames(de_df)
de_df <- de_df[, c("gene", setdiff(colnames(de_df), "gene"))]
de_df <- de_df[order(de_df$p_val_adj, de_df$p_val), , drop = FALSE]
de_df$direction <- ifelse(de_df$avg_log2FC > 0, "Higher in HA-enriched DLF states", "Higher in other DLF states")
de_df$is_significant <- !is.na(de_df$p_val_adj) & de_df$p_val_adj < padj_thr

write.csv(de_df, file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states.csv"), row.names = FALSE)

sig_df <- de_df[de_df$is_significant, , drop = FALSE]
write.csv(sig_df, file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_significant.csv"), row.names = FALSE)

iron_genes <- load_gene_set_files(iron_related_geneset_files)
iron_de_df <- de_df[toupper(de_df$gene) %in% iron_genes, , drop = FALSE]
write.csv(iron_de_df, file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_iron_related.csv"), row.names = FALSE)

selected_de_df <- de_df[de_df$gene %in% selected_genes, , drop = FALSE]
write.csv(selected_de_df, file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_selected_genes.csv"), row.names = FALSE)

cat("Top DE genes:\n")
print(head(de_df, 20))
cat("Significant genes at adjusted p < ", padj_thr, ": ", nrow(sig_df), "\n", sep = "")
cat("Iron-related genes found in DE table:", nrow(iron_de_df), "\n")
cat("Significant iron-related genes:\n")
print(iron_de_df[iron_de_df$is_significant, c("gene", "avg_log2FC", "p_val_adj", "direction"), drop = FALSE])

# volcano plot for the targeted state comparison
plot_df <- de_df[
  !is.na(de_df$avg_log2FC) &
    !is.na(de_df$p_val_adj),
  ,
  drop = FALSE
]
plot_df$neg_log10_padj <- safe_neg_log10(plot_df$p_val_adj)
plot_df$volcano_group <- "Not significant"
plot_df$volcano_group[plot_df$is_significant & plot_df$avg_log2FC > 0] <- "Higher in HA-enriched DLF states"
plot_df$volcano_group[plot_df$is_significant & plot_df$avg_log2FC < 0] <- "Higher in other DLF states"
plot_df$volcano_group <- factor(
  plot_df$volcano_group,
  levels = c("Higher in other DLF states", "Not significant", "Higher in HA-enriched DLF states")
)

top_labels <- plot_df[plot_df$is_significant, , drop = FALSE]
top_labels <- top_labels[order(top_labels$p_val_adj, -abs(top_labels$avg_log2FC)), , drop = FALSE]
top_labels <- head(top_labels, top_n_labels)

p_volcano <- ggplot(plot_df, aes(x = avg_log2FC, y = neg_log10_padj, color = volcano_group)) +
  geom_point(size = 1.7, alpha = 0.82) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.65, color = "black") +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", linewidth = 0.65, color = "black") +
  geom_text(
    data = top_labels,
    aes(label = gene),
    color = "black",
    size = 4.2,
    vjust = -0.55,
    check_overlap = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Higher in other DLF states" = unname(group_colors[other_state_group]),
      "Not significant" = "grey76",
      "Higher in HA-enriched DLF states" = unname(group_colors[ha_enriched_group])
    )
  ) +
  labs(
    x = "avg log2FC",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.margin = margin(18, 18, 18, 18)
  )

ggsave(
  filename = file.path(fig_dir, "volcano_HA_enriched_DLF_states_vs_other_DLF_states.png"),
  plot = p_volcano,
  width = 9,
  height = 7.5,
  dpi = 600
)

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

# build sample-level expression summaries for selected genes
selected_for_summary <- selected_genes_use

if (length(selected_for_summary) > 0) {
  expr_df <- FetchData(obj, vars = c(selected_for_summary, sample_col, condition_col, state_group_col))
  expr_df$sample_id <- as.character(expr_df[[sample_col]])
  expr_df$condition <- ifelse(as.character(expr_df[[condition_col]]) == "HA", "HA", "other")
  expr_df$state_group <- as.character(expr_df[[state_group_col]])

  sample_summary_rows <- list()
  row_i <- 1

  for (gene in selected_for_summary) {
    for (sample_here in unique(expr_df$sample_id)) {
      for (group_here in group_levels) {
        idx <- expr_df$sample_id == sample_here & expr_df$state_group == group_here
        if (sum(idx) == 0) next

        condition_here <- unique(expr_df$condition[idx])[1]
        values <- expr_df[[gene]][idx]

        sample_summary_rows[[row_i]] <- data.frame(
          gene = gene,
          sample_id = sample_here,
          condition = condition_here,
          state_group = group_here,
          n_cells = length(values),
          mean_expression = mean(values, na.rm = TRUE),
          median_expression = median(values, na.rm = TRUE),
          pct_expressing = mean(values > 0, na.rm = TRUE) * 100,
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1
      }
    }
  }

  sample_summary_df <- do.call(rbind, sample_summary_rows)
  sample_summary_df$condition <- factor(sample_summary_df$condition, levels = c("other", "HA"))
  sample_summary_df$state_group <- factor(sample_summary_df$state_group, levels = group_levels)
  sample_summary_df$is_low_cell_count <- sample_summary_df$n_cells < min_cells_per_sample_state

  write.csv(sample_summary_df, file.path(res_dir, "sample_level_expression_selected_genes_HA_enriched_DLF_states_vs_other_states.csv"), row.names = FALSE)

  sample_stat_rows <- list()
  stat_i <- 1

  for (gene in selected_for_summary) {
    df_gene <- sample_summary_df[
      sample_summary_df$gene == gene &
        !sample_summary_df$is_low_cell_count,
      ,
      drop = FALSE
    ]

    x <- df_gene$mean_expression[df_gene$state_group == ha_enriched_group]
    y <- df_gene$mean_expression[df_gene$state_group == other_state_group]

    p_value <- NA_real_
    if (length(x) >= 2 && length(y) >= 2) {
      p_value <- wilcox.test(x, y, exact = FALSE)$p.value
    }

    sample_stat_rows[[stat_i]] <- data.frame(
      gene = gene,
      n_samples_HA_enriched_states = length(x),
      n_samples_other_states = length(y),
      mean_HA_enriched_states = mean(x, na.rm = TRUE),
      mean_other_states = mean(y, na.rm = TRUE),
      diff_mean_HA_enriched_minus_other = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      p_wilcox = p_value,
      stringsAsFactors = FALSE
    )
    stat_i <- stat_i + 1
  }

  sample_stat_df <- do.call(rbind, sample_stat_rows)
  sample_stat_df$padj_bh <- p.adjust(sample_stat_df$p_wilcox, method = "BH")
  sample_stat_df <- sample_stat_df[order(sample_stat_df$p_wilcox), , drop = FALSE]

  write.csv(sample_stat_df, file.path(res_dir, "sample_level_expression_selected_genes_wilcoxon_HA_enriched_DLF_states_vs_other_states.csv"), row.names = FALSE)

  cat("\nSample-level selected gene summary statistics:\n")
  print(sample_stat_df)

  key_genes_for_plot <- selected_for_summary[selected_for_summary %in% c(
    "HMOX1", "NQO1", "FTL", "FTH1", "CP", "SLC40A1",
    "SFRP2", "CCL7", "CXCL1", "BIRC3"
  )]

  box_df <- sample_summary_df[
    sample_summary_df$gene %in% key_genes_for_plot &
      !sample_summary_df$is_low_cell_count,
    ,
    drop = FALSE
  ]

  if (nrow(box_df) > 0) {
    box_df$gene <- factor(box_df$gene, levels = key_genes_for_plot)

    p_box <- ggplot(box_df, aes(x = state_group, y = mean_expression, fill = state_group)) +
      geom_boxplot(width = 0.52, outlier.shape = NA, color = "black", linewidth = 0.5) +
      geom_point(aes(shape = condition), position = position_jitter(width = 0.08, height = 0), size = 2.2, color = "black") +
      facet_wrap(~ gene, scales = "free_y", ncol = 5) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      scale_shape_manual(values = c("other" = 21, "HA" = 24), drop = FALSE) +
      labs(x = NULL, y = "Mean normalized expression", fill = NULL, shape = NULL) +
      theme_classic(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 13, face = "bold", color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        plot.margin = margin(18, 18, 18, 18)
      )

    ggsave(
      filename = file.path(fig_dir, "sample_level_boxplots_key_genes_HA_enriched_DLF_states_vs_other_states.png"),
      plot = p_box,
      width = 13,
      height = 7.5,
      dpi = 600
    )
  }
}

cat("\n============================================================\n")
cat("Targeted HA-enriched destructive lining fibroblast state DE complete.\n")
cat("Comparison  : ", ha_enriched_group, " vs ", other_state_group, "\n", sep = "")
cat("Results dir : ", res_dir, "\n")
cat("Figures dir : ", fig_dir, "\n")
cat("Main table  : ", file.path(res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states.csv"), "\n")
cat("Caution     : FindMarkers is exploratory cell-level DE; sample-level summaries are included for patient-level sanity checks.\n")
cat("============================================================\n")
