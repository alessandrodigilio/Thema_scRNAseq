# ---------------------------------------
# 23_GEX_DE_HA_enriched_destructive_lining_fibroblast_states.R
# Targeted DE of HA-enriched MMP3+ lining fibroblast states
# ---------------------------------------

# Compare the two HA-enriched destructive/MMP3+ lining fibroblast
# subclusters against the remaining destructive lining fibroblast states.
# This is exploratory cell-level DE, followed by sample-level summaries
# to check whether key signals are driven by only one patient.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(readxl)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
INPUT_SUBSET_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_destructive_lining_fibroblasts_subclustered.rds")

RES_DIR <- file.path(RESULTS_DIR, "destructive_lining_fibroblast_HA_enriched_state_DE")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(FIGURES_DIR, "destructive_lining_fibroblast_HA_enriched_state_DE")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Set parameters
CLUSTER_COL <- "destructive_lining_fibroblast_subcluster"
LABEL_COL <- "destructive_lining_fibroblast_subtype"
STATE_GROUP_COL <- "destructive_lining_fibroblast_state_group"
CONDITION_COL <- "condition"
SAMPLE_COL <- "sample_id"

HA_ENRICHED_CLUSTERS <- c("2", "5")
OTHER_STATE_CLUSTERS <- c("0", "1", "3", "4")
HA_ENRICHED_GROUP <- "HA-enriched DLF states"
OTHER_STATE_GROUP <- "Other DLF states"
GROUP_LEVELS <- c(OTHER_STATE_GROUP, HA_ENRICHED_GROUP)

MIN_PCT <- 0.10
LOGFC_THR <- 0.10
PADJ_THR <- 0.05
TOP_N_LABELS <- 20
MIN_CELLS_PER_SAMPLE_STATE <- 5
GROUP_COLORS <- c(
  "Other DLF states" = "#8FA0A8",
  "HA-enriched DLF states" = "#C65A5A"
)
CONDITION_COLORS <- c(
  "other" = "#B65A5A",
  "HA" = "#5B8DB8"
)

DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS <- c(
  "0" = "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)",
  "1" = "Activated MMP3+ lining fibroblast cells (ID1+)",
  "2" = "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)",
  "3" = "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)",
  "4" = "MMP3+ lining fibroblast cells (FAM184A+)",
  "5" = "HA-enriched SFRP2+ matrix fibroblast-like cells"
)

IRON_GENE_FILES <- c(
  file.path(METADATA_DIR, "iron_genes", "ferroptosis_genes_curated.xlsx"),
  file.path(METADATA_DIR, "iron_genes", "iron_uptake_transport_genes.xlsx")
)

SELECTED_GENES <- c(
  "HMOX1", "NQO1", "FTL", "FTH1", "CP", "SLC40A1", "TFRC",
  "SFRP2", "SFRP1", "COMP", "PODN", "IGF1", "MFAP5",
  "CCL7", "CXCL1", "CCL20", "CCRL2", "BIRC3",
  "MMP3", "MMP1", "IL6", "PTGS2"
)

load_iron_genes <- function(files) {
  iron_genes <- character(0)

  for (f in files) {
    if (!file.exists(f)) next
    x <- readxl::read_excel(f)
    if (ncol(x) == 0) next

    genes_here <- trimws(as.character(x[[1]]))
    genes_here <- toupper(genes_here[!is.na(genes_here) & genes_here != ""])
    iron_genes <- c(iron_genes, genes_here)
  }

  iron_genes <- sort(unique(iron_genes))
  iron_genes <- iron_genes[iron_genes != "GENE"]
  iron_genes[iron_genes == "TRFC"] <- "TFRC"
  iron_genes
}

safe_neg_log10 <- function(x) {
  y <- -log10(x)
  if (any(is.infinite(y), na.rm = TRUE)) {
    max_finite <- max(y[is.finite(y)], na.rm = TRUE)
    if (!is.finite(max_finite)) max_finite <- 0
    y[is.infinite(y)] <- max_finite + 1
  }
  y
}

# Load destructive lining fibroblast subset object
if (!file.exists(INPUT_SUBSET_OBJECT)) {
  stop("Missing destructive lining fibroblast subset object: ", INPUT_SUBSET_OBJECT)
}

cat("Loading destructive lining fibroblast subset object...\n")
obj <- readRDS(INPUT_SUBSET_OBJECT)
cat("Cells:", ncol(obj), "\n")

required_cols <- c(CLUSTER_COL, CONDITION_COL, SAMPLE_COL)
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

# Apply final labels and define the targeted state groups
cluster_ids <- as.character(obj@meta.data[[CLUSTER_COL]])
mapped_labels <- unname(DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS[cluster_ids])
mapped_labels[is.na(mapped_labels)] <- "Unknown"

obj@meta.data[[LABEL_COL]] <- factor(
  mapped_labels,
  levels = unname(DESTRUCTIVE_LINING_FIBROBLAST_SUBCLUSTER_LABELS)
)
obj@meta.data[[LABEL_COL]] <- droplevels(obj@meta.data[[LABEL_COL]])

state_group <- rep(NA_character_, length(cluster_ids))
state_group[cluster_ids %in% HA_ENRICHED_CLUSTERS] <- HA_ENRICHED_GROUP
state_group[cluster_ids %in% OTHER_STATE_CLUSTERS] <- OTHER_STATE_GROUP

obj@meta.data[[STATE_GROUP_COL]] <- factor(state_group, levels = GROUP_LEVELS)
obj <- subset(obj, cells = rownames(obj@meta.data)[!is.na(obj@meta.data[[STATE_GROUP_COL]])])
obj@meta.data[[STATE_GROUP_COL]] <- droplevels(obj@meta.data[[STATE_GROUP_COL]])

cat("Cells by final subtype:\n")
print(table(obj@meta.data[[LABEL_COL]], useNA = "ifany"))
cat("Cells by targeted state group:\n")
print(table(obj@meta.data[[STATE_GROUP_COL]], useNA = "ifany"))
cat("Cells by sample and targeted state group:\n")
print(table(obj@meta.data[[SAMPLE_COL]], obj@meta.data[[STATE_GROUP_COL]]))

# Save sample-level composition of HA-enriched versus other DLF states
meta_df <- obj@meta.data[, c(SAMPLE_COL, CONDITION_COL, CLUSTER_COL, LABEL_COL, STATE_GROUP_COL), drop = FALSE]
colnames(meta_df) <- c("sample_id", "condition", "subcluster", "subtype", "state_group")
meta_df$condition <- ifelse(as.character(meta_df$condition) == "HA", "HA", "other")
meta_df$condition <- factor(meta_df$condition, levels = c("other", "HA"))
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$state_group <- factor(as.character(meta_df$state_group), levels = GROUP_LEVELS)

composition_df <- as.data.frame(table(meta_df$sample_id, meta_df$condition, meta_df$state_group), stringsAsFactors = FALSE)
colnames(composition_df) <- c("sample_id", "condition", "state_group", "n_cells")
composition_df <- composition_df[composition_df$n_cells > 0, , drop = FALSE]
sample_totals <- aggregate(n_cells ~ sample_id, data = composition_df, FUN = sum)
colnames(sample_totals)[2] <- "sample_total_cells"
composition_df <- merge(composition_df, sample_totals, by = "sample_id", sort = FALSE)
composition_df$fraction <- composition_df$n_cells / composition_df$sample_total_cells
composition_df$condition <- factor(as.character(composition_df$condition), levels = c("other", "HA"))
composition_df$state_group <- factor(as.character(composition_df$state_group), levels = GROUP_LEVELS)

write.csv(
  composition_df,
  file.path(RES_DIR, "sample_composition_HA_enriched_DLF_states_vs_other_states.csv"),
  row.names = FALSE
)

sample_order <- unique(composition_df[order(composition_df$condition, composition_df$sample_id), "sample_id"])
composition_df$sample_id <- factor(composition_df$sample_id, levels = sample_order)

p_composition <- ggplot(composition_df, aes(x = sample_id, y = fraction, fill = state_group)) +
  geom_col(width = 0.88, color = "black", linewidth = 0.25) +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
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
  filename = file.path(FIG_DIR, "sample_fraction_HA_enriched_DLF_states_vs_other_states.png"),
  plot = p_composition,
  width = 11,
  height = 5.8,
  dpi = 600
)

# Run exploratory cell-level differential expression
cat("\nRunning FindMarkers: HA-enriched DLF states vs other DLF states...\n")
Idents(obj) <- obj@meta.data[[STATE_GROUP_COL]]

de_df <- FindMarkers(
  object = obj,
  assay = "RNA",
  ident.1 = HA_ENRICHED_GROUP,
  ident.2 = OTHER_STATE_GROUP,
  test.use = "wilcox",
  min.pct = MIN_PCT,
  logfc.threshold = LOGFC_THR,
  only.pos = FALSE,
  verbose = FALSE
)

de_df$gene <- rownames(de_df)
de_df <- de_df[, c("gene", setdiff(colnames(de_df), "gene"))]
de_df <- de_df[order(de_df$p_val_adj, de_df$p_val), , drop = FALSE]
de_df$direction <- ifelse(de_df$avg_log2FC > 0, "Higher in HA-enriched DLF states", "Higher in other DLF states")
de_df$is_significant <- !is.na(de_df$p_val_adj) & de_df$p_val_adj < PADJ_THR

write.csv(
  de_df,
  file.path(RES_DIR, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states.csv"),
  row.names = FALSE
)

sig_df <- de_df[de_df$is_significant, , drop = FALSE]
write.csv(
  sig_df,
  file.path(RES_DIR, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_significant.csv"),
  row.names = FALSE
)

iron_genes <- load_iron_genes(IRON_GENE_FILES)
iron_de_df <- de_df[toupper(de_df$gene) %in% iron_genes, , drop = FALSE]
write.csv(
  iron_de_df,
  file.path(RES_DIR, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_iron_related.csv"),
  row.names = FALSE
)

selected_de_df <- de_df[de_df$gene %in% SELECTED_GENES, , drop = FALSE]
write.csv(
  selected_de_df,
  file.path(RES_DIR, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_selected_genes.csv"),
  row.names = FALSE
)

cat("Top DE genes:\n")
print(head(de_df, 20))
cat("Significant genes at adjusted p < ", PADJ_THR, ": ", nrow(sig_df), "\n", sep = "")
cat("Iron-related genes found in DE table:", nrow(iron_de_df), "\n")
cat("Significant iron-related genes:\n")
print(iron_de_df[iron_de_df$is_significant, c("gene", "avg_log2FC", "p_val_adj", "direction"), drop = FALSE])

# Volcano plot for the targeted state comparison
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
top_labels <- head(top_labels, TOP_N_LABELS)

p_volcano <- ggplot(plot_df, aes(x = avg_log2FC, y = neg_log10_padj, color = volcano_group)) +
  geom_point(size = 1.7, alpha = 0.82) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.65, color = "black") +
  geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed", linewidth = 0.65, color = "black") +
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
      "Higher in other DLF states" = unname(GROUP_COLORS[OTHER_STATE_GROUP]),
      "Not significant" = "grey76",
      "Higher in HA-enriched DLF states" = unname(GROUP_COLORS[HA_ENRICHED_GROUP])
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
  filename = file.path(FIG_DIR, "volcano_HA_enriched_DLF_states_vs_other_DLF_states.png"),
  plot = p_volcano,
  width = 9,
  height = 7.5,
  dpi = 600
)

# Dotplot of selected genes across the two broad state groups
selected_genes_use <- SELECTED_GENES[SELECTED_GENES %in% rownames(obj)]

if (length(selected_genes_use) > 0) {
  p_dot <- DotPlot(
    object = obj,
    features = selected_genes_use,
    group.by = STATE_GROUP_COL,
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
    filename = file.path(FIG_DIR, "dotplot_selected_genes_HA_enriched_DLF_states_vs_other_states.png"),
    plot = p_dot,
    width = 13,
    height = 4.8,
    dpi = 600
  )
}

# Build sample-level expression summaries for selected genes
selected_for_summary <- selected_genes_use

if (length(selected_for_summary) > 0) {
  expr_df <- FetchData(obj, vars = c(selected_for_summary, SAMPLE_COL, CONDITION_COL, STATE_GROUP_COL))
  expr_df$sample_id <- as.character(expr_df[[SAMPLE_COL]])
  expr_df$condition <- ifelse(as.character(expr_df[[CONDITION_COL]]) == "HA", "HA", "other")
  expr_df$state_group <- as.character(expr_df[[STATE_GROUP_COL]])

  sample_summary_rows <- list()
  row_i <- 1

  for (gene in selected_for_summary) {
    for (sample_here in unique(expr_df$sample_id)) {
      for (group_here in GROUP_LEVELS) {
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
  sample_summary_df$state_group <- factor(sample_summary_df$state_group, levels = GROUP_LEVELS)
  sample_summary_df$is_low_cell_count <- sample_summary_df$n_cells < MIN_CELLS_PER_SAMPLE_STATE

  write.csv(
    sample_summary_df,
    file.path(RES_DIR, "sample_level_expression_selected_genes_HA_enriched_DLF_states_vs_other_states.csv"),
    row.names = FALSE
  )

  sample_stat_rows <- list()
  stat_i <- 1

  for (gene in selected_for_summary) {
    df_gene <- sample_summary_df[
      sample_summary_df$gene == gene &
        !sample_summary_df$is_low_cell_count,
      ,
      drop = FALSE
    ]

    x <- df_gene$mean_expression[df_gene$state_group == HA_ENRICHED_GROUP]
    y <- df_gene$mean_expression[df_gene$state_group == OTHER_STATE_GROUP]

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

  write.csv(
    sample_stat_df,
    file.path(RES_DIR, "sample_level_expression_selected_genes_wilcoxon_HA_enriched_DLF_states_vs_other_states.csv"),
    row.names = FALSE
  )

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
      scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
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
      filename = file.path(FIG_DIR, "sample_level_boxplots_key_genes_HA_enriched_DLF_states_vs_other_states.png"),
      plot = p_box,
      width = 13,
      height = 7.5,
      dpi = 600
    )
  }
}

cat("\n============================================================\n")
cat("Targeted HA-enriched destructive lining fibroblast state DE complete.\n")
cat("Comparison  : ", HA_ENRICHED_GROUP, " vs ", OTHER_STATE_GROUP, "\n", sep = "")
cat("Results dir : ", RES_DIR, "\n")
cat("Figures dir : ", FIG_DIR, "\n")
cat("Main table  : ", file.path(RES_DIR, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states.csv"), "\n")
cat("Caution     : FindMarkers is exploratory cell-level DE; sample-level summaries are included for patient-level sanity checks.\n")
cat("============================================================\n")
