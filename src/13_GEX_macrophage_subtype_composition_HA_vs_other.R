# ===============================================================
#   13_GEX_macrophage_subtype_composition_HA_vs_other.R
# ===============================================================

# Test whether macrophage subtype proportions differ between HA and
# other samples using the sample as the statistical unit.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
INPUT_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_annotated_macrophage_states.rds")

FIG_DIR <- file.path(FIGURES_DIR, "macrophage_subtype_composition_HA_vs_other")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RES_DIR <- file.path(RESULTS_DIR, "macrophage_subtype_composition_HA_vs_other")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# Set parameters
CONDITION_COL <- "condition"
SAMPLE_COL <- "sample_id"
CELLTYPE_COL <- "macrophage_subtype"
GROUP_LEVELS <- c("other", "HA")
GROUP_COLORS <- c("other" = "#B65A5A", "HA" = "#5B8DB8")
PSEUDOCOUNT <- 1e-4
FDR_THR <- 0.05
N_FACET_COLS <- 3

# Load final annotated object
if (!file.exists(INPUT_OBJECT)) {
  stop("Missing macrophage-annotated object: ", INPUT_OBJECT)
}

cat("Loading macrophage-annotated object...\n")
obj <- readRDS(INPUT_OBJECT)
cat("Cells:", ncol(obj), "\n")

# Check columns needed
required_cols <- c(CONDITION_COL, SAMPLE_COL, CELLTYPE_COL)
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

# Keep the metadata needed for the composition analysis
meta_df <- obj@meta.data[, required_cols, drop = FALSE]
colnames(meta_df) <- c("condition", "sample_id", "cell_type")

meta_df$condition <- factor(as.character(meta_df$condition), levels = GROUP_LEVELS)
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$cell_type <- as.character(meta_df$cell_type)

meta_df <- meta_df[
  !is.na(meta_df$condition) &
    !is.na(meta_df$sample_id) & meta_df$sample_id != "" &
    !is.na(meta_df$cell_type) & meta_df$cell_type != "",
  ,
  drop = FALSE
]

cat("Cells used for composition analysis:", nrow(meta_df), "\n")
cat("Samples by condition:\n")
print(table(meta_df$sample_id, meta_df$condition))

# Count cells and calculate subtype fractions per sample
sample_condition <- unique(meta_df[, c("sample_id", "condition"), drop = FALSE])
cell_types <- unique(meta_df$cell_type)

count_table <- table(meta_df$sample_id, meta_df$cell_type)
composition_df <- as.data.frame(count_table, stringsAsFactors = FALSE)
colnames(composition_df) <- c("sample_id", "cell_type", "n_cells")

composition_df <- merge(composition_df, sample_condition, by = "sample_id", sort = FALSE)
sample_totals <- aggregate(n_cells ~ sample_id, data = composition_df, FUN = sum)
colnames(sample_totals)[2] <- "sample_total_cells"
composition_df <- merge(composition_df, sample_totals, by = "sample_id", sort = FALSE)
composition_df$fraction <- composition_df$n_cells / composition_df$sample_total_cells
composition_df$condition <- factor(as.character(composition_df$condition), levels = GROUP_LEVELS)

cat("Composition table (first 6 rows):\n")
print(head(composition_df))

write.csv(
  composition_df[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total_cells", "fraction")],
  file.path(RES_DIR, "macrophage_subtype_fraction_by_sample.csv"),
  row.names = FALSE
)

# Compare HA and other fractions across samples for each subtype
stat_rows <- list()
for (ct in cell_types) {
  df_ct <- composition_df[composition_df$cell_type == ct, , drop = FALSE]
  other_vals <- df_ct$fraction[df_ct$condition == "other"]
  ha_vals <- df_ct$fraction[df_ct$condition == "HA"]

  mean_other <- mean(other_vals)
  mean_ha <- mean(ha_vals)

  p_wilcox <- NA_real_
  if (length(other_vals) >= 2 && length(ha_vals) >= 2) {
    p_wilcox <- wilcox.test(ha_vals, other_vals, exact = FALSE)$p.value
  }

  stat_rows[[ct]] <- data.frame(
    cell_type = ct,
    n_other_samples = length(other_vals),
    n_ha_samples = length(ha_vals),
    mean_other = mean_other,
    mean_ha = mean_ha,
    median_other = median(other_vals),
    median_ha = median(ha_vals),
    diff_mean_ha_minus_other = mean_ha - mean_other,
    log2fc_mean_ha_vs_other = log2((mean_ha + PSEUDOCOUNT) / (mean_other + PSEUDOCOUNT)),
    p_wilcox = p_wilcox,
    stringsAsFactors = FALSE
  )
}

stat_df <- do.call(rbind, stat_rows)
stat_df$padj_bh <- p.adjust(stat_df$p_wilcox, method = "BH")
stat_df <- stat_df[order(stat_df$p_wilcox), , drop = FALSE]

cat("Statistics table (first 6 rows):\n")
print(head(stat_df))

write.csv(
  stat_df,
  file.path(RES_DIR, "macrophage_subtype_fraction_HA_vs_other_wilcoxon.csv"),
  row.names = FALSE
)

top_diff_df <- stat_df[order(abs(stat_df$diff_mean_ha_minus_other), decreasing = TRUE), , drop = FALSE]
cat("Top macrophage subtypes ranked by mean difference in fraction (first 6 rows):\n")
print(head(top_diff_df))

write.csv(
  top_diff_df,
  file.path(RES_DIR, "macrophage_subtype_fraction_HA_vs_other_ranked_by_effect.csv"),
  row.names = FALSE
)

num_significant <- sum(!is.na(stat_df$padj_bh) & stat_df$padj_bh < FDR_THR)
cat("Number of macrophage subtypes with significant fraction differences at FDR < ", FDR_THR, ": ", num_significant, "\n", sep = "")
cat("Significant macrophage subtypes:\n")
print(stat_df[!is.na(stat_df$padj_bh) & stat_df$padj_bh < FDR_THR, , drop = FALSE])

# Plot HA versus other fractions for each subtype
sample_order <- sample_condition[order(sample_condition$condition, sample_condition$sample_id), "sample_id"]
composition_df$sample_id <- factor(composition_df$sample_id, levels = sample_order)

annotation_df <- stat_df
annotation_df$label <- ""
annotation_df$label[!is.na(annotation_df$padj_bh) & annotation_df$padj_bh < FDR_THR] <- "*"
annotation_df$label[!is.na(annotation_df$padj_bh) & annotation_df$padj_bh < 0.01] <- "**"
annotation_df$label[!is.na(annotation_df$padj_bh) & annotation_df$padj_bh < 0.001] <- "***"
max_fraction_df <- aggregate(fraction ~ cell_type, data = composition_df, FUN = max)
annotation_df <- merge(annotation_df, max_fraction_df, by = "cell_type", sort = FALSE)
annotation_df$x <- 1.5
annotation_df$y <- ifelse(annotation_df$fraction > 0, annotation_df$fraction * 1.08, PSEUDOCOUNT * 2)

p_box <- ggplot(composition_df, aes(x = condition, y = fraction, color = condition)) +
  geom_boxplot(aes(fill = condition), width = 0.55, outlier.shape = NA, color = "black") +
  geom_point(position = position_jitter(width = 0.08, height = 0), color = "black", size = 1.8, alpha = 0.9) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = N_FACET_COLS) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(x = NULL, y = "Macrophage subtype fraction", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  geom_text(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 10
  )

print(p_box)

ggsave(
  filename = file.path(FIG_DIR, "macrophage_subtype_fraction_HA_vs_other_boxplots.png"),
  plot = p_box,
  width = 12,
  height = 9,
  dpi = 600
)

cat("\n============================================================\n")
cat("HA vs other macrophage composition analysis complete.\n")
cat("Results dir : ", RES_DIR, "\n")
cat("Figures dir : ", FIG_DIR, "\n")
cat("Main table  : ", file.path(RES_DIR, "macrophage_subtype_fraction_HA_vs_other_wilcoxon.csv"), "\n")
cat("============================================================\n")
