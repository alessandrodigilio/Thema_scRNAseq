###################################################
### Macrophage subtype composition: HA vs other ###
###################################################

# compare macrophage subtype fractions between HA and other samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")
source("src/macrophage_subclusters/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "annotated_macrophage_states.rds")
fig_dir <- file.path(figures_dir, "macrophage_subtype_composition_ha_vs_other")
res_dir <- file.path(results_dir, "macrophage_subtype_composition_ha_vs_other")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
condition_col <- "condition"
sample_col <- "sample_id"
subtype_col <- "macrophage_subtype"
group_levels <- c("other", "HA")
group_colors <- c("other" = "#B65A5A", "HA" = "#5B8DB8")
pseudocount <- 1e-4
fdr_thr <- 0.05
n_facet_cols <- 3

# load macrophage-annotated object
cat("Loading macrophage-annotated object...\n")
obj <- readRDS(input_object)
cat("Cells:", ncol(obj), "\n")

# metadata
meta_df <- obj@meta.data[, c(condition_col, sample_col, subtype_col), drop = FALSE]
colnames(meta_df) <- c("condition", "sample_id", "cell_type")
meta_df$condition <- factor(ifelse(as.character(meta_df$condition) == "HA", "HA", "other"), levels = group_levels)
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$cell_type <- as.character(meta_df$cell_type)
meta_df <- meta_df[!is.na(meta_df$sample_id) & meta_df$sample_id != "" & !is.na(meta_df$cell_type) & meta_df$cell_type != "", , drop = FALSE]

cell_types <- names(macrophage_subtype_colors)
cell_types <- cell_types[cell_types %in% unique(meta_df$cell_type)]

# composition by sample and statistics across samples
composition_df <- build_sample_composition(meta_df)
stat_df <- test_composition_ha_vs_other(composition_df, cell_types, group_levels, pseudocount)
top_diff_df <- stat_df[order(abs(stat_df$diff_mean_ha_minus_other), decreasing = TRUE), , drop = FALSE]
# save stats
write.csv(composition_df[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total_cells", "fraction")], file.path(res_dir, "macrophage_subtype_fraction_by_sample.csv"), row.names = FALSE)
write.csv(stat_df, file.path(res_dir, "macrophage_subtype_fraction_HA_vs_other_wilcoxon.csv"), row.names = FALSE)
write.csv(top_diff_df, file.path(res_dir, "macrophage_subtype_fraction_HA_vs_other_ranked_by_effect.csv"), row.names = FALSE)

# plot HA versus other fractions for each subtype
sample_condition <- unique(meta_df[, c("sample_id", "condition"), drop = FALSE])
sample_order <- sample_condition[order(sample_condition$condition, sample_condition$sample_id), "sample_id"]
composition_df$sample_id <- factor(composition_df$sample_id, levels = sample_order)

# add significance annotation
annotation_df <- stat_df
annotation_df$label <- ""
annotation_df$label[!is.na(annotation_df$padj_bh) & annotation_df$padj_bh < fdr_thr] <- "*"
annotation_df$label[!is.na(annotation_df$padj_bh) & annotation_df$padj_bh < 0.01] <- "**"
annotation_df$label[!is.na(annotation_df$padj_bh) & annotation_df$padj_bh < 0.001] <- "***"
max_fraction_df <- aggregate(fraction ~ cell_type, data = composition_df, FUN = max)
annotation_df <- merge(annotation_df, max_fraction_df, by = "cell_type", sort = FALSE)
annotation_df$x <- 1.5
annotation_df$y <- ifelse(annotation_df$fraction > 0, annotation_df$fraction * 1.08, pseudocount * 2)

# plot boxplots
p_box <- ggplot(composition_df, aes(x = condition, y = fraction, color = condition)) +
  geom_boxplot(aes(fill = condition), width = 0.55, outlier.shape = NA, color = "black") +
  geom_point(position = position_jitter(width = 0.08, height = 0), color = "black", size = 1.8, alpha = 0.9) +
  geom_text(data = annotation_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, color = "black", size = 10) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = n_facet_cols) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Macrophage subtype fraction", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(file.path(fig_dir, "macrophage_subtype_fraction_HA_vs_other_boxplots.png"), p_box, width = 12, height = 9, dpi = 600)

cat("\n============================================================\n")
cat("Macrophage composition analysis complete.\n")
cat("Results dir: ", res_dir, "\n", sep = "")
cat("Figures dir: ", fig_dir, "\n", sep = "")
cat("============================================================\n")
