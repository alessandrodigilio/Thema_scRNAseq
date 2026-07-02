################################################
### Atlas cell-type composition: HA vs other ###
################################################

# compare final atlas cell type fractions between HA and other samples

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R")

# directories
input_object <- file.path(data_dir, "integrated_object", "annotated.rds")
fig_dir <- file.path(atlas_figures_dir, "celltype_composition_ha_vs_other")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
res_dir <- file.path(atlas_results_dir, "celltype_composition_ha_vs_other")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
condition_col <- "condition"
sample_col <- "sample_id"
celltype_col <- "cell_type"
group_levels <- c("other", "HA")
group_colors <- c("other" = "#B65A5A", "HA" = "#5B8DB8")
pseudocount <- 1e-4 # to avoid zero fractions when plotting on log scale
fdr_thr <- 0.05
n_facet_cols <- 5

# load final annotated object
cat("Loading final annotated object...\n")
obj <- readRDS(input_object)
cat("Cells:", ncol(obj), "\n")

# metadata
meta_df <- obj@meta.data[, c(condition_col, sample_col, celltype_col), drop = FALSE]
colnames(meta_df) <- c("condition", "sample_id", "cell_type") # keep the needed columns
meta_df$condition <- factor(as.character(meta_df$condition), levels = group_levels)
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$cell_type <- as.character(meta_df$cell_type)
cell_types <- unique(meta_df$cell_type)

# calculate cell type fractions for each sample
composition_df <- build_sample_composition(meta_df)
composition_df$condition <- factor(as.character(composition_df$condition), levels = group_levels)
composition_df$cell_type <- factor(as.character(composition_df$cell_type), levels = cell_types)
write.csv(composition_df[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total_cells", "fraction")], file.path(res_dir, "celltype_fraction_by_sample.csv"), row.names = FALSE)

# compare sample-level fractions between HA and other
stat_df <- test_composition_ha_vs_other(
  composition_df = composition_df,
  cell_types = cell_types,
  group_levels = group_levels,
  pseudocount = pseudocount
)
# save
write.csv(stat_df, file.path(res_dir, "celltype_fraction_ha_vs_other_wilcoxon.csv"), row.names = FALSE)

# rank cell types by effect size (absolute difference in mean fractions)
top_diff_df <- stat_df[order(abs(stat_df$diff_mean_ha_minus_other), decreasing = TRUE), , drop = FALSE]
# save
write.csv(top_diff_df, file.path(res_dir, "celltype_fraction_ha_vs_other_ranked_by_effect.csv"), row.names = FALSE)

# prepare FDR stars for the boxplot facets
annotation_df <- stat_df
annotation_df$label <- ""
annotation_df$label[annotation_df$padj_bh < fdr_thr] <- "*"
annotation_df$label[annotation_df$padj_bh < 0.01] <- "**"
annotation_df$label[annotation_df$padj_bh < 0.001] <- "***"

max_fraction_df <- aggregate(fraction ~ cell_type, data = composition_df, FUN = max)
annotation_df <- merge(annotation_df, max_fraction_df, by = "cell_type", sort = FALSE)
annotation_df$x <- 1.5
annotation_df$y <- ifelse(annotation_df$fraction > 0, annotation_df$fraction * 1.08, pseudocount * 2)

# plot HA versus other fractions for each cell type
sample_order <- unique(composition_df$sample_id[order(composition_df$condition, composition_df$sample_id)])
composition_df$sample_id <- factor(composition_df$sample_id, levels = sample_order)

# plot
p_box <- ggplot(composition_df, aes(x = condition, y = fraction, color = condition)) +
  geom_boxplot(aes(fill = condition), width = 0.55, outlier.shape = NA, color = "black") +
  geom_point(position = position_jitter(width = 0.08, height = 0), color = "black", size = 1.8, alpha = 0.9) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = n_facet_cols) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Cell type fraction", fill = NULL) +
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

ggsave(file.path(fig_dir, "celltype_fraction_ha_vs_other_boxplots.png"), p_box, width = 15, height = 10, dpi = 600)

cat("\n============================================================\n")
cat("HA vs other atlas composition analysis complete.\n")
cat("Results dir: ", res_dir, "\n", sep = "")
cat("Figures dir: ", fig_dir, "\n", sep = "")
cat("Main table : ", file.path(res_dir, "celltype_fraction_ha_vs_other_wilcoxon.csv"), "\n", sep = "")
cat("============================================================\n")
