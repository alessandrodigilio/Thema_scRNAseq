#############################################################
### Macrophage subtype pseudobulk differential expression ###
#############################################################

# pseudobulk differential gene expression between HA and all other
# samples for each macrophage subtype.
# the sample is used as the statistical unit and raw RNA counts are
# summed within each sample x macrophage subtype group.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
  library(readxl)
})

setwd("~/Thema_R")
source("src/global_config.R")

# create output directories used by this step
input_object <- file.path(data_dir, "integrated_object", "annotated_macrophage_states.rds")

res_dir <- file.path(results_dir, "pseudobulk_deseq2_macrophage_subtypes")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(figures_dir, "pseudobulk_deseq2_macrophage_subtypes")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pca_dir <- file.path(fig_dir, "pca")
pca_condition_dir <- file.path(pca_dir, "pca_condition")
pca_n_cells_dir <- file.path(pca_dir, "pca_n_cells")
pca_ncount_dir <- file.path(pca_dir, "pca_nCount_RNA")
pca_nfeature_dir <- file.path(pca_dir, "pca_nFeature_RNA")
pca_mt_dir <- file.path(pca_dir, "pca_mt")
volcano_dir <- file.path(fig_dir, "volcano_plot")

dir.create(pca_condition_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_n_cells_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_ncount_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_nfeature_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_mt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
condition_col <- "condition"
sample_col <- "sample_id"
celltype_col <- "macrophage_subtype"
group_levels <- c("other", "HA")
qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
group_colors <- c("other" = "#B65A5A", "HA" = "#5B8DB8")

macrophage_subtype_colors <- c(
  "Inflammatory macrophages (KANK1+)" = "#C65A5A",
  "Inflammatory macrophages (THBS1+)" = "#b0e17b",
  "Macrophage-like state (AMTN+)" = "#D98F5C",
  "Resident macrophages (HSPA6+)" = "#8c674b",
  "Red-pulp-like resident macrophages (MERTK+)" = "#4C9F8A",
  "Mixed macrophage-like cells (RNASE1+)" = "#5f90b3",
  "Plasma-like contaminants" = "#B58ACF",
  "Low-confidence cells" = "#9FA4A9",
  "Proliferating macrophages" = "#D95FA7"
)

# pseudobulk filtering thresholds
min_cells_per_sample <- 10
min_pseudobulk_count <- 10
min_detected_cell_fraction <- 0.05
padj_thr <- 0.05
top_n_labels <- 15
log2fc_plot_limit <- 5
remove_mt_genes <- FALSE
iron_gene_files <- c(
  file.path(metadata_dir, "iron_genes", "ferroptosis_genes_curated.xlsx"),
  file.path(metadata_dir, "iron_genes", "iron_uptake_transport_genes.xlsx")
)
iron_bubble_plot_file <- file.path(fig_dir, "bubble_heatmap_iron_related_significant_genes_macrophage_subtypes.png")

# helper to sanitize subgroup names for files
sanitize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

load_iron_genes <- function(files) {
  iron_genes <- character(0)

  for (f in files) {
    if (!file.exists(f)) next
    x <- read_excel(f)
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

# load final macrophage-annotated object
if (!file.exists(input_object)) {
  stop("Missing macrophage-annotated object: ", input_object)
}

cat("Loading macrophage-annotated object...\n")
obj <- readRDS(input_object)
cat("Cells:", ncol(obj), "\n")

# check required metadata columns
required_cols <- c(condition_col, sample_col, celltype_col)
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

# use RNA assay counts for pseudobulk
DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) obj)
rna_counts <- LayerData(obj[["RNA"]], layer = "counts")

# prepare metadata
meta_cols <- unique(c(required_cols, qc_cols))
meta_cols <- meta_cols[meta_cols %in% colnames(obj@meta.data)]
meta_df <- obj@meta.data[, meta_cols, drop = FALSE]
colnames(meta_df)[colnames(meta_df) == condition_col] <- "condition"
colnames(meta_df)[colnames(meta_df) == sample_col] <- "sample_id"
colnames(meta_df)[colnames(meta_df) == celltype_col] <- "macrophage_subtype"

meta_df$condition <- ifelse(as.character(meta_df$condition) == "HA", "HA", "other")
meta_df$condition <- factor(meta_df$condition, levels = group_levels)
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$macrophage_subtype <- as.character(meta_df$macrophage_subtype)

meta_df <- meta_df[
  !is.na(meta_df$condition) &
    !is.na(meta_df$sample_id) & meta_df$sample_id != "" &
    !is.na(meta_df$macrophage_subtype) & meta_df$macrophage_subtype != "",
  ,
  drop = FALSE
]

cat("Cells used for macrophage pseudobulk analysis:", nrow(meta_df), "\n")
cat("Samples by HA vs other:\n")
print(table(meta_df$sample_id, meta_df$condition))
cat("Minimum cells per sample x macrophage subtype:", min_cells_per_sample, "\n")
cat("Gene filter: count >=", min_pseudobulk_count, "in at least smallest group size\n")
cat("Gene filter: detected in >=", min_detected_cell_fraction * 100, "% cells in at least smallest group size\n", sep = "")
cat("Remove MT genes:", remove_mt_genes, "\n")

macrophage_subtypes <- names(macrophage_subtype_colors)
macrophage_subtypes <- macrophage_subtypes[macrophage_subtypes %in% unique(meta_df$macrophage_subtype)]
summary_rows <- list()
tested_result_files <- character()

for (ct in macrophage_subtypes) {
  cat("\n------------------------------------------------------------\n")
  cat("Macrophage subtype:", ct, "\n")
  safe_ct <- sanitize_label(ct)

  cells_ct <- rownames(meta_df)[meta_df$macrophage_subtype == ct]
  meta_ct <- meta_df[cells_ct, , drop = FALSE]

  sample_cell_counts <- table(meta_ct$sample_id)
  keep_samples <- names(sample_cell_counts)[sample_cell_counts >= min_cells_per_sample]
  meta_ct <- meta_ct[meta_ct$sample_id %in% keep_samples, , drop = FALSE]

  coldata <- unique(meta_ct[, c("sample_id", "condition"), drop = FALSE])
  coldata <- coldata[order(coldata$condition, coldata$sample_id), , drop = FALSE]
  rownames(coldata) <- coldata$sample_id
  coldata$condition <- factor(as.character(coldata$condition), levels = group_levels)
  coldata$condition_deseq <- factor(
    ifelse(coldata$condition == "HA", "HA", "other"),
    levels = c("other", "HA")
  )
  coldata$n_cells_ct <- as.integer(sample_cell_counts[rownames(coldata)])
  coldata$median_nCount_RNA <- NA_real_
  coldata$median_nFeature_RNA <- NA_real_
  coldata$median_percent_mt <- NA_real_

  for (s in rownames(coldata)) {
    cells_sample <- rownames(meta_ct)[meta_ct$sample_id == s]
    meta_sample <- meta_ct[cells_sample, , drop = FALSE]
    coldata[s, "median_nCount_RNA"] <- median(meta_sample$nCount_RNA, na.rm = TRUE)
    coldata[s, "median_nFeature_RNA"] <- median(meta_sample$nFeature_RNA, na.rm = TRUE)
    coldata[s, "median_percent_mt"] <- median(meta_sample$percent.mt, na.rm = TRUE)
  }

  cat("Cells per retained sample:\n")
  print(table(meta_ct$sample_id, meta_ct$condition))
  cat("Sample-level pseudobulk metadata:\n")
  print(coldata)

  n_other <- sum(coldata$condition == "other")
  n_ha <- sum(coldata$condition == "HA")
  smallest_group_size <- min(n_other, n_ha)

  if (smallest_group_size < 2) {
    cat("Skipping: fewer than 2 retained samples in at least one group.\n")
    summary_rows[[ct]] <- data.frame(
      macrophage_subtype = ct,
      n_other_samples = n_other,
      n_ha_samples = n_ha,
      n_genes_tested = 0,
      n_sig_padj = 0,
      status = "skipped_low_samples",
      stringsAsFactors = FALSE
    )
    next
  }

  pb_counts <- matrix(
    0,
    nrow = nrow(rna_counts),
    ncol = nrow(coldata),
    dimnames = list(rownames(rna_counts), rownames(coldata))
  )

  detected_fraction <- matrix(
    0,
    nrow = nrow(rna_counts),
    ncol = nrow(coldata),
    dimnames = list(rownames(rna_counts), rownames(coldata))
  )

  for (s in rownames(coldata)) {
    cells_sample <- rownames(meta_ct)[meta_ct$sample_id == s]
    counts_sample <- rna_counts[, cells_sample, drop = FALSE]
    pb_counts[, s] <- rowSums(counts_sample)
    detected_fraction[, s] <- rowSums(counts_sample > 0) / length(cells_sample)
  }

  # keep genes with enough pseudobulk counts and enough detected cells
  keep_genes <- rowSums(pb_counts >= min_pseudobulk_count) >= smallest_group_size &
    rowSums(detected_fraction >= min_detected_cell_fraction) >= smallest_group_size

  if (remove_mt_genes) {
    keep_genes <- keep_genes & !grepl("^MT-", rownames(pb_counts))
  }

  cat("Genes before filtering:", nrow(pb_counts), "\n")
  cat("Genes after filtering :", sum(keep_genes), "\n")

  if (sum(keep_genes) == 0) {
    cat("Skipping: no genes passed pseudobulk filters.\n")
    summary_rows[[ct]] <- data.frame(
      macrophage_subtype = ct,
      n_other_samples = n_other,
      n_ha_samples = n_ha,
      n_genes_tested = 0,
      n_sig_padj = 0,
      status = "skipped_no_genes",
      stringsAsFactors = FALSE
    )
    next
  }

  pb_counts <- round(pb_counts[keep_genes, , drop = FALSE])

  dds <- DESeqDataSetFromMatrix(
    countData = pb_counts,
    colData = coldata,
    design = ~ condition_deseq
  )

  # inspect sample-level pseudobulk QC with rlog PCA
  rld <- rlog(dds, blind = TRUE)

  pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  pca_data$sample_id <- as.character(pca_data$name)
  pca_data$n_cells_ct <- coldata[pca_data$sample_id, "n_cells_ct"]
  pca_data$median_nCount_RNA <- coldata[pca_data$sample_id, "median_nCount_RNA"]
  pca_data$median_nFeature_RNA <- coldata[pca_data$sample_id, "median_nFeature_RNA"]
  pca_data$median_percent_mt <- coldata[pca_data$sample_id, "median_percent_mt"]

  p_pca_condition <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = sample_id)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.7, size = 3, color = "black", check_overlap = TRUE) +
    scale_color_manual(values = group_colors) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    coord_fixed() +
    ggtitle(paste0(ct, " | rlog PCA by HA status")) +
    theme_classic(base_size = 14)
  print(p_pca_condition)
  ggsave(
    filename = file.path(pca_condition_dir, paste0("pca_condition_", safe_ct, ".png")),
    plot = p_pca_condition,
    width = 6,
    height = 5,
    dpi = 600
  )

  p_pca_cells <- ggplot(pca_data, aes(PC1, PC2, color = n_cells_ct, label = sample_id)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.7, size = 3, color = "black", check_overlap = TRUE) +
    scale_color_gradient(low = "grey80", high = "#1F78B4") +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    coord_fixed() +
    ggtitle(paste0(ct, " | rlog PCA by n cells")) +
    theme_classic(base_size = 14)
  print(p_pca_cells)
  ggsave(
    filename = file.path(pca_n_cells_dir, paste0("pca_n_cells_", safe_ct, ".png")),
    plot = p_pca_cells,
    width = 6,
    height = 5,
    dpi = 600
  )

  p_pca_ncount <- ggplot(pca_data, aes(PC1, PC2, color = median_nCount_RNA, label = sample_id)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.7, size = 3, color = "black", check_overlap = TRUE) +
    scale_color_gradient(low = "grey80", high = "#1F78B4") +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    coord_fixed() +
    ggtitle(paste0(ct, " | rlog PCA by median nCount_RNA")) +
    theme_classic(base_size = 14)
  print(p_pca_ncount)
  ggsave(
    filename = file.path(pca_ncount_dir, paste0("pca_median_nCount_RNA_", safe_ct, ".png")),
    plot = p_pca_ncount,
    width = 6,
    height = 5,
    dpi = 600
  )

  p_pca_nfeature <- ggplot(pca_data, aes(PC1, PC2, color = median_nFeature_RNA, label = sample_id)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.7, size = 3, color = "black", check_overlap = TRUE) +
    scale_color_gradient(low = "grey80", high = "#1F78B4") +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    coord_fixed() +
    ggtitle(paste0(ct, " | rlog PCA by median nFeature_RNA")) +
    theme_classic(base_size = 14)
  print(p_pca_nfeature)
  ggsave(
    filename = file.path(pca_nfeature_dir, paste0("pca_median_nFeature_RNA_", safe_ct, ".png")),
    plot = p_pca_nfeature,
    width = 6,
    height = 5,
    dpi = 600
  )

  p_pca_mt <- ggplot(pca_data, aes(PC1, PC2, color = median_percent_mt, label = sample_id)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.7, size = 3, color = "black", check_overlap = TRUE) +
    scale_color_gradient(low = "grey80", high = "#D95F02") +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    coord_fixed() +
    ggtitle(paste0(ct, " | rlog PCA by median percent.mt")) +
    theme_classic(base_size = 14)
  print(p_pca_mt)
  ggsave(
    filename = file.path(pca_mt_dir, paste0("pca_median_percent_mt_", safe_ct, ".png")),
    plot = p_pca_mt,
    width = 6,
    height = 5,
    dpi = 600
  )

  dds <- DESeq(dds)
  cat("DESeq2 result names:\n")
  print(resultsNames(dds))

  res <- results(
    dds,
    contrast = c("condition_deseq", "HA", "other"),
    alpha = padj_thr
  )

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, res_df$pvalue), ]
  res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

  result_file <- file.path(res_dir, paste0("deseq2_HA_vs_other_", safe_ct, ".csv"))
  write.csv(res_df, result_file, row.names = FALSE)
  tested_result_files <- c(tested_result_files, result_file)

  sig_df <- res_df[
    !is.na(res_df$padj) &
      res_df$padj < padj_thr,
    ,
    drop = FALSE
  ]
  write.csv(sig_df, file.path(res_dir, paste0("deseq2_HA_vs_other_", safe_ct, "_significant.csv")), row.names = FALSE)

  cat("Top DESeq2 results:\n")
  print(head(res_df, 10))
  cat("Significant genes at padj < ", padj_thr, ": ", nrow(sig_df), "\n", sep = "")

  normalized_counts <- counts(dds, normalized = TRUE)
  write.csv(normalized_counts, file.path(res_dir, paste0("normalized_counts_", safe_ct, ".csv")))

  summary_rows[[ct]] <- data.frame(
    macrophage_subtype = ct,
    n_other_samples = n_other,
    n_ha_samples = n_ha,
    n_genes_tested = nrow(res_df),
    n_sig_padj = sum(!is.na(res_df$padj) & res_df$padj < padj_thr),
    status = "tested",
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(res_dir, "deseq2_pseudobulk_summary_by_macrophage_subtype.csv"), row.names = FALSE)

# create subtype-level volcano plots
cat("\nCreating volcano plots...\n")

for (res_file in tested_result_files) {
  res_df <- read.csv(res_file, stringsAsFactors = FALSE)
  safe_ct <- sub("^deseq2_HA_vs_other_", "", basename(res_file))
  safe_ct <- sub("\\.csv$", "", safe_ct)

  res_df$neg_log10_padj <- -log10(res_df$padj)
  res_df$volcano_group <- "Not significant"
  res_df$volcano_group[!is.na(res_df$padj) & res_df$padj < padj_thr & res_df$log2FoldChange > 0] <- "Up in HA"
  res_df$volcano_group[!is.na(res_df$padj) & res_df$padj < padj_thr & res_df$log2FoldChange < 0] <- "Up in other"
  res_df$volcano_group <- factor(res_df$volcano_group, levels = c("Up in other", "Not significant", "Up in HA"))

  plot_df <- res_df[
    !is.na(res_df$log2FoldChange) &
      !is.na(res_df$neg_log10_padj) &
      res_df$log2FoldChange >= -log2fc_plot_limit &
      res_df$log2FoldChange <= log2fc_plot_limit,
    ,
    drop = FALSE
  ]
  top_df <- plot_df[plot_df$volcano_group != "Not significant", , drop = FALSE]
  top_df <- top_df[order(top_df$padj), , drop = FALSE]
  top_df <- head(top_df, top_n_labels)

  p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = volcano_group)) +
    geom_point(size = 1.8, alpha = 0.82) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_text(
      data = top_df,
      aes(label = gene),
      color = "black",
      size = 4.2,
      vjust = -0.6,
      check_overlap = TRUE
    ) +
    scale_color_manual(
      values = c(
        "Up in other" = unname(group_colors["other"]),
        "Not significant" = "grey75",
        "Up in HA" = unname(group_colors["HA"])
      )
    ) +
    xlim(-log2fc_plot_limit, log2fc_plot_limit) +
    labs(
      title = safe_ct,
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = NULL
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 18, color = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.margin = margin(18, 18, 18, 18)
    )

  print(p_volcano)
  ggsave(
    filename = file.path(volcano_dir, paste0("volcano_HA_vs_other_", safe_ct, ".png")),
    plot = p_volcano,
    width = 9,
    height = 8,
    dpi = 600
  )
}

# create cumulative volcano plot across macrophage subtypes
cat("\nCreating cumulative volcano plot...\n")

cumulative_rows <- list()
for (res_file in tested_result_files) {
  res_df <- read.csv(res_file, stringsAsFactors = FALSE)
  safe_ct <- sub("^deseq2_HA_vs_other_", "", basename(res_file))
  safe_ct <- sub("\\.csv$", "", safe_ct)
  summary_safe_ct <- sanitize_label(summary_df$macrophage_subtype)
  cell_type_here <- summary_df$macrophage_subtype[summary_safe_ct == safe_ct][1]

  res_df$macrophage_subtype <- cell_type_here
  res_df$safe_cell_type <- safe_ct
  res_df$neg_log10_padj <- -log10(res_df$padj)
  res_df$is_significant <- !is.na(res_df$padj) & res_df$padj < padj_thr
  cumulative_rows[[safe_ct]] <- res_df
}

cumulative_df <- do.call(rbind, cumulative_rows)
cumulative_df <- cumulative_df[
  !is.na(cumulative_df$log2FoldChange) &
    !is.na(cumulative_df$neg_log10_padj) &
    cumulative_df$log2FoldChange >= -log2fc_plot_limit &
    cumulative_df$log2FoldChange <= log2fc_plot_limit,
  ,
  drop = FALSE
]

cumulative_df$plot_group <- "Not significant"
cumulative_df$plot_group[cumulative_df$is_significant] <- cumulative_df$macrophage_subtype[cumulative_df$is_significant]
cumulative_df$plot_group <- factor(cumulative_df$plot_group, levels = c("Not significant", names(macrophage_subtype_colors)))

cumulative_colors <- c("Not significant" = "grey75", macrophage_subtype_colors)
cumulative_colors <- cumulative_colors[levels(cumulative_df$plot_group)]
cumulative_colors <- cumulative_colors[!is.na(cumulative_colors)]

top_cumulative_df <- cumulative_df[cumulative_df$is_significant, , drop = FALSE]
top_cumulative_df <- top_cumulative_df[order(top_cumulative_df$padj), , drop = FALSE]
top_cumulative_df <- head(top_cumulative_df, top_n_labels)

p_cumulative_volcano <- ggplot(cumulative_df, aes(x = log2FoldChange, y = neg_log10_padj, color = plot_group)) +
  geom_point(size = 1.5, alpha = 0.72) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_text(
    data = top_cumulative_df,
    aes(label = gene),
    color = "black",
    size = 4.2,
    vjust = -0.6,
    check_overlap = TRUE
  ) +
  scale_color_manual(values = cumulative_colors, drop = TRUE) +
  xlim(-log2fc_plot_limit, log2fc_plot_limit) +
  labs(
    title = "Dysregulated genes in HA macrophage subtypes",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
    legend.text = element_text(size = 13),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(18, 18, 18, 18)
  )

print(p_cumulative_volcano)
ggsave(
  filename = file.path(volcano_dir, "volcano_cumulative_HA_vs_other_macrophage_subtypes.png"),
  plot = p_cumulative_volcano,
  width = 11,
  height = 8,
  dpi = 600
)

# save a summary of significant DEGs for each macrophage subtype
deg_summary_rows <- list()
for (ct in macrophage_subtypes) {
  safe_ct <- sanitize_label(ct)
  sig_file <- file.path(res_dir, paste0("deseq2_HA_vs_other_", safe_ct, "_significant.csv"))
  if (file.exists(sig_file)) {
    sig_df <- read.csv(sig_file, stringsAsFactors = FALSE)
    cat(ct, ": ", nrow(sig_df), " significant genes at padj < ", padj_thr, "\n", sep = "")
    deg_summary_rows[[ct]] <- data.frame(
      macrophage_subtype = ct,
      n_significant_degs = nrow(sig_df),
      significant_genes = paste(sig_df$gene, collapse = "; "),
      stringsAsFactors = FALSE
    )
  } else {
    cat(ct, ": no significant genes file found.\n")
    deg_summary_rows[[ct]] <- data.frame(
      macrophage_subtype = ct,
      n_significant_degs = 0,
      significant_genes = "",
      stringsAsFactors = FALSE
    )
  }
}

deg_summary_df <- do.call(rbind, deg_summary_rows)
write.csv(deg_summary_df, file.path(res_dir, "deseq2_significant_genes_summary_macrophage_subtypes.csv"), row.names = FALSE)

# bubble heatmap of significant iron-related DEGs across macrophage subtypes
cat("\nCreating iron-related DEG bubble heatmap...\n")

iron_genes <- load_iron_genes(iron_gene_files)
iron_bubble_rows <- list()

for (res_file in tested_result_files) {
  res_df <- read.csv(res_file, stringsAsFactors = FALSE)
  safe_ct <- sub("^deseq2_HA_vs_other_", "", basename(res_file))
  safe_ct <- sub("\\.csv$", "", safe_ct)
  summary_safe_ct <- sanitize_label(summary_df$macrophage_subtype)
  subtype_here <- summary_df$macrophage_subtype[summary_safe_ct == safe_ct][1]

  if (is.na(subtype_here) || !"gene" %in% colnames(res_df)) next

  res_df$gene_upper <- toupper(as.character(res_df$gene))
  res_df <- res_df[
    !is.na(res_df$padj) &
      res_df$padj < padj_thr &
      res_df$gene_upper %in% iron_genes,
    ,
    drop = FALSE
  ]

  if (nrow(res_df) == 0) next

  res_df$macrophage_subtype <- subtype_here
  iron_bubble_rows[[subtype_here]] <- res_df[, c("macrophage_subtype", "gene", "log2FoldChange", "padj"), drop = FALSE]
}

if (length(iron_bubble_rows) > 0) {
  iron_bubble_df <- do.call(rbind, iron_bubble_rows)
  rownames(iron_bubble_df) <- NULL
  iron_bubble_df$neg_log10_padj <- -log10(iron_bubble_df$padj)

  subtype_order <- summary_df$macrophage_subtype[summary_df$status == "tested"]
  subtype_order <- subtype_order[subtype_order %in% unique(iron_bubble_df$macrophage_subtype)]
  iron_bubble_df$macrophage_subtype <- factor(iron_bubble_df$macrophage_subtype, levels = subtype_order)

  gene_order_df <- aggregate(abs(log2FoldChange) ~ gene, data = iron_bubble_df, FUN = max)
  gene_order_df <- gene_order_df[order(gene_order_df$`abs(log2FoldChange)`, decreasing = TRUE), , drop = FALSE]
  iron_bubble_df$gene <- factor(iron_bubble_df$gene, levels = rev(gene_order_df$gene))

  max_abs_lfc <- max(abs(iron_bubble_df$log2FoldChange), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) {
    max_abs_lfc <- 1
  }

  p_iron_bubble <- ggplot(iron_bubble_df, aes(x = macrophage_subtype, y = gene)) +
    geom_point(aes(size = neg_log10_padj, fill = log2FoldChange), shape = 21, color = "black", stroke = 0.25) +
    scale_fill_gradient2(
      low = "#4C78A8",
      mid = "white",
      high = "#7A1F2B",
      midpoint = 0,
      limits = c(-max_abs_lfc, max_abs_lfc),
      name = "log2FC"
    ) +
    scale_size_continuous(
      name = "-log10 adj p",
      range = c(2.5, 9)
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.position = "top",
      legend.key.width = grid::unit(0.45, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(2.6, "cm"),
        barheight = grid::unit(0.35, "cm")
      ),
      size = guide_legend(title.position = "top")
    )

  ggsave(
    filename = iron_bubble_plot_file,
    plot = p_iron_bubble,
    width = 12,
    height = max(6, 0.28 * length(levels(iron_bubble_df$gene)) + 2),
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Macrophage subtype pseudobulk DESeq2 analysis complete.\n")
cat("Results dir : ", res_dir, "\n")
cat("Figures dir : ", fig_dir, "\n")
cat("Summary     : ", file.path(res_dir, "deseq2_pseudobulk_summary_by_macrophage_subtype.csv"), "\n")
cat("DEG summary : ", file.path(res_dir, "deseq2_significant_genes_summary_macrophage_subtypes.csv"), "\n")
if (length(iron_bubble_rows) > 0) {
  cat("Iron plot   : ", iron_bubble_plot_file, "\n")
}
cat("============================================================\n")
