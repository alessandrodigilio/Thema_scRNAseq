#########################
### utility functions ###
#########################

# -------------------------------------------------------------------------
# filtering helpers used by 01_filter_samples.R
# -------------------------------------------------------------------------

# filtering: threshold helpers
# these functions keep the threshold code readable in the main loop

# apply a lower QC threshold; missing thresholds are ignored
apply_min <- function(x, thr_val) {
  if (is.na(thr_val)) return(rep(TRUE, length(x)))
  x >= thr_val
}

# apply an upper QC threshold; missing thresholds are ignored
apply_max <- function(x, thr_val) {
  if (is.na(thr_val)) return(rep(TRUE, length(x)))
  x <= thr_val
}

# filtering: 10x input
# each sample is expected to have one raw feature-barcode folder

# read one 10x raw feature-barcode matrix from data/raw_counts
read_sample_counts <- function(sample_id, sample_path) {
  matrix_dir <- file.path(raw_data_dir, sample_path)

  if (!dir.exists(matrix_dir)) stop("Missing 10x input folder: ", matrix_dir)

  counts <- Read10X(data.dir = matrix_dir)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  cat("Loaded sample", sample_id, "from", matrix_dir, "\n")
  counts
}

# filtering: SoupX ambient RNA correction
# this block estimates and removes ambient RNA before the final Seurat object

# read the SoupX contamination estimate after autoEstCont
extract_soupx_rho <- function(soupx_obj) {
  rho <- soupx_obj$fit$rhoEst
  if (is.null(rho) || length(rho) == 0) return(NA_real_)
  as.numeric(rho[1])
}

# run SoupX on a filtered count matrix using the matching raw droplet matrix
run_soupx_correction <- function(sample_id, raw_counts, filtered_counts) {
  umi_before <- sum(filtered_counts)
  genes_before <- colSums(filtered_counts > 0)

  common_genes <- intersect(rownames(raw_counts), rownames(filtered_counts))
  if (length(common_genes) == 0) stop("No shared genes between raw and filtered matrices for ", sample_id)

  cat("Running SoupX for sample:", sample_id, "\n")

  filtered_counts <- filtered_counts[common_genes, , drop = FALSE]
  raw_counts <- raw_counts[common_genes, , drop = FALSE]

  soupx_obj <- SoupChannel(tod = raw_counts, toc = filtered_counts)

  # temporary clustering gives SoupX cell groups for contamination estimation
  tmp <- CreateSeuratObject(counts = filtered_counts, assay = "RNA", project = paste0(sample_id, "_SoupX"))
  tmp <- NormalizeData(tmp, verbose = FALSE)
  tmp <- FindVariableFeatures(tmp, verbose = FALSE)
  tmp <- ScaleData(tmp, verbose = FALSE)
  tmp <- RunPCA(tmp, npcs = max(10, scdblfinder_dims), verbose = FALSE)

  dims_use <- seq_len(min(scdblfinder_dims, ncol(Embeddings(tmp, "pca"))))
  tmp <- FindNeighbors(tmp, dims = dims_use, verbose = FALSE)
  tmp <- FindClusters(tmp, resolution = 0.2, verbose = FALSE)

  cluster_map <- setNames(as.character(tmp$seurat_clusters), colnames(tmp))
  soupx_obj <- setClusters(soupx_obj, cluster_map)
  soupx_obj <- autoEstCont(soupx_obj, doPlot = FALSE)

  corrected <- adjustCounts(soupx_obj, roundToInt = TRUE)
  corrected <- corrected[, colnames(filtered_counts), drop = FALSE]

  umi_after <- sum(corrected)
  genes_after <- colSums(corrected > 0)
  rho_est <- extract_soupx_rho(soupx_obj)
  removed_umis <- colSums(filtered_counts) - colSums(corrected)
  removed_fraction <- ifelse(
    colSums(filtered_counts) > 0,
    removed_umis / colSums(filtered_counts),
    NA_real_
  )

  cat(sprintf(
    "SoupX summary | sample=%s | rho=%.4f | umi_before=%d | umi_after=%d | umis_removed=%d | frac_removed=%.4f | median_genes_before=%.1f | median_genes_after=%.1f\n",
    sample_id,
    ifelse(is.na(rho_est), -1, rho_est),
    umi_before,
    umi_after,
    umi_before - umi_after,
    ifelse(umi_before > 0, (umi_before - umi_after) / umi_before, NA_real_),
    median(genes_before),
    median(genes_after)
  ))

  rm(tmp, soupx_obj)
  gc()

  list(
    counts = corrected,
    removed_umis = removed_umis,
    removed_fraction = removed_fraction,
    stats = list(
      soupx_status = "applied",
      soupx_rho = rho_est,
      soupx_umi_before = umi_before,
      soupx_umi_after = umi_after,
      soupx_umis_removed = umi_before - umi_after,
      soupx_fraction_removed = ifelse(umi_before > 0, (umi_before - umi_after) / umi_before, NA_real_),
      soupx_median_genes_before = median(genes_before),
      soupx_median_genes_after = median(genes_after)
    )
  )
}

# filtering: scDblFinder doublet removal
# this block labels cells, removes predicted doublets and saves a score plot

# run scDblFinder and keep only singlet cells
run_scdblfinder_filter <- function(obj, sample_id) {
  cells_before <- ncol(obj)

  cat("Running scDblFinder for sample:", sample_id, "\n")

  sce <- SingleCellExperiment(
    assays = list(counts = GetAssayData(obj, assay = "RNA", layer = "counts"))
  )

  sce <- scDblFinder(
    sce,
    dims = scdblfinder_dims,
    verbose = FALSE
  )

  dbl_class <- as.character(colData(sce)$scDblFinder.class)
  dbl_score <- as.numeric(colData(sce)$scDblFinder.score)

  obj$scDblFinder.class <- dbl_class
  obj$scDblFinder.score <- dbl_score

  set.seed(10010101)
  sce <- runPCA(sce, ncomponents = scdblfinder_dims, exprs_values = "counts")
  sce <- runTSNE(sce, dimred = "PCA")
  plot_score <- plotTSNE(sce, colour_by = "scDblFinder.score") +
    ggtitle(paste0(sample_id, " scDblFinder score"))

  singlets <- colnames(sce)[dbl_class == "singlet"]
  n_singlets <- length(singlets)
  n_doublets_removed <- cells_before - n_singlets

  cat(sprintf(
    "scDblFinder summary | sample=%s | cells_in=%d | singlets=%d | doublets_removed=%d\n",
    sample_id,
    cells_before,
    n_singlets,
    n_doublets_removed
  ))

  obj <- subset(obj, cells = singlets)

  list(
    obj = obj,
    plot_score = plot_score,
    stats = list(
      scdblfinder_status = "applied",
      scdblfinder_singlets = n_singlets,
      scdblfinder_doublets_removed = n_doublets_removed
    )
  )
}

# filtering: QC plot
# this panel summarizes ambient RNA removal and doublet scores for each sample

# save one QC panel with SoupX removal and scDblFinder score
save_filtering_panel <- function(obj_qc, plot_score, sample_id) {
  soupx_qc_fraction <- obj_qc$soupx_removed_fraction

  plot_soupx <- ggplot(
    data.frame(soupx_removed_fraction = soupx_qc_fraction),
    aes(x = soupx_removed_fraction)
  ) +
    geom_histogram(bins = 40, fill = "grey70", color = "white") +
    labs(
      title = paste0(sample_id, " SoupX removed RNA fraction"),
      x = "Removed fraction",
      y = "Cells"
    ) +
    theme_classic()

  panel_plot <- wrap_plots(plot_soupx, plot_score, ncol = 2)

  panel_path <- file.path(fig_filtering_dir, paste0(sample_id, "_soupx_scdblfinder_panel.png"))
  ggsave(
    filename = panel_path,
    plot = panel_plot,
    width = 12,
    height = 5,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  panel_path
}

# -------------------------------------------------------------------------
# integration helpers used by 02_integration.R
# -------------------------------------------------------------------------

# convert numeric clustering resolutions into stable metadata column names
res_to_colname <- function(resolution) {
  paste0("snn_res.", format(resolution, nsmall = 1, trim = TRUE))
}

# -------------------------------------------------------------------------
# manual annotation helpers used by 04_annotation.R
# -------------------------------------------------------------------------

# sort cluster labels numerically when possible
sort_cluster_levels <- function(x) {
  x <- unique(as.character(x))
  suppressWarnings(x_num <- as.integer(x))
  if (all(!is.na(x_num))) return(as.character(sort(x_num)))
  sort(x)
}

# get the colors used by the cell types present in the object
build_celltype_colors <- function(labels, color_map) {
  labels <- unique(as.character(labels))
  color_map[labels]
}

# convert a marker list into a simple cell type/gene table
build_marker_table <- function(marker_list, genes_available) {
  marker_rows <- list()

  for (cell_type in names(marker_list)) {
    genes_here <- unique(marker_list[[cell_type]])
    genes_here <- genes_here[genes_here %in% genes_available]
    marker_rows[[cell_type]] <- data.frame(cell_type = cell_type, gene = genes_here, stringsAsFactors = FALSE)
  }

  do.call(rbind, marker_rows)
}

# count cell type ratios inside each sample or condition
build_ratio_plot_data <- function(meta_df, x_col, fill_col) {
  plot_df <- meta_df[, c(x_col, fill_col), drop = FALSE]
  plot_df[[x_col]] <- as.character(plot_df[[x_col]])
  plot_df[[fill_col]] <- as.character(plot_df[[fill_col]])

  count_df <- as.data.frame(table(plot_df[[x_col]], plot_df[[fill_col]]), stringsAsFactors = FALSE)
  colnames(count_df) <- c("group", "cell_type", "n_cells")
  count_df <- count_df[count_df$n_cells > 0, , drop = FALSE]

  totals <- aggregate(n_cells ~ group, data = count_df, FUN = sum)
  count_df <- merge(count_df, totals, by = "group", suffixes = c("", "_total"), sort = FALSE)
  count_df$ratio <- count_df$n_cells / count_df$n_cells_total
  count_df$group <- factor(count_df$group, levels = unique(count_df$group))
  count_df$cell_type <- factor(count_df$cell_type, levels = unique(plot_df[[fill_col]]))
  count_df
}

# make a stacked barplot of cell type ratios
make_ratio_plot <- function(plot_df, fill_colors) {
  ggplot(plot_df, aes(x = group, y = ratio, fill = cell_type)) +
    geom_col(width = 0.92, color = NA) +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Ratio", fill = NULL) +
    theme_classic(base_size = 24) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 24, color = "black"),
      axis.text.y = element_text(size = 24, color = "black"),
      axis.title.y = element_text(size = 24, color = "black"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length = grid::unit(0.22, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.position = "right",
      legend.key.size = grid::unit(0.45, "cm"),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# -------------------------------------------------------------------------
# composition helpers used by 05_celltype_composition_ha_vs_other.R
# -------------------------------------------------------------------------

# calculate cell type fractions inside each sample
build_sample_composition <- function(meta_df) {
  sample_condition <- unique(meta_df[, c("sample_id", "condition"), drop = FALSE])

  composition_df <- as.data.frame(table(meta_df$sample_id, meta_df$cell_type), stringsAsFactors = FALSE)
  colnames(composition_df) <- c("sample_id", "cell_type", "n_cells")

  composition_df <- merge(composition_df, sample_condition, by = "sample_id", sort = FALSE)
  sample_totals <- aggregate(n_cells ~ sample_id, data = composition_df, FUN = sum)
  colnames(sample_totals)[2] <- "sample_total_cells"

  composition_df <- merge(composition_df, sample_totals, by = "sample_id", sort = FALSE)
  composition_df$fraction <- composition_df$n_cells / composition_df$sample_total_cells
  composition_df
}

# test HA versus other cell type fractions across samples
test_composition_ha_vs_other <- function(composition_df, cell_types, group_levels, pseudocount) {
  stat_rows <- list()

  for (cell_type in cell_types) {
    df_ct <- composition_df[composition_df$cell_type == cell_type, , drop = FALSE]
    other_vals <- df_ct$fraction[df_ct$condition == group_levels[1]]
    ha_vals <- df_ct$fraction[df_ct$condition == group_levels[2]]

    mean_other <- mean(other_vals)
    mean_ha <- mean(ha_vals)

    stat_rows[[cell_type]] <- data.frame(
      cell_type = cell_type,
      n_other_samples = length(other_vals),
      n_ha_samples = length(ha_vals),
      mean_other = mean_other,
      mean_ha = mean_ha,
      median_other = median(other_vals),
      median_ha = median(ha_vals),
      diff_mean_ha_minus_other = mean_ha - mean_other,
      log2fc_mean_ha_vs_other = log2((mean_ha + pseudocount) / (mean_other + pseudocount)),
      p_wilcox = wilcox.test(ha_vals, other_vals, exact = FALSE)$p.value,
      stringsAsFactors = FALSE
    )
  }

  stat_df <- do.call(rbind, stat_rows)
  stat_df$padj_bh <- p.adjust(stat_df$p_wilcox, method = "BH")
  stat_df[order(stat_df$p_wilcox), , drop = FALSE]
}

# -------------------------------------------------------------------------
# pseudobulk helpers used by 06_pseudobulk_deseq2_ha_vs_other.R
# -------------------------------------------------------------------------

# make a compact file-safe label
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# build the sample-level metadata used by DESeq2
make_pseudobulk_coldata <- function(meta_ct, group_levels) {
  sample_cell_counts <- table(meta_ct$sample_id)

  coldata <- unique(meta_ct[, c("sample_id", "condition"), drop = FALSE])
  coldata <- coldata[order(coldata$condition, coldata$sample_id), , drop = FALSE]
  rownames(coldata) <- coldata$sample_id
  coldata$condition <- factor(as.character(coldata$condition), levels = group_levels)
  coldata$condition_deseq <- factor(ifelse(coldata$condition == "HA", "HA", "other"), levels = c("other", "HA"))
  coldata$n_cells_ct <- as.integer(sample_cell_counts[rownames(coldata)])

  coldata$median_nCount_RNA <- NA_real_
  coldata$median_nFeature_RNA <- NA_real_
  coldata$median_percent_mt <- NA_real_

  for (sample_id in rownames(coldata)) {
    cells_sample <- rownames(meta_ct)[meta_ct$sample_id == sample_id]
    meta_sample <- meta_ct[cells_sample, , drop = FALSE]
    coldata[sample_id, "median_nCount_RNA"] <- median(meta_sample$nCount_RNA, na.rm = TRUE)
    coldata[sample_id, "median_nFeature_RNA"] <- median(meta_sample$nFeature_RNA, na.rm = TRUE)
    coldata[sample_id, "median_percent_mt"] <- median(meta_sample$percent.mt, na.rm = TRUE)
  }

  coldata
}

# sum raw counts inside each sample and filter genes before DESeq2
make_pseudobulk_counts <- function(rna_counts,
                                   meta_ct,
                                   coldata,
                                   smallest_group_size,
                                   min_pseudobulk_count,
                                   min_detected_cell_fraction,
                                   remove_mt_genes = FALSE) {
  pb_counts <- matrix(0, nrow = nrow(rna_counts), ncol = nrow(coldata), dimnames = list(rownames(rna_counts), rownames(coldata)))
  detected_fraction <- pb_counts

  for (sample_id in rownames(coldata)) {
    cells_sample <- rownames(meta_ct)[meta_ct$sample_id == sample_id]
    counts_sample <- rna_counts[, cells_sample, drop = FALSE]
    pb_counts[, sample_id] <- rowSums(counts_sample)
    detected_fraction[, sample_id] <- rowSums(counts_sample > 0) / length(cells_sample)
  }

  keep_genes <- rowSums(pb_counts >= min_pseudobulk_count) >= smallest_group_size &
    rowSums(detected_fraction >= min_detected_cell_fraction) >= smallest_group_size

  if (remove_mt_genes) keep_genes <- keep_genes & !grepl("^MT-", rownames(pb_counts))

  list(
    counts = round(pb_counts[keep_genes, , drop = FALSE]),
    n_genes_before = nrow(pb_counts),
    n_genes_after = sum(keep_genes)
  )
}

# run pseudobulk DESeq2 for one cell type
run_pseudobulk_deseq2_celltype <- function(cell_type,
                                           rna_counts,
                                           meta_df,
                                           res_dir,
                                           group_levels,
                                           group_colors,
                                           min_cells_per_sample,
                                           min_pseudobulk_count,
                                           min_detected_cell_fraction,
                                           padj_thr,
                                           remove_mt_genes = FALSE) {
  safe_ct <- safe_name(cell_type)
  meta_ct <- meta_df[meta_df$cell_type == cell_type, , drop = FALSE]

  sample_cell_counts <- table(meta_ct$sample_id)
  keep_samples <- names(sample_cell_counts)[sample_cell_counts >= min_cells_per_sample]
  meta_ct <- meta_ct[meta_ct$sample_id %in% keep_samples, , drop = FALSE]
  coldata <- make_pseudobulk_coldata(meta_ct, group_levels)

  n_other <- sum(coldata$condition == "other")
  n_ha <- sum(coldata$condition == "HA")
  smallest_group_size <- min(n_other, n_ha)

  if (smallest_group_size < 2) {
    return(list(
      result_file = NA_character_,
      summary = data.frame(cell_type = cell_type, n_other_samples = n_other, n_ha_samples = n_ha, n_genes_tested = 0, n_sig_padj = 0, status = "skipped_low_samples", stringsAsFactors = FALSE)
    ))
  }

  pb <- make_pseudobulk_counts(
    rna_counts = rna_counts,
    meta_ct = meta_ct,
    coldata = coldata,
    smallest_group_size = smallest_group_size,
    min_pseudobulk_count = min_pseudobulk_count,
    min_detected_cell_fraction = min_detected_cell_fraction,
    remove_mt_genes = remove_mt_genes
  )

  if (pb$n_genes_after == 0) {
    return(list(
      result_file = NA_character_,
      summary = data.frame(cell_type = cell_type, n_other_samples = n_other, n_ha_samples = n_ha, n_genes_tested = 0, n_sig_padj = 0, status = "skipped_no_genes", stringsAsFactors = FALSE)
    ))
  }

  dds <- DESeqDataSetFromMatrix(countData = pb$counts, colData = coldata, design = ~ condition_deseq)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition_deseq", "HA", "other"), alpha = padj_thr)

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]
  res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

  result_file <- file.path(res_dir, paste0("deseq2_HA_vs_other_", safe_ct, ".csv"))
  write.csv(res_df, result_file, row.names = FALSE)

  sig_df <- res_df[!is.na(res_df$padj) & res_df$padj < padj_thr, , drop = FALSE]
  write.csv(sig_df, file.path(res_dir, paste0("deseq2_HA_vs_other_", safe_ct, "_significant.csv")), row.names = FALSE)
  write.csv(counts(dds, normalized = TRUE), file.path(res_dir, paste0("normalized_counts_", safe_ct, ".csv")))

  list(
    result_file = result_file,
    summary = data.frame(
      cell_type = cell_type,
      n_other_samples = n_other,
      n_ha_samples = n_ha,
      n_genes_tested = nrow(res_df),
      n_sig_padj = nrow(sig_df),
      status = "tested",
      stringsAsFactors = FALSE
    )
  )
}

# make a volcano plot from one DESeq2 result table
make_pseudobulk_volcano <- function(res_df, title, group_colors, padj_thr, top_n_labels, log2fc_plot_limit) {
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
  top_df <- head(top_df[order(top_df$padj), , drop = FALSE], top_n_labels)

  ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = volcano_group)) +
    geom_point(size = 1.8, alpha = 0.82) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_text(data = top_df, aes(label = gene), color = "black", size = 4.2, vjust = -0.6, check_overlap = TRUE) +
    scale_color_manual(values = c("Up in other" = unname(group_colors["other"]), "Not significant" = "grey75", "Up in HA" = unname(group_colors["HA"]))) +
    xlim(-log2fc_plot_limit, log2fc_plot_limit) +
    labs(title = title, x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
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
}

# save individual volcano plots
save_pseudobulk_volcanoes <- function(result_files, volcano_dir, group_colors, padj_thr, top_n_labels, log2fc_plot_limit) {
  for (result_file in result_files) {
    res_df <- read.csv(result_file, stringsAsFactors = FALSE)
    safe_ct <- sub("^deseq2_HA_vs_other_", "", basename(result_file))
    safe_ct <- sub("\\.csv$", "", safe_ct)

    p <- make_pseudobulk_volcano(res_df, safe_ct, group_colors, padj_thr, top_n_labels, log2fc_plot_limit)
    ggsave(file.path(volcano_dir, paste0("volcano_HA_vs_other_", safe_ct, ".png")), p, width = 9, height = 8, dpi = 600)
  }
}

# save the cumulative volcano plot across all tested cell types
save_cumulative_pseudobulk_volcano <- function(result_files,
                                               summary_df,
                                               volcano_dir,
                                               cluster_name_colors,
                                               padj_thr,
                                               top_n_labels,
                                               log2fc_plot_limit) {
  cumulative_rows <- list()
  summary_safe_ct <- safe_name(summary_df$cell_type)

  for (result_file in result_files) {
    res_df <- read.csv(result_file, stringsAsFactors = FALSE)
    safe_ct <- sub("^deseq2_HA_vs_other_", "", basename(result_file))
    safe_ct <- sub("\\.csv$", "", safe_ct)

    res_df$cell_type <- summary_df$cell_type[summary_safe_ct == safe_ct][1]
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
  cumulative_df$plot_group[cumulative_df$is_significant] <- cumulative_df$cell_type[cumulative_df$is_significant]
  cumulative_df$plot_group <- factor(cumulative_df$plot_group, levels = c("Not significant", names(cluster_name_colors)))

  cumulative_colors <- c("Not significant" = "grey75", cluster_name_colors)
  cumulative_colors <- cumulative_colors[levels(cumulative_df$plot_group)]
  cumulative_colors <- cumulative_colors[!is.na(cumulative_colors)]

  top_df <- cumulative_df[cumulative_df$is_significant, , drop = FALSE]
  top_df <- head(top_df[order(top_df$padj), , drop = FALSE], top_n_labels)

  p <- ggplot(cumulative_df, aes(x = log2FoldChange, y = neg_log10_padj, color = plot_group)) +
    geom_point(size = 1.5, alpha = 0.72) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_text(data = top_df, aes(label = gene), color = "black", size = 4.2, vjust = -0.6, check_overlap = TRUE) +
    scale_color_manual(values = cumulative_colors, drop = TRUE) +
    xlim(-log2fc_plot_limit, log2fc_plot_limit) +
    labs(title = "Dysregulated genes in HA synovium", x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
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

  ggsave(file.path(volcano_dir, "volcano_cumulative_HA_vs_other.png"), p, width = 11, height = 8, dpi = 600)
}

# summarize significant genes from saved DESeq2 result tables
summarize_pseudobulk_significant_genes <- function(summary_df, res_dir) {
  deg_summary_rows <- list()

  for (cell_type in summary_df$cell_type) {
    safe_ct <- safe_name(cell_type)
    sig_file <- file.path(res_dir, paste0("deseq2_HA_vs_other_", safe_ct, "_significant.csv"))

    if (file.exists(sig_file)) {
      sig_df <- read.csv(sig_file, stringsAsFactors = FALSE)
      deg_summary_rows[[cell_type]] <- data.frame(cell_type = cell_type, n_significant_degs = nrow(sig_df), significant_genes = paste(sig_df$gene, collapse = "; "), stringsAsFactors = FALSE)
    } else {
      deg_summary_rows[[cell_type]] <- data.frame(cell_type = cell_type, n_significant_degs = 0, significant_genes = "", stringsAsFactors = FALSE)
    }
  }

do.call(rbind, deg_summary_rows)
}

# -------------------------------------------------------------------------
# iron and ferroptosis helpers used by 08_ferroptosis_ha_vs_other.R
# -------------------------------------------------------------------------

# read gene symbols from one or more Excel files
load_excel_gene_list <- function(files, include_column_names = TRUE) {
  genes <- character(0)

  for (f in files) {
    x <- read.xlsx(f)
    if (ncol(x) == 0) next

    genes_here <- x[[1]]
    if (include_column_names) genes_here <- c(colnames(x)[1], genes_here)

    genes <- c(genes, genes_here)
  }

  genes <- toupper(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & genes != "" & genes != "GENE"]
  genes[genes == "TRFC"] <- "TFRC"
  sort(unique(genes))
}

# make sure the RNA assay has normalized data before module scoring
prepare_rna_assay_for_scoring <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj, assay = "RNA")

  rna_data <- tryCatch(LayerData(obj[["RNA"]], layer = "data"), error = function(e) NULL)
  has_data <- !is.null(rna_data)

  if (has_data && inherits(rna_data, "dgCMatrix")) {
    has_data <- length(rna_data@x) > 0
  }

  if (!has_data) {
    cat("Normalizing RNA assay...\n")
    obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }

  obj
}

# add the curated ferroptosis module score to cell metadata
add_ferroptosis_score <- function(obj, genes_use, score_col = "ferroptosis_score") {
  obj <- AddModuleScore(
    object = obj,
    features = list(genes_use),
    assay = "RNA",
    name = "Ferroptosis_Score",
    ctrl = 100,
    seed = 1234
  )

  obj[[score_col]] <- obj$Ferroptosis_Score1
  obj
}

# summarize one score by any set of metadata columns
summarize_score_by <- function(meta_df, score_col, group_cols) {
  formula_txt <- paste(score_col, "~", paste(group_cols, collapse = " + "))
  score_summary <- aggregate(
    as.formula(formula_txt),
    data = meta_df,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = length(x))
  )

  out <- score_summary[, group_cols, drop = FALSE]
  score_mat <- score_summary[[score_col]]
  out$mean_score <- score_mat[, "mean"]
  out$median_score <- score_mat[, "median"]
  out$sd_score <- score_mat[, "sd"]
  out$n_cells <- score_mat[, "n"]
  out
}

# draw a UMAP colored by ferroptosis score
make_score_umap <- function(obj, score_col, reduction_name, score_colors, score_max, split_by = NULL) {
  FeaturePlot(
    object = obj,
    features = score_col,
    reduction = reduction_name,
    split.by = split_by,
    raster = FALSE,
    order = TRUE,
    cols = score_colors,
    min.cutoff = 0,
    max.cutoff = score_max
  )
}

# draw ferroptosis score violins by HA/other
make_score_violin_group <- function(meta_df, score_col, group_col, fill_colors, y_label) {
  ggplot(meta_df, aes(x = .data[[group_col]], y = .data[[score_col]], fill = .data[[group_col]])) +
    geom_violin(trim = TRUE, scale = "width", color = NA) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.4) +
    scale_fill_manual(values = fill_colors) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      legend.position = "none"
    )
}

# draw ferroptosis score violins by cell type and condition
make_score_violin_celltype <- function(meta_df, score_col, celltype_col, group_col, fill_colors, y_label) {
  ggplot(meta_df, aes(x = .data[[celltype_col]], y = .data[[score_col]], fill = .data[[group_col]])) +
    geom_violin(trim = TRUE, scale = "width", position = position_dodge(width = 0.8), color = NA) +
    scale_fill_manual(values = fill_colors) +
    labs(x = NULL, y = y_label, fill = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      legend.position = "top"
    )
}

# put each word on a separate line in facet labels
wrap_each_word <- function(x) {
  x <- gsub(" / ", "/", as.character(x), fixed = TRUE)
  vapply(strsplit(x, " "), function(words) paste(words, collapse = "\n"), character(1))
}

# draw average ferroptosis gene expression by cell type and condition
make_ferroptosis_gene_heatmap <- function(obj, genes_use, heatmap_celltypes) {
  dot_features <- genes_use[genes_use %in% rownames(obj)]
  heatmap_celltypes <- intersect(heatmap_celltypes, levels(obj$cell_type))

  if (length(dot_features) == 0 || length(heatmap_celltypes) == 0) return(NULL)

  keep_cells <- rownames(obj@meta.data)[as.character(obj$cell_type) %in% heatmap_celltypes]
  obj_heatmap <- subset(obj, cells = keep_cells)
  obj_heatmap$cell_type <- factor(obj_heatmap$cell_type, levels = heatmap_celltypes)
  obj_heatmap$celltype_group <- paste(obj_heatmap$cell_type, obj_heatmap$HA_vs_other, sep = " | ")

  group_levels <- as.vector(rbind(
    paste(heatmap_celltypes, "HA", sep = " | "),
    paste(heatmap_celltypes, "other", sep = " | ")
  ))

  obj_heatmap$celltype_group <- factor(obj_heatmap$celltype_group, levels = group_levels)

  avg_expr <- AverageExpression(
    object = obj_heatmap,
    assays = "RNA",
    features = dot_features,
    group.by = "celltype_group",
    slot = "data",
    verbose = FALSE
  )$RNA

  avg_expr <- avg_expr[, group_levels[group_levels %in% colnames(avg_expr)], drop = FALSE]
  if (ncol(avg_expr) == 0) return(NULL)

  scaled_expr <- t(scale(t(as.matrix(avg_expr))))
  scaled_expr[is.na(scaled_expr)] <- 0

  heatmap_df <- as.data.frame(as.table(scaled_expr), stringsAsFactors = FALSE)
  colnames(heatmap_df) <- c("gene", "group", "z_score")
  heatmap_df$gene <- factor(heatmap_df$gene, levels = rev(dot_features))
  heatmap_df$group <- factor(heatmap_df$group, levels = colnames(scaled_expr))

  heatmap_groups <- strsplit(as.character(heatmap_df$group), " | ", fixed = TRUE)
  heatmap_df$cell_type <- factor(vapply(heatmap_groups, `[`, character(1), 1), levels = heatmap_celltypes)
  heatmap_df$HA_vs_other <- factor(vapply(heatmap_groups, `[`, character(1), 2), levels = c("HA", "other"))

  p <- ggplot(heatmap_df, aes(x = HA_vs_other, y = gene, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(cols = vars(cell_type), scales = "free_x", space = "free_x", labeller = labeller(cell_type = wrap_each_word)) +
    scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#7A1F2B", midpoint = 0, name = "Scaled\nexpression") +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.spacing.x = grid::unit(0.22, "cm"),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 8.5, color = "black", lineheight = 0.95),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11),
      legend.position = "top",
      legend.key.width = grid::unit(0.5, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.08, "cm"),
      plot.margin = margin(20, 35, 35, 25)
    )

  list(plot = p, n_celltypes = length(unique(heatmap_df$cell_type)), n_genes = nrow(scaled_expr))
}

# collect iron-related genes found in each atlas pseudobulk DESeq2 result table
scan_iron_genes_in_pseudobulk <- function(summary_df, pseudobulk_res_dir, iron_genes, padj_thr) {
  all_hits <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(summary_df))) {
    cell_type <- summary_df$cell_type[i]
    safe_ct <- safe_name(cell_type)
    deg_file <- file.path(pseudobulk_res_dir, paste0("deseq2_HA_vs_other_", safe_ct, ".csv"))
    if (!file.exists(deg_file)) next

    deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)
    deg_df$gene_upper <- toupper(as.character(deg_df$gene))
    deg_df <- deg_df[deg_df$gene_upper %in% iron_genes, , drop = FALSE]
    if (nrow(deg_df) == 0) next

    deg_df$cell_type <- cell_type
    deg_df$direction_in_HA <- ifelse(deg_df$log2FoldChange > 0, "Up in HA", ifelse(deg_df$log2FoldChange < 0, "Down in HA", "No change"))
    deg_df$is_significant <- !is.na(deg_df$padj) & deg_df$padj < padj_thr

    out <- deg_df[, c("cell_type", "gene", "log2FoldChange", "padj", "direction_in_HA", "is_significant"), drop = FALSE]
    out <- out[order(out$padj, -abs(out$log2FoldChange)), , drop = FALSE]

    all_hits[[cell_type]] <- out
    summary_rows[[cell_type]] <- data.frame(cell_type = cell_type, n_iron_genes_found = nrow(out), n_significant_padj_lt_0.05 = sum(out$is_significant, na.rm = TRUE), stringsAsFactors = FALSE)
  }

  if (length(all_hits) == 0) {
    empty_hits <- data.frame(cell_type = character(0), gene = character(0), log2FoldChange = numeric(0), padj = numeric(0), direction_in_HA = character(0), is_significant = logical(0))
    empty_summary <- data.frame(cell_type = character(0), n_iron_genes_found = integer(0), n_significant_padj_lt_0.05 = integer(0))
    return(list(all_hits = empty_hits, significant_hits = empty_hits, summary = empty_summary, by_celltype = all_hits))
  }

  all_hits_df <- do.call(rbind, all_hits)
  summary_df_out <- do.call(rbind, summary_rows)

  rownames(all_hits_df) <- NULL
  rownames(summary_df_out) <- NULL

  list(
    all_hits = all_hits_df,
    significant_hits = all_hits_df[all_hits_df$is_significant, , drop = FALSE],
    summary = summary_df_out,
    by_celltype = all_hits
  )
}

# save iron-related DEGs to one workbook
save_iron_gene_workbook <- function(iron_results, output_xlsx) {
  wb <- createWorkbook()

  addWorksheet(wb, "significant_padj_lt_0.05")
  writeData(wb, "significant_padj_lt_0.05", iron_results$significant_hits)

  addWorksheet(wb, "summary")
  writeData(wb, "summary", iron_results$summary)

  for (nm in names(iron_results$by_celltype)) {
    df_here <- iron_results$by_celltype[[nm]]
    df_here <- df_here[df_here$is_significant, , drop = FALSE]
    if (nrow(df_here) == 0) next

    sheet_name <- substr(safe_name(nm), 1, 31)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df_here)
  }

  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
}

# draw a bubble heatmap of significant iron-related pseudobulk DEGs
make_iron_gene_bubble_plot <- function(sig_hits_df, celltype_order) {
  if (nrow(sig_hits_df) == 0) return(NULL)

  plot_df <- sig_hits_df
  plot_df$neg_log10_padj <- -log10(plot_df$padj)

  celltype_order <- celltype_order[celltype_order %in% unique(plot_df$cell_type)]
  plot_df$cell_type <- factor(plot_df$cell_type, levels = celltype_order)

  gene_order_df <- aggregate(abs(log2FoldChange) ~ gene, data = plot_df, FUN = max)
  gene_order_df <- gene_order_df[order(gene_order_df$`abs(log2FoldChange)`, decreasing = TRUE), , drop = FALSE]
  plot_df$gene <- factor(plot_df$gene, levels = rev(gene_order_df$gene))

  max_abs_lfc <- max(abs(plot_df$log2FoldChange), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) max_abs_lfc <- 1

  ggplot(plot_df, aes(x = cell_type, y = gene)) +
    geom_point(aes(size = neg_log10_padj, fill = log2FoldChange), shape = 21, color = "black", stroke = 0.25) +
    scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#7A1F2B", midpoint = 0, limits = c(-max_abs_lfc, max_abs_lfc), name = "log2FC") +
    scale_size_continuous(name = "-log10 adj p", range = c(2.5, 9)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title = element_blank(),
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
      fill = guide_colorbar(title.position = "top", barwidth = grid::unit(2.6, "cm"), barheight = grid::unit(0.35, "cm")),
      size = guide_legend(title.position = "top")
    )
}

# -------------------------------------------------------------------------
# enrichment helpers used by 07_enrichment_ha_vs_other.R
# -------------------------------------------------------------------------

# load one MSigDB table using either the new or old msigdbr argument names
load_msigdb_table <- function(species, collection, subcollection) {
  out <- tryCatch(
    msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) > 0) return(as.data.frame(out))

  out <- tryCatch(
    msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) > 0) return(as.data.frame(out))
  NULL
}

# build GO BP, KEGG and Reactome pathway lists for fgsea
load_msigdb_pathways <- function(msigdb_species) {
  pathway_tables <- list()

  gobp_df <- load_msigdb_table(msigdb_species, "C5", "GO:BP")
  gobp_df$database <- "GO_BP"
  pathway_tables[["GO_BP"]] <- gobp_df

  kegg_tables <- list(
    load_msigdb_table(msigdb_species, "C2", "CP:KEGG"),
    load_msigdb_table(msigdb_species, "C2", "CP:KEGG_LEGACY"),
    load_msigdb_table(msigdb_species, "C2", "CP:KEGG_MEDICUS")
  )
  kegg_tables <- kegg_tables[!vapply(kegg_tables, is.null, logical(1))]
  kegg_df <- do.call(rbind, kegg_tables)
  kegg_df$database <- "KEGG"
  pathway_tables[["KEGG"]] <- kegg_df

  reactome_df <- load_msigdb_table(msigdb_species, "C2", "CP:REACTOME")
  reactome_df$database <- "Reactome"
  pathway_tables[["Reactome"]] <- reactome_df

  pathway_df <- do.call(rbind, pathway_tables)
  pathway_df$gene_symbol <- toupper(as.character(pathway_df$gene_symbol))
  pathway_df$pathway_id <- paste(pathway_df$database, pathway_df$gs_name, sep = "__")
  pathway_df <- pathway_df[!is.na(pathway_df$gene_symbol) & pathway_df$gene_symbol != "", , drop = FALSE]

  pathways <- split(pathway_df$gene_symbol, pathway_df$pathway_id)
  pathways <- lapply(pathways, function(x) sort(unique(x)))
  pathways <- pathways[vapply(pathways, length, integer(1)) > 0]

  pathway_info <- unique(pathway_df[, c("pathway_id", "database", "gs_name"), drop = FALSE])
  rownames(pathway_info) <- pathway_info$pathway_id

  list(pathways = pathways, pathway_info = pathway_info)
}

# build signed ranks from DESeq2 results
build_rank_vector <- function(deg_df) {
  if ("stat" %in% colnames(deg_df)) {
    rank_df <- deg_df[!is.na(deg_df$stat), c("gene", "stat"), drop = FALSE]
    colnames(rank_df) <- c("gene", "rank_value")
  } else {
    rank_df <- deg_df[!is.na(deg_df$log2FoldChange) & !is.na(deg_df$pvalue) & deg_df$pvalue > 0, c("gene", "log2FoldChange", "pvalue"), drop = FALSE]
    rank_df$rank_value <- sign(rank_df$log2FoldChange) * -log10(rank_df$pvalue)
    rank_df <- rank_df[, c("gene", "rank_value"), drop = FALSE]
  }

  rank_df$gene <- toupper(as.character(rank_df$gene))
  rank_df <- rank_df[!is.na(rank_df$gene) & rank_df$gene != "" & !is.na(rank_df$rank_value), , drop = FALSE]
  rank_df <- rank_df[order(abs(rank_df$rank_value), decreasing = TRUE), , drop = FALSE]
  rank_df <- rank_df[!duplicated(rank_df$gene), , drop = FALSE]

  ranks <- rank_df$rank_value
  names(ranks) <- rank_df$gene
  sort(ranks, decreasing = TRUE)
}

# test whether pathway names match selected biological terms
matches_any_pattern <- function(x, patterns, ignore.case = TRUE) {
  grepl(paste(patterns, collapse = "|"), x, ignore.case = ignore.case, perl = TRUE)
}

# run fgsea for all tested cell types
run_fgsea_pseudobulk <- function(summary_df,
                                 pseudobulk_res_dir,
                                 fgsea_table_dir,
                                 fgsea_iron_table_dir,
                                 fgsea_curve_dir,
                                 pathways,
                                 pathway_info,
                                 iron_related_patterns,
                                 padj_thr,
                                 min_size,
                                 max_size,
                                 fgsea_seed) {
  all_sig_results <- list()
  all_iron_sig_results <- list()

  for (cell_type in summary_df$cell_type) {
    safe_ct <- safe_name(cell_type)
    deg_file <- file.path(pseudobulk_res_dir, paste0("deseq2_HA_vs_other_", safe_ct, ".csv"))
    deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)
    ranks <- build_rank_vector(deg_df)

    set.seed(fgsea_seed)
    fgsea_res <- fgseaMultilevel(pathways = pathways, stats = ranks, minSize = min_size, maxSize = max_size, eps = 0, nproc = 1)
    fgsea_res <- as.data.frame(fgsea_res)
    fgsea_res <- fgsea_res[order(fgsea_res$padj, fgsea_res$pval), , drop = FALSE]
    fgsea_res$cell_type <- cell_type
    fgsea_res$database <- as.character(pathway_info[fgsea_res$pathway, "database", drop = TRUE])
    fgsea_res$pathway_name <- as.character(pathway_info[fgsea_res$pathway, "gs_name", drop = TRUE])
    fgsea_res$is_iron_related <- matches_any_pattern(fgsea_res$pathway_name, iron_related_patterns)
    fgsea_res$leadingEdge <- vapply(fgsea_res$leadingEdge, paste, character(1), collapse = ";")
    fgsea_res <- fgsea_res[, c("cell_type", "database", "pathway_name", "is_iron_related", "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge"), drop = FALSE]

    sig_res <- fgsea_res[!is.na(fgsea_res$padj) & fgsea_res$padj < padj_thr, , drop = FALSE]
    all_sig_results[[cell_type]] <- sig_res

    for (db in unique(sig_res$database)) {
      sig_db <- sig_res[sig_res$database == db, , drop = FALSE]
      write.csv(sig_db, file.path(fgsea_table_dir, paste0("fgsea_ha_vs_other_", safe_ct, "_", db, "_significant.csv")), row.names = FALSE)
    }

    iron_sig_res <- sig_res[sig_res$is_iron_related, , drop = FALSE]
    all_iron_sig_results[[cell_type]] <- iron_sig_res

    if (nrow(iron_sig_res) > 0) {
      write.csv(iron_sig_res, file.path(fgsea_iron_table_dir, paste0("fgsea_ha_vs_other_", safe_ct, "_iron_related_significant.csv")), row.names = FALSE)

      for (i in seq_len(nrow(iron_sig_res))) {
        pathway_id <- iron_sig_res$pathway[i]
        pathway_name <- iron_sig_res$pathway_name[i]
        curve_name <- paste0(safe_name(iron_sig_res$database[i]), "_", sprintf("%03d", i), "_", substr(safe_name(pathway_name), 1, 80))

        p_curve <- plotEnrichment(pathways[[pathway_id]], ranks) +
          labs(
            title = pathway_name,
            subtitle = paste0(cell_type, " | ", iron_sig_res$database[i], " | NES = ", round(iron_sig_res$NES[i], 2), " | padj = ", signif(iron_sig_res$padj[i], 3)),
            x = "Ranked genes",
            y = "Enrichment score"
          ) +
          theme_classic(base_size = 16) +
          theme(
            axis.text = element_text(size = 14, color = "black"),
            axis.title = element_text(size = 16, color = "black"),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
            plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
            panel.grid = element_blank(),
            plot.margin = margin(18, 18, 18, 18)
          )

        ggsave(file.path(fgsea_curve_dir, paste0("gsea_curve_", safe_ct, "_", curve_name, ".png")), p_curve, width = 8, height = 6, dpi = 600)
      }
    }
  }

  all_sig_df <- do.call(rbind, all_sig_results)
  rownames(all_sig_df) <- NULL
  write.csv(all_sig_df, file.path(dirname(fgsea_table_dir), "fgsea_ha_vs_other_all_celltypes_significant.csv"), row.names = FALSE)

  all_iron_sig_df <- do.call(rbind, all_iron_sig_results)
  rownames(all_iron_sig_df) <- NULL
  write.csv(all_iron_sig_df, file.path(dirname(fgsea_table_dir), "fgsea_ha_vs_other_all_celltypes_iron_related_significant.csv"), row.names = FALSE)
}

# clean Enrichr term labels
clean_enrichr_terms <- function(df) {
  df$Term <- gsub("\\s+R-HSA-[0-9]+$", "", df$Term)
  df$Term <- gsub("\\s+KEGG_[0-9]+$", "", df$Term)
  df$Term <- gsub("\\s*\\([^\\)]+\\)", "", df$Term)
  df$Term <- trimws(df$Term)
  df
}

# wrap long pathway labels
wrap_pathway_label <- function(x, n_words = 5) {
  vapply(strsplit(x, " ", fixed = TRUE), function(words) {
    paste(tapply(words, ceiling(seq_along(words) / n_words), paste, collapse = " "), collapse = "\n")
  }, character(1))
}

# prepare Enrichr table for plotting
prepare_enrichr_plot_df <- function(df, input_gene_count, padj_cutoff, top_n_terms, excluded_term_patterns) {
  df <- clean_enrichr_terms(as.data.frame(df))
  df <- df[!matches_any_pattern(df$Term, excluded_term_patterns, ignore.case = FALSE), , drop = FALSE]
  df$num_genes <- sapply(strsplit(df$Genes, ";"), length)
  df$gene_ratio <- df$num_genes / input_gene_count
  df <- df[!is.na(df$Adjusted.P.value) & df$Adjusted.P.value < padj_cutoff, , drop = FALSE]
  df <- head(df[order(df$Adjusted.P.value), , drop = FALSE], top_n_terms)
  df$Term <- wrap_pathway_label(df$Term)
  df <- df[order(df$gene_ratio), , drop = FALSE]
  df$Term <- factor(df$Term, levels = df$Term)
  df
}

# plot Enrichr pathway dotplot
make_enrichr_dotplot <- function(plot_df) {
  ggplot(plot_df, aes(x = gene_ratio, y = Term)) +
    geom_point(aes(size = num_genes, color = Adjusted.P.value), alpha = 0.9) +
    scale_color_continuous(low = "#8B3E2F", high = "coral", name = "adj p-value") +
    scale_size_continuous(name = "Genes", range = c(3, 9)) +
    labs(x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
}

# make tile plot with genes from selected Enrichr pathways
make_selected_pathway_gene_tile <- function(pathway_df, deg_df) {
  term_levels <- levels(pathway_df$Term)
  gene_tiles <- list()

  for (i in seq_len(nrow(pathway_df))) {
    genes <- unique(trimws(strsplit(pathway_df$Genes[i], ";", fixed = TRUE)[[1]]))
    genes <- genes[genes != ""]
    gene_tiles[[i]] <- data.frame(Term = pathway_df$Term[i], gene = genes, gene_index = seq_along(genes), stringsAsFactors = FALSE)
  }

  gene_tiles <- do.call(rbind, gene_tiles)
  gene_tiles$gene_row <- ifelse(gene_tiles$gene_index %% 2 == 0, 2, 1)
  gene_tiles$gene_col <- ceiling(gene_tiles$gene_index / 2)
  gene_tiles <- merge(gene_tiles, unique(deg_df[, c("gene", "log2FoldChange")]), by = "gene", all.x = TRUE, sort = FALSE)
  gene_tiles$Term <- factor(gene_tiles$Term, levels = term_levels)
  gene_tiles$Term_index <- as.numeric(gene_tiles$Term)
  gene_tiles$y_center <- gene_tiles$Term_index + ifelse(gene_tiles$gene_row == 1, -0.22, 0.22)
  gene_tiles$x_center <- gene_tiles$gene_col

  max_cols <- max(gene_tiles$gene_col, na.rm = TRUE)
  max_abs_lfc <- max(abs(gene_tiles$log2FoldChange), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) max_abs_lfc <- 1

  ggplot(gene_tiles, aes(x = x_center, y = y_center, fill = log2FoldChange)) +
    geom_tile(width = 0.96, height = 0.38, color = "black", linewidth = 0.35) +
    geom_text(aes(label = gene), size = 3.4, color = "black") +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, limits = c(-max_abs_lfc, max_abs_lfc), oob = scales::squish, name = "log2FC") +
    scale_x_continuous(limits = c(0.5, max_cols + 0.5), breaks = seq_len(max_cols), labels = rep("", max_cols), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(term_levels), labels = rep("", length(term_levels)), expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank(),
      legend.position = "top"
    )
}

# save selected Enrichr pathway tables and gene panels
save_selected_enrichr_pathways <- function(enrich_res,
                                           deg_df,
                                           cell_type,
                                           safe_ct,
                                           term_patterns,
                                           panel_name,
                                           selected_celltypes,
                                           padj_cutoff,
                                           db_labels,
                                           table_dir,
                                           fig_dir) {
  if (!(cell_type %in% selected_celltypes)) return(invisible(NULL))

  plot_rows <- list()
  selected_rows <- list()

  for (db in names(enrich_res)) {
    df <- clean_enrichr_terms(as.data.frame(enrich_res[[db]]))
    df$num_genes <- sapply(strsplit(df$Genes, ";"), length)
    df$gene_ratio <- df$num_genes / max(1, length(unique(deg_df$gene)))
    df <- df[matches_any_pattern(df$Term, term_patterns) & !is.na(df$Adjusted.P.value) & df$Adjusted.P.value < padj_cutoff, , drop = FALSE]
    if (nrow(df) == 0) next

    df <- df[order(df$gene_ratio, df$Adjusted.P.value), , drop = FALSE]
    df$Term <- paste0(df$Term, " (", db_labels[db], ")")
    df$Term <- wrap_pathway_label(df$Term)
    plot_rows[[db]] <- df

    for (i in seq_len(nrow(df))) {
      genes_here <- unique(trimws(strsplit(df$Genes[i], ";", fixed = TRUE)[[1]]))
      genes_here <- genes_here[genes_here != ""]
      deg_subset <- deg_df[deg_df$gene %in% genes_here, c("gene", "log2FoldChange", "padj"), drop = FALSE]
      deg_subset$Term <- df$Term[i]
      deg_subset$Adjusted.P.value <- df$Adjusted.P.value[i]
      selected_rows[[paste(db, i, sep = "_")]] <- deg_subset
    }
  }

  if (length(plot_rows) == 0 || length(selected_rows) == 0) return(invisible(NULL))

  pathway_plot_df <- do.call(rbind, plot_rows)
  pathway_plot_df <- pathway_plot_df[order(pathway_plot_df$gene_ratio, pathway_plot_df$Adjusted.P.value), , drop = FALSE]
  pathway_plot_df$Term <- factor(pathway_plot_df$Term, levels = unique(pathway_plot_df$Term))

  selected_df <- unique(do.call(rbind, selected_rows))
  write.csv(selected_df, file.path(table_dir, paste0("selected_", panel_name, "_pathway_genes_", safe_ct, ".csv")), row.names = FALSE)

  p_dot <- make_enrichr_dotplot(pathway_plot_df)
  p_genes <- make_selected_pathway_gene_tile(pathway_plot_df, deg_df)
  combined_plot <- (p_dot + p_genes + patchwork::plot_layout(widths = c(1, 1.6), guides = "collect")) &
    theme(legend.position = "top", legend.box = "horizontal")

  n_terms <- length(unique(selected_df$Term))
  panel_width <- ifelse(n_terms == 1, 12, 18)
  panel_height <- ifelse(n_terms == 1, 3.9, max(4.8, 1.5 * n_terms))

  ggsave(file.path(fig_dir, paste0("selected_", panel_name, "_pathway_gene_panel_", safe_ct, ".png")), combined_plot, width = panel_width, height = panel_height, dpi = 300)
}

# run Enrichr over-representation analysis for all suitable cell types
run_enrichr_pseudobulk <- function(summary_df,
                                   pseudobulk_res_dir,
                                   enrichr_table_dir,
                                   selected_table_dir,
                                   selected_fig_dir,
                                   dbs,
                                   db_labels,
                                   min_degs_per_celltype,
                                   padj_cutoff,
                                   selected_celltypes,
                                   iron_related_patterns,
                                   tgf_beta_patterns) {
  summary_df <- summary_df[summary_df$status == "tested" & summary_df$n_sig_padj >= min_degs_per_celltype, , drop = FALSE]

  for (cell_type in summary_df$cell_type) {
    safe_ct <- safe_name(cell_type)
    deg_file <- file.path(pseudobulk_res_dir, paste0("deseq2_HA_vs_other_", safe_ct, "_significant.csv"))
    deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)
    deg_df <- deg_df[!is.na(deg_df$padj) & deg_df$padj < padj_cutoff, , drop = FALSE]
    deg_genes <- unique(deg_df$gene)
    deg_genes <- deg_genes[!is.na(deg_genes) & deg_genes != ""]

    enrich_res <- enrichr(deg_genes, dbs)

    for (db in names(enrich_res)) {
      write.csv(enrich_res[[db]], file.path(enrichr_table_dir, paste0("enrichr_", safe_ct, "_", db, ".csv")), row.names = FALSE)
    }

    save_selected_enrichr_pathways(enrich_res, deg_df, cell_type, safe_ct, iron_related_patterns, "iron_related", selected_celltypes, padj_cutoff, db_labels, selected_table_dir, selected_fig_dir)
    save_selected_enrichr_pathways(enrich_res, deg_df, cell_type, safe_ct, tgf_beta_patterns, "tgf_beta_related", selected_celltypes, padj_cutoff, db_labels, selected_table_dir, selected_fig_dir)
  }
}
