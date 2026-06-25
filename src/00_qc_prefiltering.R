#######################
### QC prefiltering ###
#######################

# scRNA-seq samples listed in metadata/subject_info.xlsx
# script useful to choose QC thresholds for filtering low-quality cells and doublets

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

# set wd
setwd("~/Thema_R")
source("src/global_config.R")

# dir
qc_fig_dir <- file.path(figures_dir, "qc_prefiltering")
dir.create(qc_fig_dir, recursive = TRUE, showWarnings = FALSE)
qc_table_dir <- file.path(qc_dir, "pre_filtering")
dir.create(qc_table_dir, recursive = TRUE, showWarnings = FALSE)

# load sample metadata
subject_info_path <- file.path(metadata_dir, "subject_info.xlsx")
if (!file.exists(subject_info_path)) stop("Missing: ", subject_info_path)
# clean column names
subject_info <- as.data.frame(read_excel(subject_info_path))
colnames(subject_info) <- trimws(colnames(subject_info))
subject_info[] <- lapply(subject_info, function(x) if (is.character(x)) trimws(x) else x)
# match metadata with raw count folders
if (!"sample_name" %in% colnames(subject_info)) {
  stop("subject_info.xlsx needs a 'sample_name' column")
}

# remove empty or duplicated sample
subject_info$sample_name <- trimws(as.character(subject_info$sample_name))
subject_info <- subject_info[subject_info$sample_name != "" & !is.na(subject_info$sample_name), , drop = FALSE]
subject_info <- subject_info[!duplicated(subject_info$sample_name), , drop = FALSE]
rownames(subject_info) <- subject_info$sample_name

# use only the samples listed in subject_info.xlsx
samples <- subject_info$sample_name

cat("Samples in subject_info.xlsx:", length(samples), "\n")
cat("Samples to QC:", paste(samples, collapse = ", "), "\n")

# calculate the same QC quantiles for every metric
qc_quantiles_vec <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, 7))
  as.numeric(stats::quantile(
    x,
    probs = c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99),
    na.rm = TRUE
  ))
}

# calculate the median
qc_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

# create one small quantile table for nCount, nFeature and mt RNA
qc_quantile_table <- function(df, group_label) {
  ncount_q <- qc_quantiles_vec(df$nCount_RNA)
  nfeature_q <- qc_quantiles_vec(df$nFeature_RNA)
  mito_q <- qc_quantiles_vec(df$percent.mt)

  data.frame(
    group = group_label,
    metric = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    q01 = c(ncount_q[1], nfeature_q[1], mito_q[1]),
    q05 = c(ncount_q[2], nfeature_q[2], mito_q[2]),
    q10 = c(ncount_q[3], nfeature_q[3], mito_q[3]),
    q50 = c(ncount_q[4], nfeature_q[4], mito_q[4]),
    q90 = c(ncount_q[5], nfeature_q[5], mito_q[5]),
    q95 = c(ncount_q[6], nfeature_q[6], mito_q[6]),
    q99 = c(ncount_q[7], nfeature_q[7], mito_q[7])
  )
}

# find mt genes using human or mouse-style gene prefixes
mitochondrial_features <- function(feature_names) {
  grep("^(MT-|mt-|Mt-)", feature_names, value = TRUE)
}
# add percent.mt to the Seurat metadata
add_percent_mt <- function(obj, sample_id) {
  mt_features <- mitochondrial_features(rownames(obj))
  if (length(mt_features) == 0) {
    warning("No mitochondrial features matching ^MT-, ^mt- or ^Mt- found for ", sample_id)
    obj$percent.mt <- NA_real_
  } else {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_features)
  }
  obj
}

# read one 10x matrix from the raw feature-barcode folder
read_sample_counts <- function(sample_id) {
  matrix_dir <- file.path(raw_data_dir, paste0(sample_id, "_raw_feature_bc_matrix"))

  if (!dir.exists(matrix_dir)) stop("Missing 10x input folder: ", matrix_dir)

  counts <- Read10X(data.dir = matrix_dir)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      warning("Gene Expression matrix not named explicitly for ", sample_id, "; using first matrix: ", names(counts)[1])
      counts <- counts[[1]]
    }
  }

  cat("Loaded sample", sample_id, "from", matrix_dir, "\n")
  counts
}

# save one histogram panel per sample for the main QC metrics
save_hist_plot <- function(df, sample_id) {
  # log-scaled histograms cannot use zero values
  df_nonzero <- df[df$nCount_RNA > 0 & df$nFeature_RNA > 0, , drop = FALSE]
  if (nrow(df_nonzero) == 0) df_nonzero <- df

  # percent.mt can be missing if mitochondrial genes were not found
  df_mito <- df_nonzero[is.finite(df_nonzero$percent.mt), , drop = FALSE]

  p1 <- ggplot(df_nonzero, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 80, fill = "grey80", color = "white") +
    scale_x_log10() +
    scale_y_sqrt() +
    labs(title = "nFeature_RNA", x = NULL, y = "Cells") +
    theme_classic()

  p2 <- ggplot(df_nonzero, aes(x = nCount_RNA)) +
    geom_histogram(bins = 80, fill = "grey80", color = "white") +
    scale_x_log10() +
    scale_y_sqrt() +
    labs(title = "nCount_RNA", x = NULL, y = "Cells") +
    theme_classic()

  if (nrow(df_mito) == 0) {
    p3 <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No MT-* features found") +
      labs(title = "percent.mt", x = NULL, y = NULL) +
      theme_classic() +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  } else {
    p3 <- ggplot(df_mito, aes(x = percent.mt)) +
      geom_histogram(bins = 80, fill = "grey80", color = "white") +
      scale_y_sqrt() +
      labs(title = "percent.mt", x = NULL, y = "Cells") +
      theme_classic()
  }

  hist_plot <- (p1 | p2 | p3) +
    plot_annotation(
      title = paste0(sample_id, " | RNA QC histograms"),
      subtitle = "Non-zero barcodes only"
    )

  ggsave(
    file.path(qc_fig_dir, paste0(sample_id, "_QC_hist.png")),
    plot = hist_plot,
    width = 13,
    height = 4.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

summary_rows <- list()

# process each sample independently
for (s in samples) {
  cat("\n==============================\n")
  cat("QC prefiltering for sample:", s, "\n")
  cat("==============================\n")

  # load raw counts and create a temporary Seurat object
  counts <- read_sample_counts(sample_id = s)

  obj <- CreateSeuratObject(counts = counts, assay = "RNA", project = s)
  obj$sample_id <- s

  # copy sample-level metadata into the Seurat object
  si <- subject_info[s, , drop = FALSE]
  for (cn in colnames(si)) {
    if (cn == "sample_name") next
    obj[[cn]] <- si[[cn]][1]
  }

  # calculate mt percentage and extract cell-level QC metadata
  obj <- add_percent_mt(obj, sample_id = s)
  md <- obj@meta.data

  cat("\n--- Basic counts (raw RNA) ---\n")
  cat("Cells:", ncol(obj), "\n")
  cat("Genes (RNA features):", nrow(obj[["RNA"]]), "\n")
  cat("Non-zero barcodes:", sum(md$nCount_RNA > 0 & md$nFeature_RNA > 0), "\n")

  # summarize QC both before and after removing zero-count barcodes
  md_nonzero <- md[md$nCount_RNA > 0 & md$nFeature_RNA > 0, , drop = FALSE]
  qtab <- rbind(
    qc_quantile_table(md, "all_raw_barcodes"),
    qc_quantile_table(md_nonzero, "nonzero_barcodes")
  )

  cat("\n--- RNA QC quantiles (1%, 5%, 10%, 50%, 90%, 95%, 99%) ---\n")
  print(qtab)
  write.csv(qtab, file.path(qc_table_dir, paste0(s, "_QC_quantiles.csv")), row.names = FALSE)
  write.csv(qtab, file.path(qc_fig_dir, paste0(s, "_QC_quantiles.csv")), row.names = FALSE)

  # scatter plot helps detect low-complexity cells and high-count outliers
  p_scatter <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle(paste0(s, " | RNA: nCount_RNA vs nFeature_RNA")) +
    theme_classic()

  ggsave(
    file.path(qc_fig_dir, paste0(s, "_QC_scatter.png")),
    plot = p_scatter,
    width = 6,
    height = 5,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  save_hist_plot(md, sample_id = s)

  # keep one row per sample for the final summary table
  summary_rows[[s]] <- data.frame(
    sample_id = s,
    n_cells_raw = ncol(obj),
    n_barcodes_nonzero = sum(md$nCount_RNA > 0 & md$nFeature_RNA > 0),
    n_genes = nrow(obj[["RNA"]]),
    median_nCount_RNA = qc_median(md$nCount_RNA),
    median_nFeature_RNA = qc_median(md$nFeature_RNA),
    median_percent_mt = qc_median(md$percent.mt)
  )

  # remove the temporary obj before loading the next sample
  rm(obj, counts)
  gc()
}

# save the global QC summary across samples
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(qc_table_dir, "QC_prefiltering_summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(qc_fig_dir, "QC_prefiltering_summary.csv"), row.names = FALSE)

cat("\n========================================\n")
cat("RNA QC prefiltering completed.\n")
cat("QC plots:          ", qc_fig_dir, "\n")
cat("QC tables:         ", qc_table_dir, "\n")
cat("========================================\n")
