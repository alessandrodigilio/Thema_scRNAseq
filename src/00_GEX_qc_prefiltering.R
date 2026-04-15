# ---------------------------------------
# 00_GEX_qc_prefiltering.R
# RNA QC before filtering
# ---------------------------------------

# Inspect raw scRNA-seq GEX samples listed in subject_info.xlsx before
# choosing filtering thresholds.

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

# Work from the project root.
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Output directories
QC_FIG_DIR <- file.path(FIGURES_DIR, "qc_prefiltering")
dir.create(QC_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

QC_TABLE_DIR <- file.path(QC_DIR, "pre_filtering")
dir.create(QC_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Metadata
subject_info_path <- file.path(METADATA_DIR, "subject_info.xlsx")
if (!file.exists(subject_info_path)) stop("Missing: ", subject_info_path)

subject_info <- as.data.frame(readxl::read_excel(subject_info_path))
colnames(subject_info) <- trimws(colnames(subject_info))
subject_info[] <- lapply(subject_info, function(x) if (is.character(x)) trimws(x) else x)

if (!"sample_name" %in% colnames(subject_info)) {
  stop("subject_info.xlsx needs a 'sample_name' column")
}

subject_info$sample_name <- trimws(as.character(subject_info$sample_name))
subject_info <- subject_info[subject_info$sample_name != "" & !is.na(subject_info$sample_name), , drop = FALSE]
subject_info <- subject_info[!duplicated(subject_info$sample_name), , drop = FALSE]
rownames(subject_info) <- subject_info$sample_name

SAMPLES <- setdiff(subject_info$sample_name, EXCLUDED_SAMPLES)

cat("Samples requested in subject_info.xlsx:", length(subject_info$sample_name), "\n")
cat("Samples to QC:", length(SAMPLES), "\n")
cat("Included samples:", paste(SAMPLES, collapse = ", "), "\n")

# Report raw-count inputs that are present but not used by this project.
raw_entries <- list.files(GEX_RAW_DATA_DIR, full.names = FALSE)
raw_sample_ids <- unique(gsub("(_raw_feature_bc_matrix)?(\\.h5)?$", "", raw_entries))
unused_raw <- setdiff(raw_sample_ids, SAMPLES)
if (length(unused_raw) > 0) {
  cat("Raw-count samples not in subject_info.xlsx and skipped:", paste(sort(unused_raw), collapse = ", "), "\n")
}

# Helpers
qc_quantiles_vec <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, 7))
  as.numeric(stats::quantile(
    x,
    probs = c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99),
    na.rm = TRUE
  ))
}

qc_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

qc_quantile_table <- function(df, group_label) {
  data.frame(
    group = group_label,
    metric = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    q01 = c(qc_quantiles_vec(df$nCount_RNA)[1], qc_quantiles_vec(df$nFeature_RNA)[1], qc_quantiles_vec(df$percent.mt)[1]),
    q05 = c(qc_quantiles_vec(df$nCount_RNA)[2], qc_quantiles_vec(df$nFeature_RNA)[2], qc_quantiles_vec(df$percent.mt)[2]),
    q10 = c(qc_quantiles_vec(df$nCount_RNA)[3], qc_quantiles_vec(df$nFeature_RNA)[3], qc_quantiles_vec(df$percent.mt)[3]),
    q50 = c(qc_quantiles_vec(df$nCount_RNA)[4], qc_quantiles_vec(df$nFeature_RNA)[4], qc_quantiles_vec(df$percent.mt)[4]),
    q90 = c(qc_quantiles_vec(df$nCount_RNA)[5], qc_quantiles_vec(df$nFeature_RNA)[5], qc_quantiles_vec(df$percent.mt)[5]),
    q95 = c(qc_quantiles_vec(df$nCount_RNA)[6], qc_quantiles_vec(df$nFeature_RNA)[6], qc_quantiles_vec(df$percent.mt)[6]),
    q99 = c(qc_quantiles_vec(df$nCount_RNA)[7], qc_quantiles_vec(df$nFeature_RNA)[7], qc_quantiles_vec(df$percent.mt)[7])
  )
}

mitochondrial_features <- function(feature_names) {
  grep("^(MT-|mt-|Mt-)", feature_names, value = TRUE)
}

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

read_sample_counts <- function(sample_id) {
  sample_path <- paste0(sample_id, "_raw_feature_bc_matrix")
  matrix_dir <- file.path(GEX_RAW_DATA_DIR, sample_path)
  h5_path <- paste0(matrix_dir, ".h5")

  counts <- NULL
  source_used <- NULL

  if (dir.exists(matrix_dir)) {
    counts <- Read10X(data.dir = matrix_dir)
    source_used <- matrix_dir
  } else if (file.exists(h5_path)) {
    counts <- Read10X_h5(filename = h5_path)
    source_used <- h5_path
  } else {
    fallback_dir <- file.path(GEX_RAW_DATA_DIR, paste0(sample_id, "_raw_feature_bc_matrix"))
    fallback_h5 <- paste0(fallback_dir, ".h5")

    if (dir.exists(fallback_dir)) {
      counts <- Read10X(data.dir = fallback_dir)
      source_used <- fallback_dir
    } else if (file.exists(fallback_h5)) {
      counts <- Read10X_h5(filename = fallback_h5)
      source_used <- fallback_h5
    } else {
      stop(
        "Missing 10x input for sample ", sample_id,
        ". Checked: ", matrix_dir, ", ", h5_path, ", ", fallback_dir, " and ", fallback_h5
      )
    }
  }

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      warning("Gene Expression matrix not named explicitly for ", sample_id, "; using first matrix: ", names(counts)[1])
      counts <- counts[[1]]
    }
  }

  cat("Loaded sample", sample_id, "from", source_used, "\n")
  counts
}

save_hist_plot <- function(df, sample_id) {
  df_nonzero <- df[df$nCount_RNA > 0 & df$nFeature_RNA > 0, , drop = FALSE]
  if (nrow(df_nonzero) == 0) df_nonzero <- df
  df_mito <- df_nonzero[is.finite(df_nonzero$percent.mt), , drop = FALSE]

  p1 <- ggplot2::ggplot(df_nonzero, ggplot2::aes(x = nFeature_RNA)) +
    ggplot2::geom_histogram(bins = 80, fill = "grey80", color = "white") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_sqrt() +
    ggplot2::labs(title = "nFeature_RNA", x = NULL, y = "Cells") +
    ggplot2::theme_classic()

  p2 <- ggplot2::ggplot(df_nonzero, ggplot2::aes(x = nCount_RNA)) +
    ggplot2::geom_histogram(bins = 80, fill = "grey80", color = "white") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_sqrt() +
    ggplot2::labs(title = "nCount_RNA", x = NULL, y = "Cells") +
    ggplot2::theme_classic()

  if (nrow(df_mito) == 0) {
    p3 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No MT-* features found") +
      ggplot2::labs(title = "percent.mt", x = NULL, y = NULL) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  } else {
    p3 <- ggplot2::ggplot(df_mito, ggplot2::aes(x = percent.mt)) +
      ggplot2::geom_histogram(bins = 80, fill = "grey80", color = "white") +
      ggplot2::scale_y_sqrt() +
      ggplot2::labs(title = "percent.mt", x = NULL, y = "Cells") +
      ggplot2::theme_classic()
  }

  hist_plot <- (p1 | p2 | p3) +
    patchwork::plot_annotation(
      title = paste0(sample_id, " | RNA QC histograms"),
      subtitle = "Non-zero barcodes only"
    )

  ggplot2::ggsave(
    file.path(QC_FIG_DIR, paste0(sample_id, "_GEX_QC_hist.png")),
    plot = hist_plot,
    width = 13,
    height = 4.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

summary_rows <- list()

for (s in SAMPLES) {
  cat("\n==============================\n")
  cat("QC prefiltering for sample:", s, "\n")
  cat("==============================\n")

  counts <- read_sample_counts(sample_id = s)

  obj <- CreateSeuratObject(counts = counts, assay = "RNA", project = s)
  obj$sample_id <- s

  if (s %in% rownames(subject_info)) {
    si <- subject_info[s, , drop = FALSE]
    for (cn in colnames(si)) {
      if (cn == "sample_name") next
      obj[[cn]] <- si[[cn]][1]
    }
  }

  obj <- add_percent_mt(obj, sample_id = s)
  md <- obj@meta.data

  cat("\n--- Basic counts (raw GEX) ---\n")
  cat("Cells:", ncol(obj), "\n")
  cat("Genes (RNA features):", nrow(obj[["RNA"]]), "\n")
  cat("Non-zero barcodes:", sum(md$nCount_RNA > 0 & md$nFeature_RNA > 0), "\n")

  md_nonzero <- md[md$nCount_RNA > 0 & md$nFeature_RNA > 0, , drop = FALSE]
  qtab <- rbind(
    qc_quantile_table(md, "all_raw_barcodes"),
    qc_quantile_table(md_nonzero, "nonzero_barcodes")
  )

  cat("\n--- GEX QC quantiles (1%, 5%, 10%, 50%, 90%, 95%, 99%) ---\n")
  print(qtab)
  write.csv(qtab, file.path(QC_TABLE_DIR, paste0(s, "_GEX_QC_quantiles.csv")), row.names = FALSE)
  write.csv(qtab, file.path(QC_FIG_DIR, paste0(s, "_GEX_QC_quantiles.csv")), row.names = FALSE)

  p_scatter <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggplot2::ggtitle(paste0(s, " | GEX: nCount_RNA vs nFeature_RNA")) +
    ggplot2::theme_classic()

  ggplot2::ggsave(
    file.path(QC_FIG_DIR, paste0(s, "_GEX_QC_scatter.png")),
    plot = p_scatter,
    width = 6,
    height = 5,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  md_plot <- md
  md_plot$sample_id <- s
  save_hist_plot(md_plot, sample_id = s)

  summary_rows[[s]] <- data.frame(
    sample_id = s,
    n_cells_raw = ncol(obj),
    n_barcodes_nonzero = sum(md$nCount_RNA > 0 & md$nFeature_RNA > 0),
    n_genes = nrow(obj[["RNA"]]),
    median_nCount_RNA = qc_median(md$nCount_RNA),
    median_nFeature_RNA = qc_median(md$nFeature_RNA),
    median_percent_mt = qc_median(md$percent.mt)
  )

  rm(obj, counts)
  gc()
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(QC_TABLE_DIR, "GEX_QC_prefiltering_summary.csv"), row.names = FALSE)
write.csv(summary_df, file.path(QC_FIG_DIR, "GEX_QC_prefiltering_summary.csv"), row.names = FALSE)

cat("\n========================================\n")
cat("GEX QC prefiltering completed.\n")
cat("QC plots:          ", QC_FIG_DIR, "\n")
cat("QC tables:         ", QC_TABLE_DIR, "\n")
cat("========================================\n")
