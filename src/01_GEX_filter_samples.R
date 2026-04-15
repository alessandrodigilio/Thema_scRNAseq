# ---------------------------------------
# 01_GEX_filter_samples.R
# Filter scRNA-seq GEX samples
# ---------------------------------------

# Load scRNA-seq samples listed in subject_info.xlsx, clean RNA with SoupX,
# filter cells by per-sample QC thresholds, remove RNA doublets with
# scDblFinder, and save filtered Seurat objects plus a merged object.

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(SoupX)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scater)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

# Work from the project root.
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
OBJECTS_DIR <- FILTERED_DATA_DIR
dir.create(OBJECTS_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_DIR <- FILTERING_LOG_DIR
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_FILTERING_DIR <- file.path(FIGURES_DIR, "filtering")
dir.create(FIG_FILTERING_DIR, recursive = TRUE, showWarnings = FALSE)

# Load per-sample filtering thresholds
thr_path <- file.path(METADATA_DIR, "filters.xlsx")
if (!file.exists(thr_path)) stop("Missing thresholds Excel: ", thr_path)

thr <- as.data.frame(readxl::read_excel(thr_path))
colnames(thr) <- trimws(colnames(thr))
thr[] <- lapply(thr, function(x) if (is.character(x)) trimws(x) else x)

required_thr_cols <- c(
  "sample_name", "sample_path", "min_genes", "max_genes",
  "min_counts", "max_counts", "max_pct_mt"
)
if (!all(required_thr_cols %in% colnames(thr))) {
  stop(
    "filters.xlsx is missing columns: ",
    paste(setdiff(required_thr_cols, colnames(thr)), collapse = ", ")
  )
}

thr$sample_name <- trimws(as.character(thr$sample_name))

# Load sample annotations and define the only samples to process
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

SAMPLES <- subject_info$sample_name

missing_thr <- setdiff(SAMPLES, thr$sample_name)
if (length(missing_thr) > 0) {
  stop("Missing filtering thresholds for samples: ", paste(missing_thr, collapse = ", "))
}

SAMPLES <- setdiff(SAMPLES, EXCLUDED_SAMPLES)
SKIP_COMPLETED_SAMPLES <- c("10TA")

cat("Samples requested in subject_info.xlsx:", length(subject_info$sample_name), "\n")
cat("Samples to process:", length(SAMPLES), "\n")
cat("Included samples:", paste(SAMPLES, collapse = ", "), "\n")
cat("Completed samples to reuse:", paste(SKIP_COMPLETED_SAMPLES, collapse = ", "), "\n")

# Apply lower thresholds while skipping missing cutoffs
apply_min <- function(x, thr_val) {
  if (is.na(thr_val)) return(rep(TRUE, length(x)))
  x >= thr_val
}

# Apply upper thresholds while skipping missing cutoffs
apply_max <- function(x, thr_val) {
  if (is.na(thr_val)) return(rep(TRUE, length(x)))
  x <= thr_val
}

# Extract the first non-missing numeric value from optional SoupX slots
first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(x[[1]])
}

# Read the estimated contamination fraction from the SoupX object
extract_soupx_rho <- function(soupx_obj) {
  first_non_na(c(
    tryCatch(soupx_obj$metaData$rho[1], error = function(e) NA_real_),
    tryCatch(soupx_obj$fit$rhoEst, error = function(e) NA_real_),
    tryCatch(soupx_obj$fit$rho, error = function(e) NA_real_)
  ))
}

# Read a 10x sample from either a matrix directory or a .h5 file
read_sample_counts <- function(sample_id, sample_path) {
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
        ". Checked: ", matrix_dir, " and ", h5_path
      )
    }
  }

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  cat("Loaded sample", sample_id, "from", source_used, "\n")
  counts
}

# Read only feature names from the 10x matrix directory when available.
read_sample_features <- function(sample_id, sample_path) {
  matrix_dir <- file.path(GEX_RAW_DATA_DIR, sample_path)
  feature_file <- file.path(matrix_dir, "features.tsv.gz")

  if (!file.exists(feature_file)) {
    fallback_dir <- file.path(GEX_RAW_DATA_DIR, paste0(sample_id, "_raw_feature_bc_matrix"))
    feature_file <- file.path(fallback_dir, "features.tsv.gz")
  }

  if (!file.exists(feature_file)) return(NULL)

  features <- utils::read.table(
    gzfile(feature_file),
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  make.unique(as.character(features[[2]]))
}

# Run SoupX on RNA counts and return corrected counts plus summary metrics
run_soupx_if_available <- function(sample_id, raw_counts, filtered_counts) {
  umi_before <- sum(filtered_counts)
  genes_before <- Matrix::colSums(filtered_counts > 0)
  zero_removed <- stats::setNames(rep(0, ncol(filtered_counts)), colnames(filtered_counts))

  if (!isTRUE(RUN_SOUPX)) {
    return(list(
      counts = filtered_counts,
      removed_umis = zero_removed,
      removed_fraction = zero_removed,
      stats = list(
        soupx_status = "disabled",
        soupx_rho = NA_real_,
        soupx_umi_before = umi_before,
        soupx_umi_after = umi_before,
        soupx_umis_removed = 0,
        soupx_fraction_removed = 0,
        soupx_median_genes_before = median(genes_before),
        soupx_median_genes_after = median(genes_before)
      )
    ))
  }

  if (ncol(filtered_counts) == 0) {
    return(list(
      counts = filtered_counts,
      removed_umis = zero_removed,
      removed_fraction = zero_removed,
      stats = list(
        soupx_status = "no_cells_after_prefilter",
        soupx_rho = NA_real_,
        soupx_umi_before = umi_before,
        soupx_umi_after = umi_before,
        soupx_umis_removed = 0,
        soupx_fraction_removed = 0,
        soupx_median_genes_before = NA_real_,
        soupx_median_genes_after = NA_real_
      )
    ))
  }

  common_genes <- intersect(rownames(raw_counts), rownames(filtered_counts))
  if (length(common_genes) == 0) {
    warning("SoupX skipped for ", sample_id, ": no shared genes between raw and filtered RNA matrices")
    return(list(
      counts = filtered_counts,
      removed_umis = zero_removed,
      removed_fraction = zero_removed,
      stats = list(
        soupx_status = "no_shared_genes",
        soupx_rho = NA_real_,
        soupx_umi_before = umi_before,
        soupx_umi_after = umi_before,
        soupx_umis_removed = 0,
        soupx_fraction_removed = 0,
        soupx_median_genes_before = median(genes_before),
        soupx_median_genes_after = median(genes_before)
      )
    ))
  }

  cat("Running SoupX for sample:", sample_id, "\n")

  filtered_counts <- filtered_counts[common_genes, , drop = FALSE]
  raw_counts <- raw_counts[common_genes, , drop = FALSE]

  soupx_obj <- SoupX::SoupChannel(tod = raw_counts, toc = filtered_counts)

  tmp <- CreateSeuratObject(counts = filtered_counts, assay = "RNA", project = paste0(sample_id, "_SoupX"))
  tmp <- NormalizeData(tmp, verbose = FALSE)
  tmp <- FindVariableFeatures(tmp, verbose = FALSE)
  tmp <- ScaleData(tmp, verbose = FALSE)
  tmp <- RunPCA(tmp, npcs = max(10, SCDBLFINDER_DIMS), verbose = FALSE)

  dims_use <- seq_len(min(SCDBLFINDER_DIMS, ncol(Embeddings(tmp, "pca"))))
  tmp <- FindNeighbors(tmp, dims = dims_use, verbose = FALSE)
  tmp <- FindClusters(tmp, resolution = 0.2, verbose = FALSE)

  cluster_map <- setNames(as.character(tmp$seurat_clusters), colnames(tmp))
  soupx_obj <- SoupX::setClusters(soupx_obj, cluster_map)
  soupx_obj <- SoupX::autoEstCont(soupx_obj, doPlot = FALSE)

  corrected <- SoupX::adjustCounts(soupx_obj, roundToInt = TRUE)
  corrected <- corrected[, colnames(filtered_counts), drop = FALSE]

  umi_after <- sum(corrected)
  genes_after <- Matrix::colSums(corrected > 0)
  rho_est <- extract_soupx_rho(soupx_obj)
  removed_umis <- Matrix::colSums(filtered_counts) - Matrix::colSums(corrected)
  removed_fraction <- ifelse(
    Matrix::colSums(filtered_counts) > 0,
    removed_umis / Matrix::colSums(filtered_counts),
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

# Run scDblFinder on RNA counts after QC filtering
run_scdblfinder_if_available <- function(obj, sample_id) {
  cells_before <- ncol(obj)

  if (!isTRUE(RUN_SCDBLFINDER)) {
    obj$scDblFinder.class <- NA_character_
    obj$scDblFinder.score <- NA_real_
    return(list(
      obj = obj,
      plot_score = NULL,
      stats = list(
        scdblfinder_status = "disabled",
        scdblfinder_error = NA_character_,
        scdblfinder_singlets = cells_before,
        scdblfinder_doublets_removed = 0
      )
    ))
  }

  if (cells_before < 100) {
    warning("Too few cells for scDblFinder in sample ", sample_id, "; skipping doublet detection")
    obj$scDblFinder.class <- NA_character_
    obj$scDblFinder.score <- NA_real_
    return(list(
      obj = obj,
      plot_score = NULL,
      stats = list(
        scdblfinder_status = "too_few_cells",
        scdblfinder_error = NA_character_,
        scdblfinder_singlets = cells_before,
        scdblfinder_doublets_removed = 0
      )
    ))
  }

  cat("Running scDblFinder for sample:", sample_id, "\n")

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Seurat::GetAssayData(obj, assay = "RNA", layer = "counts"))
  )

  dbl_try <- tryCatch({
    sce <- scDblFinder::scDblFinder(
      sce,
      dims = SCDBLFINDER_DIMS,
      verbose = FALSE
    )

    dbl_class <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
    dbl_score <- as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score)

    obj$scDblFinder.class <- dbl_class
    obj$scDblFinder.score <- dbl_score

    set.seed(10010101)
    sce <- scater::runPCA(sce, ncomponents = SCDBLFINDER_DIMS, exprs_values = "counts")
    sce <- scater::runTSNE(sce, dimred = "PCA")
    plot_score <- scater::plotTSNE(sce, colour_by = "scDblFinder.score") +
      ggplot2::ggtitle(paste0(sample_id, " scDblFinder score"))

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
        scdblfinder_error = NA_character_,
        scdblfinder_singlets = n_singlets,
        scdblfinder_doublets_removed = n_doublets_removed
      )
    )
  }, error = function(e) {
    warning("scDblFinder failed for sample ", sample_id, ": ", conditionMessage(e))
    obj$scDblFinder.class <- NA_character_
    obj$scDblFinder.score <- NA_real_
    list(
      obj = obj,
      plot_score = NULL,
      stats = list(
        scdblfinder_status = "error",
        scdblfinder_error = conditionMessage(e),
        scdblfinder_singlets = cells_before,
        scdblfinder_doublets_removed = 0
      )
    )
  })

  dbl_try
}

# Build and save a two-panel SoupX/scDblFinder summary for one sample
save_filtering_panel <- function(obj_qc, plot_score, sample_id) {
  soupx_qc_fraction <- obj_qc$soupx_removed_fraction
  if (is.null(soupx_qc_fraction)) return(NA_character_)

  plot_soupx <- ggplot2::ggplot(
    data.frame(soupx_removed_fraction = soupx_qc_fraction),
    ggplot2::aes(x = soupx_removed_fraction)
  ) +
    ggplot2::geom_histogram(bins = 40, fill = "grey70", color = "white") +
    ggplot2::labs(
      title = paste0(sample_id, " SoupX removed RNA fraction"),
      x = "Removed fraction",
      y = "Cells"
    ) +
    ggplot2::theme_classic()

  if (is.null(plot_score)) {
    panel_plot <- plot_soupx
  } else {
    panel_plot <- patchwork::wrap_plots(plot_soupx, plot_score, ncol = 2)
  }

  panel_path <- file.path(FIG_FILTERING_DIR, paste0(sample_id, "_soupx_scdblfinder_panel.png"))
  ggplot2::ggsave(
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

# Store filtered Seurat objects before the final merge
objs_filt <- list()

# Accumulate one summary row per sample for the final CSV log
log_rows <- list()

cat("\nComputing common RNA features across samples...\n")
feature_sets <- list()
for (s in SAMPLES) {
  thr_row <- thr[thr$sample_name == s, , drop = FALSE]
  features <- read_sample_features(
    sample_id = s,
    sample_path = thr_row$sample_path[1]
  )

  if (!is.null(features)) {
    feature_sets[[s]] <- features
    cat(s, "features:", length(features), "\n")
  }
}

COMMON_FEATURES <- NULL
if (length(feature_sets) == length(SAMPLES)) {
  COMMON_FEATURES <- Reduce(intersect, feature_sets)
  cat("Common features retained:", length(COMMON_FEATURES), "\n")
} else {
  warning("Could not read features.tsv.gz for all samples; feature intersection will be skipped")
}

for (s in SAMPLES) {
  cat("\n==============================\n")
  cat("Filtering sample:", s, "\n")
  cat("==============================\n")

  out_rds <- file.path(OBJECTS_DIR, paste0(s, "_filtered.rds"))

  if (s %in% SKIP_COMPLETED_SAMPLES) {
    if (!file.exists(out_rds)) {
      stop("Cannot reuse completed sample ", s, ": missing ", out_rds)
    }

    cat("Reusing existing filtered object:", out_rds, "\n")
    obj_f <- readRDS(out_rds)
    objs_filt[[s]] <- obj_f

    panel_path <- file.path(FIG_FILTERING_DIR, paste0(s, "_soupx_scdblfinder_panel.png"))
    if (!file.exists(panel_path)) panel_path <- NA_character_

    log_rows[[s]] <- data.frame(
      sample_id = s,
      status = "reused_existing",
      exclusion_reason = NA_character_,
      n_before = NA_integer_,
      n_after_qc = NA_integer_,
      n_after_scdblfinder = ncol(obj_f),
      frac_kept = NA_real_,
      soupx_status = "reused_existing",
      soupx_rho = NA_real_,
      soupx_umi_before = NA_real_,
      soupx_umi_after = NA_real_,
      soupx_umis_removed = NA_real_,
      soupx_fraction_removed = NA_real_,
      soupx_median_genes_before = NA_real_,
      soupx_median_genes_after = NA_real_,
      scdblfinder_status = "reused_existing",
      scdblfinder_error = NA_character_,
      scdblfinder_singlets = ncol(obj_f),
      scdblfinder_doublets_removed = NA_integer_,
      panel_plot_path = panel_path
    )

    next
  }

  thr_row <- thr[thr$sample_name == s, , drop = FALSE]
  if (nrow(thr_row) != 1) stop("Thresholds not found or duplicated for sample: ", s)

  counts <- read_sample_counts(
    sample_id = s,
    sample_path = thr_row$sample_path[1]
  )

  if (!is.null(COMMON_FEATURES)) {
    counts <- counts[COMMON_FEATURES, , drop = FALSE]
  }

  # Build a temporary object to compute RNA QC and define filtered cells
  obj_raw <- CreateSeuratObject(counts = counts, assay = "RNA", project = s)
  obj_raw$sample_id <- s
  obj_raw[["percent.mt"]] <- PercentageFeatureSet(obj_raw, pattern = "^MT-")

  md_raw <- obj_raw@meta.data
  cells_keep <- rownames(md_raw)[
    apply_min(md_raw$nCount_RNA, as.numeric(thr_row$min_counts)) &
      apply_max(md_raw$nCount_RNA, as.numeric(thr_row$max_counts)) &
      apply_min(md_raw$nFeature_RNA, as.numeric(thr_row$min_genes)) &
      apply_max(md_raw$nFeature_RNA, as.numeric(thr_row$max_genes)) &
      apply_max(md_raw$percent.mt, as.numeric(thr_row$max_pct_mt))
  ]

  filtered_counts <- counts[, cells_keep, drop = FALSE]

  # Clean RNA ambient signal before building the final object
  soupx_res <- run_soupx_if_available(
    sample_id = s,
    raw_counts = counts,
    filtered_counts = filtered_counts
  )
  counts_rna <- soupx_res$counts
  soupx_removed_umis <- soupx_res$removed_umis
  soupx_removed_fraction <- soupx_res$removed_fraction
  soupx_stats <- soupx_res$stats

  # Create the filtered GEX object and attach sample metadata
  obj <- CreateSeuratObject(
    counts = counts_rna,
    assay = "RNA",
    project = s
  )
  obj$sample_id <- s
  obj$soupx_removed_umis <- soupx_removed_umis[colnames(obj)]
  obj$soupx_removed_fraction <- soupx_removed_fraction[colnames(obj)]

  if (s %in% rownames(subject_info)) {
    si <- subject_info[s, , drop = FALSE]
    for (cn in colnames(si)) {
      if (cn == "sample_name") next
      obj[[cn]] <- si[[cn]][1]
    }
  }

  # Compute RNA QC metrics used in filtering on the cleaned object
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  md <- obj@meta.data
  cells_keep_qc <- rownames(md)[
    apply_min(md$nCount_RNA, as.numeric(thr_row$min_counts)) &
      apply_max(md$nCount_RNA, as.numeric(thr_row$max_counts)) &
      apply_min(md$nFeature_RNA, as.numeric(thr_row$min_genes)) &
      apply_max(md$nFeature_RNA, as.numeric(thr_row$max_genes)) &
      apply_max(md$percent.mt, as.numeric(thr_row$max_pct_mt))
  ]

  n_before <- ncol(obj_raw)
  obj_f <- subset(obj, cells = cells_keep_qc)
  n_after_qc <- ncol(obj_f)

  # Detect RNA doublets on the post-QC object
  dbl_res <- run_scdblfinder_if_available(obj_f, sample_id = s)
  panel_path <- save_filtering_panel(obj_f, dbl_res$plot_score, sample_id = s)
  obj_f <- dbl_res$obj
  dbl_stats <- dbl_res$stats
  n_after_scdblfinder <- ncol(obj_f)

  cat("Cells before filter:", n_before, "\n")
  cat("Cells after QC filter:", n_after_qc, "\n")
  cat("Cells after scDblFinder:", n_after_scdblfinder, "\n")

  # Save the filtered sample object
  saveRDS(obj_f, out_rds)

  # Keep the object for the final merge
  objs_filt[[s]] <- obj_f

  # Log QC, SoupX, and scDblFinder results for this sample
  log_rows[[s]] <- data.frame(
    sample_id = s,
    status = "kept",
    exclusion_reason = NA_character_,
    n_before = n_before,
    n_after_qc = n_after_qc,
    n_after_scdblfinder = n_after_scdblfinder,
    frac_kept = round(n_after_scdblfinder / n_before, 4),
    soupx_status = soupx_stats$soupx_status,
    soupx_rho = soupx_stats$soupx_rho,
    soupx_umi_before = soupx_stats$soupx_umi_before,
    soupx_umi_after = soupx_stats$soupx_umi_after,
    soupx_umis_removed = soupx_stats$soupx_umis_removed,
    soupx_fraction_removed = soupx_stats$soupx_fraction_removed,
    soupx_median_genes_before = soupx_stats$soupx_median_genes_before,
    soupx_median_genes_after = soupx_stats$soupx_median_genes_after,
    scdblfinder_status = dbl_stats$scdblfinder_status,
    scdblfinder_error = dbl_stats$scdblfinder_error,
    scdblfinder_singlets = dbl_stats$scdblfinder_singlets,
    scdblfinder_doublets_removed = dbl_stats$scdblfinder_doublets_removed,
    panel_plot_path = panel_path
  )

  rm(obj_raw, obj, obj_f, counts, counts_rna)
  gc()
}

if (length(objs_filt) == 0) {
  stop("No filtered objects were produced")
}

# Save the filtering summary table
log_df <- do.call(rbind, log_rows)
write.csv(log_df, file.path(LOG_DIR, "gex_filtering_log.csv"), row.names = FALSE)
write.csv(log_df, file.path(POST_FILTER_QC_DIR, "filtering_summary.csv"), row.names = FALSE)
write.csv(log_df, file.path(FILTERING_QC_DIR, "filtering_summary.csv"), row.names = FALSE)

# Merge filtered objects into one Seurat object for downstream integration
cat("\nMerging filtered samples...\n")

merged_obj <- if (length(objs_filt) == 1) {
  objs_filt[[1]]
} else {
  merge(
    x = objs_filt[[1]],
    y = objs_filt[2:length(objs_filt)],
    add.cell.ids = names(objs_filt),
    project = "Thema_R_GEX"
  )
}

merged_path <- file.path(OBJECTS_DIR, "gex_filtered_samples.rds")
saveRDS(merged_obj, merged_path)

cat("\n============================================================\n")
cat("Filtering complete.\n")
cat("Merged object saved : ", merged_path, "\n")
cat("Per-sample objects  : ", OBJECTS_DIR, "\n")
cat("Filtering log       : ", file.path(LOG_DIR, "gex_filtering_log.csv"), "\n")
cat("Samples kept        : ", length(objs_filt), "\n")
cat("============================================================\n")
