############################################
### Sample filtering and doublet removal ###
############################################

# load scRNA-seq samples listed in subject_info.xlsx, clean RNA with SoupX,
# filter cells by per-sample QC thresholds, remove RNA doublets with
# scDblFinder, and save filtered Seurat objects plus a merged object

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

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/atlas/utils.R") # helper functions

# directories
objects_dir <- filtered_data_dir
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
log_dir <- filtering_log_dir
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
fig_filtering_dir <- file.path(atlas_figures_dir, "filtering")
dir.create(fig_filtering_dir, recursive = TRUE, showWarnings = FALSE)

# load per-sample filtering thresholds
thr_path <- file.path(metadata_dir, "filters.xlsx")
if (!file.exists(thr_path)) stop("Missing thresholds Excel: ", thr_path)
thr <- as.data.frame(read_excel(thr_path))
colnames(thr) <- trimws(colnames(thr))
thr[] <- lapply(thr, function(x) if (is.character(x)) trimws(x) else x)

required_thr_cols <- c(
  "sample_name", "sample_path", "min_genes", "max_genes",
  "min_counts", "max_counts", "max_pct_mt"
)
# check col names
if (!all(required_thr_cols %in% colnames(thr))) {
  stop(
    "filters.xlsx is missing columns: ",
    paste(setdiff(required_thr_cols, colnames(thr)), collapse = ", ")
  )
}
thr$sample_name <- trimws(as.character(thr$sample_name))

# load sample annotations and define the only samples to process
subject_info_path <- file.path(metadata_dir, "subject_info.xlsx")
if (!file.exists(subject_info_path)) stop("Missing: ", subject_info_path)
# read and clean subject_info.xlsx
subject_info <- as.data.frame(read_excel(subject_info_path))
colnames(subject_info) <- trimws(colnames(subject_info))
subject_info[] <- lapply(subject_info, function(x) if (is.character(x)) trimws(x) else x)
# check
if (!"sample_name" %in% colnames(subject_info)) {
  stop("subject_info.xlsx needs a 'sample_name' column")
}
subject_info$sample_name <- trimws(as.character(subject_info$sample_name))
subject_info <- subject_info[subject_info$sample_name != "" & !is.na(subject_info$sample_name), , drop = FALSE]
subject_info <- subject_info[!duplicated(subject_info$sample_name), , drop = FALSE]
rownames(subject_info) <- subject_info$sample_name

samples <- subject_info$sample_name

cat("Samples in subject_info.xlsx:", length(samples), "\n")
cat("Samples to process:", length(samples), "\n")
cat("Sample IDs:", paste(samples, collapse = ", "), "\n")

# store filtered Seurat objects before the final merge
objs_filt <- list()

# accumulate one summary row per sample
log_rows <- list()

# keep only features present in all samples, matching the original filtering workflow
cat("\nComputing common RNA features across samples...\n")
feature_sets <- list()
for (s in samples) {
  thr_row <- thr[thr$sample_name == s, , drop = FALSE]
  if (nrow(thr_row) != 1) stop("Thresholds not found or duplicated for sample: ", s)

  features <- read_sample_features(
    sample_id = s,
    sample_path = thr_row$sample_path[1]
  )

  feature_sets[[s]] <- features
  cat(s, "features:", length(features), "\n")
}

common_features <- Reduce(intersect, feature_sets)
cat("Common features retained:", length(common_features), "\n")

for (s in samples) {
  cat("Filtering sample:", s, "\n")
  out_rds <- file.path(objects_dir, paste0(s, "_filtered.rds"))
  thr_row <- thr[thr$sample_name == s, , drop = FALSE]
  if (nrow(thr_row) != 1) stop("Thresholds not found or duplicated for sample: ", s)
  counts <- read_sample_counts(
    sample_id = s,
    sample_path = thr_row$sample_path[1]
  )
  counts <- counts[common_features, , drop = FALSE]

  # build a temp obj to compute RNA QC and define filtered cells
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

  # clean RNA ambient signal before building the final object
  soupx_res <- run_soupx_correction(
    sample_id = s,
    raw_counts = counts,
    filtered_counts = filtered_counts
  )
  counts_rna <- soupx_res$counts
  soupx_removed_umis <- soupx_res$removed_umis
  soupx_removed_fraction <- soupx_res$removed_fraction
  soupx_stats <- soupx_res$stats

  # create the filtered RNA object and attach sample metadata
  obj <- CreateSeuratObject(
    counts = counts_rna,
    assay = "RNA",
    project = s
  )
  obj$sample_id <- s
  obj$soupx_removed_umis <- soupx_removed_umis[colnames(obj)]
  obj$soupx_removed_fraction <- soupx_removed_fraction[colnames(obj)]

  si <- subject_info[s, , drop = FALSE]
  for (cn in colnames(si)) {
    if (cn == "sample_name") next
    obj[[cn]] <- si[[cn]][1]
  }

  # compute RNA QC metrics used in filtering on the cleaned object
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

  # detect RNA doublets on the post-QC object
  dbl_res <- run_scdblfinder_filter(obj_f, sample_id = s)
  panel_path <- save_filtering_panel(obj_f, dbl_res$plot_score, sample_id = s)
  obj_f <- dbl_res$obj
  dbl_stats <- dbl_res$stats
  n_after_scdblfinder <- ncol(obj_f)

  cat("Cells before filter:", n_before, "\n")
  cat("Cells after QC filter:", n_after_qc, "\n")
  cat("Cells after scDblFinder:", n_after_scdblfinder, "\n")

  # save the filtered sample object
  saveRDS(obj_f, out_rds)

  # keep the object for the final merge
  objs_filt[[s]] <- obj_f

  # log QC, SoupX, and scDblFinder results for this sample
  log_rows[[s]] <- data.frame(
    sample_id = s,
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
    scdblfinder_singlets = dbl_stats$scdblfinder_singlets,
    scdblfinder_doublets_removed = dbl_stats$scdblfinder_doublets_removed,
    panel_plot_path = panel_path
  )

  rm(obj_raw, obj, obj_f, counts, counts_rna)
  gc()
}

# save the filtering summary table
log_df <- do.call(rbind, log_rows)
write.csv(log_df, file.path(log_dir, "filtering_log.csv"), row.names = FALSE)
write.csv(log_df, file.path(post_filter_qc_dir, "filtering_summary.csv"), row.names = FALSE)
write.csv(log_df, file.path(filtering_qc_dir, "filtering_summary.csv"), row.names = FALSE)

# merge filtered objects into one Seurat object for downstream integration
cat("\nMerging filtered samples...\n")

merged_obj <- if (length(objs_filt) == 1) {
  objs_filt[[1]]
} else {
  merge(
    x = objs_filt[[1]],
    y = objs_filt[2:length(objs_filt)],
    add.cell.ids = names(objs_filt),
    project = "Thema_R"
  )
}

merged_path <- file.path(objects_dir, "filtered_samples.rds")
saveRDS(merged_obj, merged_path)

cat("\n============================================================\n")
cat("Filtering complete.\n")
cat("Merged object saved : ", merged_path, "\n")
cat("Per-sample objects  : ", objects_dir, "\n")
cat("Filtering log       : ", file.path(log_dir, "filtering_log.csv"), "\n")
cat("Samples kept        : ", length(objs_filt), "\n")
cat("============================================================\n")
