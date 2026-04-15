# ---------------------------------------
# global_config.R
# Project-wide configuration file
# ---------------------------------------

# Project base directory
PROJECT_ROOT <- "/data/home/alessandro.digilio/Thema_R"

# Define key directories
DATA_DIR     <- file.path(PROJECT_ROOT, "data")
METADATA_DIR <- file.path(PROJECT_ROOT, "metadata")
SRC_DIR      <- file.path(PROJECT_ROOT, "src")
RESULTS_DIR  <- file.path(PROJECT_ROOT, "results")
QC_DIR       <- file.path(RESULTS_DIR, "qc")
FIGURES_DIR  <- file.path(PROJECT_ROOT, "figures")
LOGS_DIR     <- file.path(PROJECT_ROOT, "logs")

# Create directories if missing
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)

# GEX-specific directories
GEX_RAW_DATA_DIR    <- file.path(DATA_DIR, "raw_counts")
FILTERED_DATA_DIR   <- file.path(DATA_DIR, "filtered_data")
INTEGRATED_DATA_DIR <- file.path(DATA_DIR, "integrated_object")
POST_FILTER_QC_DIR  <- file.path(QC_DIR, "post_filtering")
FILTERING_QC_DIR    <- file.path(QC_DIR, "filtering")
FILTERING_LOG_DIR   <- file.path(RESULTS_DIR, "logs")

dir.create(FILTERED_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INTEGRATED_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(POST_FILTER_QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FILTERING_QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FILTERING_LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# Project-level sample exclusions
EXCLUDED_SAMPLES <- character(0)

# Optional cell-level filtering settings
RUN_SOUPX <- TRUE
RUN_SCDBLFINDER <- TRUE
SCDBLFINDER_DIMS <- 20

# Gene set files
FERROPTOSIS_GENESET_FILE <- file.path(METADATA_DIR, "iron_genes", "ferroptosis_genes_curated.xlsx")

# Integration settings
BATCH_VAR <- "sample_id"
N_PCS <- 25
HARMONY_DIMS <- 1:25
UMAP_NEIGHBORS <- 50
UMAP_MIN_DIST <- 0.4
UMAP_SPREAD <- 1.0
GRAPH_K_PARAM <- 30
CLUSTER_RES_GRID <- seq(0.2, 1.0, by = 0.1)
SELECTED_CLUSTER_RES <- 0.4

# Final manual cluster labels
GEX_CLUSTER_CELLTYPE <- c(
  "1"  = "Sublining fibroblasts (SFRP2+)",
  "2"  = "Lining fibroblasts (PRG4+)",
  "3"  = "Inflammatory fibroblasts (ADAM12+)",
  "4"  = "Activated / cytotoxic T cells",
  "5"  = "Remodeling fibroblasts (HTRA1+)",
  "6"  = "Inflammatory macrophages (IL1B+)",
  "7"  = "Activated endothelial cells",
  "8"  = "Destructive lining fibroblasts (MMP3+)",
  "9"  = "Resident macrophages (C1QC+)",
  "10" = "Pericytes / vascular smooth muscle cells",
  "11" = "Stress-response cells",
  "12" = "B cells",
  "13" = "cDC2",
  "14" = "Mixed-lineage cells"
)

# Final cell type colors for annotated UMAPs
cluster_name_colors <- c(
  "Sublining fibroblasts (SFRP2+)" = "#1F77B4",
  "Lining fibroblasts (PRG4+)" = "#D62728",
  "Inflammatory fibroblasts (ADAM12+)" = "#FF7F0E",
  "Activated / cytotoxic T cells" = "#FFBB78",
  "Remodeling fibroblasts (HTRA1+)" = "#2CA02C",
  "Inflammatory macrophages (IL1B+)" = "#98DF8A",
  "Activated endothelial cells" = "#AEC7E8",
  "Destructive lining fibroblasts (MMP3+)" = "#FF9896",
  "Resident macrophages (C1QC+)" = "#9467BD",
  "Pericytes / vascular smooth muscle cells" = "#C5B0D5",
  "Stress-response cells" = "#8C564B",
  "B cells" = "#C49C94",
  "cDC2" = "#9EDAE5",
  "Mixed-lineage cells" = "#F7B6D2"
)

# Three canonical markers per final cell type for the annotation dotplot
marker_genes <- list(
  "Sublining fibroblasts (SFRP2+)" = c("SFRP2", "COL1A1", "MFAP5"),
  "Lining fibroblasts (PRG4+)" = c("PRG4", "CLIC5", "DEFB1"),
  "Inflammatory fibroblasts (ADAM12+)" = c("ADAM12", "FGF2", "LIF"),
  "Activated / cytotoxic T cells" = c("GZMA", "BCL11B", "CD3D"),
  "Remodeling fibroblasts (HTRA1+)" = c("HTRA1", "CRTAC1", "NDUFA4L2"),
  "Inflammatory macrophages (IL1B+)" = c("IL1B", "IL1RN", "NLRP3"),
  "Activated endothelial cells" = c("SELE", "EMCN", "CLEC14A"),
  "Destructive lining fibroblasts (MMP3+)" = c("MMP3", "MMP1", "CXCL1"),
  "Resident macrophages (C1QC+)" = c("C1QC", "FOLR2", "TIMD4"),
  "Pericytes / vascular smooth muscle cells" = c("RGS5", "ACTA2", "MYH11"),
  "B cells" = c("IGKC", "CD79A", "IKZF3"),
  "Stress-response cells" = c("HSPA6", "HSPA1A", "DNAJB1"),
  "cDC2" = c("CD1C", "FCER1A", "CLEC10A"),
  "Mixed-lineage cells" = c("RNASE1", "C1QA", "CSN1S1")
)
