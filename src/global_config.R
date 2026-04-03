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
FILTERING_LOG_DIR   <- file.path(RESULTS_DIR, "logs")

dir.create(FILTERED_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INTEGRATED_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(POST_FILTER_QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FILTERING_LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# Project-level sample exclusions
EXCLUDED_SAMPLES <- character(0)

# Optional cell-level filtering settings
RUN_SOUPX <- TRUE
RUN_SCDBLFINDER <- TRUE
SCDBLFINDER_DIMS <- 20

# Integration settings
BATCH_VAR <- "sample_id"
N_PCS <- 25
HARMONY_DIMS <- 1:25
UMAP_NEIGHBORS <- 50
UMAP_MIN_DIST <- 0.4
UMAP_SPREAD <- 1.0
GRAPH_K_PARAM <- 30
CLUSTER_RES_GRID <- seq(0.2, 1.0, by = 0.1)
SELECTED_CLUSTER_RES <- 0.3
