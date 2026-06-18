# ===============================================================
#   25_GEX_quick_fibroblast_TGFbeta_featureplots.R
# ===============================================================

# Quick FeaturePlot of TGF-beta-related genes in destructive lining fibroblast subclusters.
# This script is intended for fast exploratory plotting from console.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Input and output
INPUT_OBJECT <- file.path(DATA_DIR, "integrated_object", "gex_destructive_lining_fibroblasts_subclustered.rds")
FIG_DIR <- file.path(FIGURES_DIR, "quick_fibroblast_TGFbeta_featureplots")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Set parameters
REDUCTION_NAME <- "umap.harmony.destructive.lining.fibroblast"
PT_SIZE <- 0.35
N_COLS <- 4

# Edit this list to test other TGF-beta-related genes.
GENES_TO_PLOT <- c(
  "TGFB1", "TGFB2", "TGFB3",
  "TGFBR1", "TGFBR2", "TGFBR3",
  "SMAD2", "SMAD3", "SMAD4", "SMAD7",
  "SERPINE1", "CTGF", "COL1A1", "COL3A1"
)

OUTPUT_NAME <- paste0(
  "featureplot_destructive_lining_fibroblasts_TGFbeta_",
  paste(GENES_TO_PLOT, collapse = "_"),
  ".png"
)

# Load object
if (!file.exists(INPUT_OBJECT)) {
  stop("Missing destructive lining fibroblast object: ", INPUT_OBJECT)
}

cat("Loading destructive lining fibroblast subcluster object...\n")
obj_fibro <- readRDS(INPUT_OBJECT)
cat("Destructive lining fibroblast cells loaded:", ncol(obj_fibro), "\n")

if (!REDUCTION_NAME %in% names(obj_fibro@reductions)) {
  cat("Requested reduction not found:", REDUCTION_NAME, "\n")
  cat("Available reductions:\n")
  print(names(obj_fibro@reductions))

  umap_candidates <- grep("umap", names(obj_fibro@reductions), ignore.case = TRUE, value = TRUE)
  if (length(umap_candidates) == 0) {
    stop("No UMAP-like reduction found in object.")
  }

  REDUCTION_NAME <- umap_candidates[1]
  cat("Using fallback reduction:", REDUCTION_NAME, "\n")
}

# Make sure RNA normalized data are available.
DefaultAssay(obj_fibro) <- "RNA"
obj_fibro <- tryCatch(JoinLayers(obj_fibro, assay = "RNA"), error = function(e) obj_fibro)

rna_data_layer <- tryCatch(
  LayerData(obj_fibro[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  cat("Normalizing RNA assay...\n")
  obj_fibro <- NormalizeData(obj_fibro, assay = "RNA", verbose = FALSE)
}

# Check selected genes.
genes_use <- GENES_TO_PLOT[GENES_TO_PLOT %in% rownames(obj_fibro)]
missing_genes <- setdiff(GENES_TO_PLOT, genes_use)

cat("Genes found:\n")
print(genes_use)

if (length(missing_genes) > 0) {
  cat("Genes not found:\n")
  print(missing_genes)
}

if (length(genes_use) == 0) {
  stop("None of the selected genes were found in the fibroblast object.")
}

# FeaturePlot on the destructive lining fibroblast subcluster UMAP.
p_feature <- FeaturePlot(
  object = obj_fibro,
  features = genes_use,
  reduction = REDUCTION_NAME,
  cols = c("#F3F1EC", "#7A1F2B"),
  order = TRUE,
  raster = FALSE,
  pt.size = PT_SIZE,
  ncol = N_COLS
) &
  theme_classic(base_size = 18) &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(16, 16, 16, 16)
  )

print(p_feature)

n_rows <- ceiling(length(genes_use) / N_COLS)
output_file <- file.path(FIG_DIR, OUTPUT_NAME)

ggsave(
  filename = output_file,
  plot = p_feature,
  width = 18,
  height = 4.2 * n_rows,
  dpi = 600
)

cat("Saved plot:\n")
cat(output_file, "\n")
