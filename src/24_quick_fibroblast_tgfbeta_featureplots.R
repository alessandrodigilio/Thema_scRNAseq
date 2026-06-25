##############################################################
### Quick TGF-beta plots in destructive lining fibroblasts ###
##############################################################

# quick FeaturePlot of TGF-beta-related genes in destructive lining fibroblast subclusters.
# this script is intended for fast exploratory plotting from console.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

setwd("~/Thema_R")
source("src/global_config.R")

# input and output
input_object <- file.path(data_dir, "integrated_object", "destructive_lining_fibroblasts_subclustered.rds")
fig_dir <- file.path(figures_dir, "quick_fibroblast_TGFbeta_featureplots")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
reduction_name <- "umap.harmony.destructive.lining.fibroblast"
pt_size <- 0.35
n_cols <- 4

# edit this list to test other TGF-beta-related genes.
genes_to_plot <- c(
  "TGFB1", "TGFB2", "TGFB3",
  "TGFBR1", "TGFBR2", "TGFBR3",
  "SMAD2", "SMAD3", "SMAD4", "SMAD7",
  "SERPINE1", "CTGF", "COL1A1", "COL3A1"
)

output_name <- paste0(
  "featureplot_destructive_lining_fibroblasts_TGFbeta_",
  paste(genes_to_plot, collapse = "_"),
  ".png"
)

# load object
if (!file.exists(input_object)) {
  stop("Missing destructive lining fibroblast object: ", input_object)
}

cat("Loading destructive lining fibroblast subcluster object...\n")
obj_fibro <- readRDS(input_object)
cat("Destructive lining fibroblast cells loaded:", ncol(obj_fibro), "\n")

if (!reduction_name %in% names(obj_fibro@reductions)) {
  cat("Requested reduction not found:", reduction_name, "\n")
  cat("Available reductions:\n")
  print(names(obj_fibro@reductions))

  umap_candidates <- grep("umap", names(obj_fibro@reductions), ignore.case = TRUE, value = TRUE)
  if (length(umap_candidates) == 0) {
    stop("No UMAP-like reduction found in object.")
  }

  reduction_name <- umap_candidates[1]
  cat("Using fallback reduction:", reduction_name, "\n")
}

# make sure RNA normalized data are available.
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

# check selected genes.
genes_use <- genes_to_plot[genes_to_plot %in% rownames(obj_fibro)]
missing_genes <- setdiff(genes_to_plot, genes_use)

cat("Genes found:\n")
print(genes_use)

if (length(missing_genes) > 0) {
  cat("Genes not found:\n")
  print(missing_genes)
}

if (length(genes_use) == 0) {
  stop("None of the selected genes were found in the fibroblast object.")
}

# featurePlot on the destructive lining fibroblast subcluster UMAP.
p_feature <- FeaturePlot(
  object = obj_fibro,
  features = genes_use,
  reduction = reduction_name,
  cols = c("#F3F1EC", "#7A1F2B"),
  order = TRUE,
  raster = FALSE,
  pt.size = pt_size,
  ncol = n_cols
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

n_rows <- ceiling(length(genes_use) / n_cols)
output_file <- file.path(fig_dir, output_name)

ggsave(
  filename = output_file,
  plot = p_feature,
  width = 18,
  height = 4.2 * n_rows,
  dpi = 600
)

cat("Saved plot:\n")
cat(output_file, "\n")
