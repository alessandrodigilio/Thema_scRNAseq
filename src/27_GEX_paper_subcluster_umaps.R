# ---------------------------------------
# 28_GEX_paper_subcluster_umaps.R
# Publication-style UMAPs for reclustered cell compartments
# ---------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(readxl)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Save all publication panels in one dedicated folder
FIG_DIR <- file.path(FIGURES_DIR, "paper")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Plot settings for small publication panels
FIG_WIDTH <- 6
FIG_HEIGHT <- 6
FIG_DPI <- 600

PT_SIZE_MACROPHAGE <- 0.55
PT_SIZE_ENDOTHELIAL <- 0.70
PT_SIZE_FIBROBLAST <- 0.85

FERROPTOSIS_LOW_COLOR <- "#F5E6E8"
FERROPTOSIS_HIGH_COLOR <- "#7A1F2B"

# Keep full labels, but wrap long names in the legend
wrap_legend_label <- function(x) {
  vapply(strwrap(x, width = 34, simplify = FALSE), paste, collapse = "\n", FUN.VALUE = character(1))
}

# Common theme for all UMAPs
paper_umap_theme <- theme_classic(base_size = 18) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 18, color = "black"),
    axis.line = element_line(linewidth = 0.9, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid = element_blank(),
    aspect.ratio = 1,
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black", lineheight = 0.9),
    legend.position = "bottom",
    legend.key.size = grid::unit(0.48, "cm"),
    legend.spacing.y = grid::unit(0.10, "cm"),
    legend.box = "vertical",
    plot.margin = margin(10, 10, 10, 10)
  )


# ---------------------------------------
# 1. Destructive lining fibroblast subtypes
# ---------------------------------------

fibro_obj <- readRDS(file.path(DATA_DIR, "integrated_object", "gex_destructive_lining_fibroblasts_subclustered.rds"))

fibro_cluster_col <- "destructive_lining_fibroblast_subcluster"
fibro_label_col <- "destructive_lining_fibroblast_subtype"
fibro_reduction <- "umap.harmony.destructive.lining.fibroblast"

fibro_labels <- c(
  "0" = "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)",
  "1" = "Activated MMP3+ lining fibroblast cells (ID1+)",
  "2" = "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)",
  "3" = "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)",
  "4" = "MMP3+ lining fibroblast cells (FAM184A+)",
  "5" = "HA-enriched SFRP2+ matrix fibroblast-like cells"
)

fibro_colors <- c(
  "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)" = "#5F7EA6",
  "Activated MMP3+ lining fibroblast cells (ID1+)" = "#D98F5C",
  "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)" = "#C65A5A",
  "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)" = "#4C9F8A",
  "MMP3+ lining fibroblast cells (FAM184A+)" = "#9B7AAE",
  "HA-enriched SFRP2+ matrix fibroblast-like cells" = "#C7B24A"
)

if (!fibro_cluster_col %in% colnames(fibro_obj@meta.data)) {
  stop("Missing column: ", fibro_cluster_col)
}

if (!fibro_reduction %in% names(fibro_obj@reductions)) {
  stop("Missing reduction: ", fibro_reduction)
}

fibro_obj[[fibro_label_col]] <- unname(fibro_labels[as.character(fibro_obj[[fibro_cluster_col]][, 1])])
fibro_obj[[fibro_label_col]] <- factor(fibro_obj[[fibro_label_col]][, 1], levels = unname(fibro_labels))

p_fibro <- DimPlot(
  object = fibro_obj,
  reduction = fibro_reduction,
  group.by = fibro_label_col,
  raster = FALSE,
  pt.size = PT_SIZE_FIBROBLAST
) +
  scale_color_manual(values = fibro_colors, labels = wrap_legend_label, drop = FALSE) +
  coord_fixed() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  paper_umap_theme +
  guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 5)))

ggsave(
  filename = file.path(FIG_DIR, "paper_umap_destructive_lining_fibroblast_subtypes.png"),
  plot = p_fibro,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)

# Same fibroblast UMAP split by HA vs other, colored by fibroblast subtype
if ("condition" %in% colnames(fibro_obj@meta.data)) {
  fibro_obj$condition <- ifelse(as.character(fibro_obj$condition) == "HA", "HA", "other")
  fibro_obj$condition <- factor(fibro_obj$condition, levels = c("HA", "other"))

  p_fibro_ha_other <- DimPlot(
    object = fibro_obj,
    reduction = fibro_reduction,
    group.by = fibro_label_col,
    split.by = "condition",
    raster = FALSE,
    pt.size = PT_SIZE_FIBROBLAST
  ) +
    scale_color_manual(values = fibro_colors, labels = wrap_legend_label, drop = FALSE) +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2") +
    paper_umap_theme +
    theme(
      strip.text = element_text(size = 18, face = "bold", color = "black"),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 5)))

  ggsave(
    filename = file.path(FIG_DIR, "paper_umap_destructive_lining_fibroblast_subtypes_HA_vs_other.png"),
    plot = p_fibro_ha_other,
    width = 8,
    height = 8,
    dpi = FIG_DPI
  )
}


# ---------------------------------------
# 2. Macrophage subtypes
# ---------------------------------------

macro_obj <- readRDS(file.path(DATA_DIR, "integrated_object", "gex_macrophages_subclustered.rds"))

macro_cluster_col <- "macrophage_subcluster"
macro_label_col <- "macrophage_subtype"
macro_reduction <- "umap.harmony.macrophage"

macro_labels <- c(
  "0" = "Inflammatory macrophages (KANK1+)",
  "1" = "Inflammatory macrophages (THBS1+)",
  "2" = "Macrophage-like state (AMTN+)",
  "3" = "Resident macrophages (HSPA6+)",
  "4" = "Red-pulp-like resident macrophages (MERTK+)",
  "5" = "Mixed macrophage-like cells (RNASE1+)",
  "6" = "Plasma-like contaminants",
  "7" = "Low-confidence cells",
  "8" = "Proliferating macrophages"
)

macro_colors <- c(
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

if (!macro_cluster_col %in% colnames(macro_obj@meta.data)) {
  stop("Missing column: ", macro_cluster_col)
}

if (!macro_reduction %in% names(macro_obj@reductions)) {
  stop("Missing reduction: ", macro_reduction)
}

macro_obj[[macro_label_col]] <- unname(macro_labels[as.character(macro_obj[[macro_cluster_col]][, 1])])
macro_obj[[macro_label_col]] <- factor(macro_obj[[macro_label_col]][, 1], levels = unname(macro_labels))

p_macro <- DimPlot(
  object = macro_obj,
  reduction = macro_reduction,
  group.by = macro_label_col,
  raster = FALSE,
  pt.size = PT_SIZE_MACROPHAGE
) +
  scale_color_manual(values = macro_colors, labels = wrap_legend_label, drop = FALSE) +
  coord_fixed() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  paper_umap_theme +
  guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 5)))

ggsave(
  filename = file.path(FIG_DIR, "paper_umap_macrophage_subtypes.png"),
  plot = p_macro,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)


# ---------------------------------------
# 3. Activated endothelial subtypes
# ---------------------------------------

endo_obj <- readRDS(file.path(DATA_DIR, "integrated_object", "gex_activated_endothelial_subclustered.rds"))

endo_cluster_col <- "endothelial_subcluster"
endo_label_col <- "endothelial_subtype"
endo_reduction <- "umap.harmony.endothelial"

endo_labels <- c(
  "0" = "Endothelial cells (PLXNA4+)",
  "1" = "Stress-response endothelial cells (HSPA6+)",
  "2" = "Activated endothelial cells (IL6+)",
  "3" = "Endothelial cells (ZNF385B+)",
  "4" = "Endothelial cells (EDNRB+)",
  "5" = "Arterial-like endothelial cells (GJA5+)",
  "6" = "Endothelial cells (SLC2A14+)",
  "7" = "Mixed stromal-like cells",
  "8" = "Mural-like cells"
)

endo_colors <- c(
  "Endothelial cells (PLXNA4+)" = "#D46A6A",
  "Stress-response endothelial cells (HSPA6+)" = "#C7B24A",
  "Activated endothelial cells (IL6+)" = "#84B547",
  "Endothelial cells (ZNF385B+)" = "#7C7FD4",
  "Endothelial cells (EDNRB+)" = "#49B6A5",
  "Arterial-like endothelial cells (GJA5+)" = "#4E88C7",
  "Endothelial cells (SLC2A14+)" = "#57C6D9",
  "Mixed stromal-like cells" = "#C58A8A",
  "Mural-like cells" = "#C86DD7"
)

if (!endo_cluster_col %in% colnames(endo_obj@meta.data)) {
  stop("Missing column: ", endo_cluster_col)
}

if (!endo_reduction %in% names(endo_obj@reductions)) {
  stop("Missing reduction: ", endo_reduction)
}

endo_obj[[endo_label_col]] <- unname(endo_labels[as.character(endo_obj[[endo_cluster_col]][, 1])])
endo_obj[[endo_label_col]] <- factor(endo_obj[[endo_label_col]][, 1], levels = unname(endo_labels))

p_endo <- DimPlot(
  object = endo_obj,
  reduction = endo_reduction,
  group.by = endo_label_col,
  raster = FALSE,
  pt.size = PT_SIZE_ENDOTHELIAL
) +
  scale_color_manual(values = endo_colors, labels = wrap_legend_label, drop = FALSE) +
  coord_fixed() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  paper_umap_theme +
  guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 5)))

ggsave(
  filename = file.path(FIG_DIR, "paper_umap_endothelial_subtypes.png"),
  plot = p_endo,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)


# ---------------------------------------
# 4. Endothelial ferroptosis score: HA vs other
# ---------------------------------------

# Prepare RNA assay for module score calculation
DefaultAssay(endo_obj) <- "RNA"
endo_obj <- JoinLayers(endo_obj, assay = "RNA")

rna_data_layer <- tryCatch(
  LayerData(endo_obj[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  endo_obj <- NormalizeData(endo_obj, assay = "RNA", verbose = FALSE)
}

# Read the ferroptosis gene set used in the previous ferroptosis analysis
ferroptosis_genes_df <- as.data.frame(readxl::read_excel(FERROPTOSIS_GENESET_FILE))
ferroptosis_genes <- unique(trimws(as.character(ferroptosis_genes_df[[1]])))
ferroptosis_genes <- ferroptosis_genes[!is.na(ferroptosis_genes) & ferroptosis_genes != ""]
ferroptosis_genes[ferroptosis_genes == "TRFC"] <- "TFRC"
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(endo_obj))

if (length(ferroptosis_genes) == 0) {
  stop("No ferroptosis genes were found in the endothelial object")
}

# Add one score per cell from the ferroptosis gene set
endo_obj <- AddModuleScore(
  object = endo_obj,
  features = list(ferroptosis_genes),
  assay = "RNA",
  name = "Ferroptosis_Score",
  ctrl = 100,
  seed = 1234
)

endo_obj$ferroptosis_score <- endo_obj$Ferroptosis_Score1
endo_obj$HA_vs_other <- ifelse(as.character(endo_obj$condition) == "HA", "HA", "other")
endo_obj$HA_vs_other <- factor(endo_obj$HA_vs_other, levels = c("HA", "other"))

score_max <- max(endo_obj$ferroptosis_score, na.rm = TRUE)

# Plot ferroptosis score on the endothelial UMAP, split by condition
p_endo_ferroptosis_ha_other <- FeaturePlot(
  object = endo_obj,
  features = "ferroptosis_score",
  reduction = endo_reduction,
  split.by = "HA_vs_other",
  raster = FALSE,
  order = TRUE,
  pt.size = PT_SIZE_ENDOTHELIAL,
  cols = c(FERROPTOSIS_LOW_COLOR, FERROPTOSIS_HIGH_COLOR),
  min.cutoff = 0,
  max.cutoff = score_max
) +
  coord_fixed() +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Ferroptosis\nscore") +
  paper_umap_theme +
  theme(
    strip.text = element_text(size = 18, face = "bold", color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold", color = "black", angle = 90),
    legend.text = element_text(size = 12, color = "black")
  ) +
  guides(color = guide_colorbar(barheight = grid::unit(3.0, "cm"), barwidth = grid::unit(0.35, "cm")))

ggsave(
  filename = file.path(FIG_DIR, "paper_umap_endothelial_ferroptosis_score_HA_vs_other.png"),
  plot = p_endo_ferroptosis_ha_other,
  width = 8,
  height = 8,
  dpi = FIG_DPI
)

cat("Saved paper UMAPs in: ", FIG_DIR, "\n", sep = "")
