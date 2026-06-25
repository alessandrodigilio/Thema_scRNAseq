######################################################
### Publication UMAPs for reclustered compartments ###
######################################################

# prepare RNA assay for module score calculation
DefaultAssay(endo_obj) <- "RNA"
endo_obj <- JoinLayers(endo_obj, assay = "RNA")

rna_data_layer <- tryCatch(
  LayerData(endo_obj[["RNA"]], layer = "data"),
  error = function(e) NULL
)

if (is.null(rna_data_layer) || length(rna_data_layer@x) == 0) {
  endo_obj <- NormalizeData(endo_obj, assay = "RNA", verbose = FALSE)
}

# read the ferroptosis gene set used in the previous ferroptosis analysis
ferroptosis_genes_df <- as.data.frame(read_excel(ferroptosis_geneset_file))
ferroptosis_genes <- unique(trimws(as.character(ferroptosis_genes_df[[1]])))
ferroptosis_genes <- ferroptosis_genes[!is.na(ferroptosis_genes) & ferroptosis_genes != ""]
ferroptosis_genes[ferroptosis_genes == "TRFC"] <- "TFRC"
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(endo_obj))

if (length(ferroptosis_genes) == 0) {
  stop("No ferroptosis genes were found in the endothelial object")
}

# add one score per cell from the ferroptosis gene set
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

# plot ferroptosis score on the endothelial UMAP, split by condition
p_endo_ferroptosis_ha_other <- FeaturePlot(
  object = endo_obj,
  features = "ferroptosis_score",
  reduction = endo_reduction,
  split.by = "HA_vs_other",
  raster = FALSE,
  order = TRUE,
  pt.size = pt_size_endothelial,
  cols = c(ferroptosis_low_color, ferroptosis_high_color),
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
  filename = file.path(fig_dir, "paper_umap_endothelial_ferroptosis_score_HA_vs_other.png"),
  plot = p_endo_ferroptosis_ha_other,
  width = 8,
  height = 8,
  dpi = fig_dpi
)

cat("Saved paper UMAPs in: ", fig_dir, "\n", sep = "")
