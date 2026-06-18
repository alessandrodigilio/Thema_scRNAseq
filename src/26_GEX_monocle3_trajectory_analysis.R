# ---------------------------------------
# 26_GEX_monocle3_trajectory_analysis.R
# Monocle3 trajectory and pseudotime analysis
# ---------------------------------------

# Goal:
# Run Monocle3 trajectory analysis on already reclustered compartments.
#
# Important:
# The macrophage, endothelial and destructive lining fibroblast labels were
# generated from subset-specific Seurat/Harmony UMAPs in scripts 11, 16 and 21.
# Therefore this script follows the tutorial logic:
# - use the existing Seurat subset UMAP
# - use the existing Seurat subclusters/states
# - let Monocle learn only the trajectory graph and pseudotime
#
# We do NOT recalculate PCA/UMAP in Monocle here, because that would create
# a new geometry different from the one used for the manual annotation.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(ggplot2)
  library(patchwork)
})

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
OUT_DIR <- file.path(DATA_DIR, "integrated_object")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(FIGURES_DIR, "monocle3_trajectory_analysis")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

RES_DIR <- file.path(RESULTS_DIR, "monocle3_trajectory_analysis")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# General parameters
ASSAY_USE <- "RNA"
MIN_CELLS <- 200
RUN_GRAPH_TEST <- TRUE
GRAPH_TEST_CORES <- 4
N_GRAPH_TEST_GENES <- 2000
PT_SIZE <- 0.65

set.seed(1234)

# Palettes copied from the corresponding final annotation scripts.
macrophage_subtype_colors <- c(
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

endothelial_subtype_colors <- c(
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

destructive_lining_fibroblast_subtype_colors <- c(
  "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)" = "#5F7EA6",
  "Activated MMP3+ lining fibroblast cells (ID1+)" = "#D98F5C",
  "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)" = "#C65A5A",
  "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)" = "#4C9F8A",
  "MMP3+ lining fibroblast cells (FAM184A+)" = "#9B7AAE",
  "HA-enriched SFRP2+ matrix fibroblast-like cells" = "#C7B24A"
)

condition_colors <- c(
  "other" = "#B65A5A",
  "HA" = "#5B8DB8"
)

build_label_colors <- function(labels, color_map) {
  labels <- unique(as.character(labels))
  labels <- labels[!is.na(labels)]
  cols <- color_map[labels]
  missing_labels <- labels[is.na(cols)]

  if (length(missing_labels) > 0) {
    fallback <- grDevices::hcl.colors(length(missing_labels), palette = "Set 3")
    names(fallback) <- missing_labels
    cols[names(fallback)] <- fallback
  }

  cols[labels]
}

# ---------------------------------------
# Define the trajectories to run
# ---------------------------------------

# These are the same subset objects created during subclustering.
# They contain the subset-specific UMAP reductions:
# - umap.harmony.macrophage
# - umap.harmony.endothelial
# - umap.harmony.destructive.lining.fibroblast

trajectory_list <- list(
  destructive_lining_fibroblast_states = list(
    input_object = file.path(DATA_DIR, "integrated_object", "gex_destructive_lining_fibroblasts_subclustered.rds"),
    reduction_name = "umap.harmony.destructive.lining.fibroblast",
    cluster_col = "destructive_lining_fibroblast_subcluster",
    state_col = "destructive_lining_fibroblast_subtype",
    state_colors = destructive_lining_fibroblast_subtype_colors,
    cluster_to_state = c(
      "0" = "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)",
      "1" = "Activated MMP3+ lining fibroblast cells (ID1+)",
      "2" = "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)",
      "3" = "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)",
      "4" = "MMP3+ lining fibroblast cells (FAM184A+)",
      "5" = "HA-enriched SFRP2+ matrix fibroblast-like cells"
    ),
    states_keep = c(
      "HLA-II MMP3+ lining fibroblasts (HLA-DRA+)",
      "Activated MMP3+ lining fibroblast cells (ID1+)",
      "HA-enriched inflammatory MMP3+ lining fibroblasts (CCL7+/CXCL1+)",
      "Matrix-adhesion MMP3+ lining fibroblast cells (ITGB8+)",
      "MMP3+ lining fibroblast cells (FAM184A+)",
      "HA-enriched SFRP2+ matrix fibroblast-like cells"
    ),
    root_states = c(
      "MMP3+ lining fibroblast cells (FAM184A+)",
      "Activated MMP3+ lining fibroblast cells (ID1+)"
    ),
    marker_genes = c(
      "FAM184A", "ID1", "HLA-DRA", "CD74", "CCL7", "CXCL1",
      "CCL20", "ITGB8", "SFRP2", "MMP3", "MMP1", "HMOX1", "SLC40A1"
    )
  ),

  macrophage_states = list(
    input_object = file.path(DATA_DIR, "integrated_object", "gex_macrophages_subclustered.rds"),
    reduction_name = "umap.harmony.macrophage",
    cluster_col = "macrophage_subcluster",
    state_col = "macrophage_subtype",
    state_colors = macrophage_subtype_colors,
    cluster_to_state = c(
      "0" = "Inflammatory macrophages (KANK1+)",
      "1" = "Inflammatory macrophages (THBS1+)",
      "2" = "Macrophage-like state (AMTN+)",
      "3" = "Resident macrophages (HSPA6+)",
      "4" = "Red-pulp-like resident macrophages (MERTK+)",
      "5" = "Mixed macrophage-like cells (RNASE1+)",
      "6" = "Plasma-like contaminants",
      "7" = "Low-confidence cells",
      "8" = "Proliferating macrophages"
    ),
    states_keep = c(
      "Inflammatory macrophages (KANK1+)",
      "Inflammatory macrophages (THBS1+)",
      "Macrophage-like state (AMTN+)",
      "Resident macrophages (HSPA6+)",
      "Red-pulp-like resident macrophages (MERTK+)",
      "Proliferating macrophages"
    ),
    root_states = c(
      "Red-pulp-like resident macrophages (MERTK+)",
      "Resident macrophages (HSPA6+)"
    ),
    marker_genes = c(
      "MERTK", "CD163", "FCGR3A", "SPIC", "HMOX1", "SLC40A1",
      "KANK1", "THBS1", "IL1B", "AMTN", "SULF1", "MKI67"
    )
  ),

  endothelial_states = list(
    input_object = file.path(DATA_DIR, "integrated_object", "gex_activated_endothelial_subclustered.rds"),
    reduction_name = "umap.harmony.endothelial",
    cluster_col = "endothelial_subcluster",
    state_col = "endothelial_subtype",
    state_colors = endothelial_subtype_colors,
    cluster_to_state = c(
      "0" = "Endothelial cells (PLXNA4+)",
      "1" = "Stress-response endothelial cells (HSPA6+)",
      "2" = "Activated endothelial cells (IL6+)",
      "3" = "Endothelial cells (ZNF385B+)",
      "4" = "Endothelial cells (EDNRB+)",
      "5" = "Arterial-like endothelial cells (GJA5+)",
      "6" = "Endothelial cells (SLC2A14+)",
      "7" = "Mixed stromal-like cells",
      "8" = "Mural-like cells"
    ),
    states_keep = c(
      "Endothelial cells (PLXNA4+)",
      "Stress-response endothelial cells (HSPA6+)",
      "Activated endothelial cells (IL6+)",
      "Endothelial cells (ZNF385B+)",
      "Endothelial cells (EDNRB+)",
      "Arterial-like endothelial cells (GJA5+)",
      "Endothelial cells (SLC2A14+)"
    ),
    root_states = c(
      "Endothelial cells (PLXNA4+)",
      "Endothelial cells (ZNF385B+)",
      "Endothelial cells (EDNRB+)"
    ),
    marker_genes = c(
      "PLXNA4", "ZNF385B", "EDNRB", "GJA5", "SLC2A14", "IL6",
      "RGS16", "HES1", "HSPA6", "MMP1", "HMOX1", "SLC40A1"
    )
  )
)

# ---------------------------------------
# Run each trajectory
# ---------------------------------------

run_summary <- data.frame()

for (analysis_id in names(trajectory_list)) {
  cat("\n============================================================\n")
  cat("Running trajectory:", analysis_id, "\n")
  cat("============================================================\n")

  params <- trajectory_list[[analysis_id]]

  # 1. Load the already subclustered Seurat subset object.
  if (!file.exists(params$input_object)) {
    warning("Missing input object, skipping: ", params$input_object)
    next
  }

  obj <- readRDS(params$input_object)
  DefaultAssay(obj) <- ASSAY_USE

  # 2. Check that the expected subcluster column and UMAP reduction exist.
  if (!params$cluster_col %in% colnames(obj@meta.data)) {
    stop("Metadata column not found: ", params$cluster_col)
  }

  if (!params$reduction_name %in% names(obj@reductions)) {
    stop("Reduction not found: ", params$reduction_name)
  }

  # 3. Add final state labels directly to the subset object.
  #    The subset objects contain subcluster IDs. The final scripts mapped
  #    those IDs to biological labels; we repeat that simple mapping here.
  cluster_ids <- as.character(obj@meta.data[[params$cluster_col]])
  state_labels <- unname(params$cluster_to_state[cluster_ids])
  state_labels[is.na(state_labels)] <- "Unknown"
  obj@meta.data[[params$state_col]] <- state_labels

  # 4. Keep only the states relevant for trajectory.
  #    This excludes contaminant, mixed or low-confidence states where defined.
  cells_keep <- rownames(obj@meta.data)[
    as.character(obj@meta.data[[params$state_col]]) %in% params$states_keep
  ]

  if (length(cells_keep) < MIN_CELLS) {
    warning("Too few cells for ", analysis_id, ": ", length(cells_keep))
    next
  }

  obj <- subset(obj, cells = cells_keep)

  # 5. Make sure the RNA data layer is usable.
  obj <- tryCatch(
    JoinLayers(obj, assay = ASSAY_USE),
    error = function(e) obj
  )

  rna_data <- tryCatch(
    LayerData(obj[[ASSAY_USE]], layer = "data"),
    error = function(e) NULL
  )

  if (is.null(rna_data) || length(rna_data) == 0) {
    obj <- NormalizeData(obj, assay = ASSAY_USE, verbose = FALSE)
  }

  # 6. Convert Seurat object to Monocle3 cell_data_set.
  cds <- as.cell_data_set(obj)
  rowData(cds)$gene_short_name <- rownames(cds)

  # 7. Copy the existing Seurat subset UMAP into Monocle.
  #    This is the key point: we keep the same geometry used for annotation.
  seurat_umap <- Embeddings(obj, reduction = params$reduction_name)
  seurat_umap <- seurat_umap[colnames(cds), , drop = FALSE]
  colnames(seurat_umap) <- c("UMAP_1", "UMAP_2")
  cds@int_colData@listData$reducedDims$UMAP <- seurat_umap

  # 8. Give Monocle cluster/partition information.
  #    One partition means one shared trajectory in the selected compartment.
  monocle_partition <- factor(rep(1, ncol(cds)))
  names(monocle_partition) <- colnames(cds)
  cds@clusters$UMAP$partitions <- monocle_partition

  monocle_clusters <- factor(as.character(colData(cds)[[params$state_col]]))
  names(monocle_clusters) <- colnames(cds)
  cds@clusters$UMAP$clusters <- monocle_clusters

  # 9. Learn the trajectory graph on the existing UMAP.
  cds <- learn_graph(
    cds,
    use_partition = FALSE,
    verbose = FALSE
  )

  # 10. Choose root cells from the biologically selected starting states.
  root_cells <- rownames(obj@meta.data)[
    as.character(obj@meta.data[[params$state_col]]) %in% params$root_states
  ]

  root_cells <- intersect(root_cells, colnames(cds))

  if (length(root_cells) == 0) {
    stop("No root cells found for: ", analysis_id)
  }

  # 11. Order cells in pseudotime.
  cds <- order_cells(
    cds,
    reduction_method = "UMAP",
    root_cells = root_cells
  )

  # 12. Extract pseudotime and save it with metadata.
  pseudotime_values <- pseudotime(cds)
  pseudotime_values[!is.finite(pseudotime_values)] <- NA_real_
  colData(cds)$pseudotime <- pseudotime_values[colnames(cds)]

  cell_metadata <- as.data.frame(colData(cds))
  cell_metadata$cell_id <- rownames(cell_metadata)
  cell_metadata$pseudotime <- as.numeric(cell_metadata$pseudotime)
  cell_metadata$is_root_cell <- cell_metadata$cell_id %in% root_cells

  write.csv(
    cell_metadata,
    file.path(RES_DIR, paste0(analysis_id, "_cell_pseudotime_metadata.csv")),
    row.names = FALSE
  )

  # 13. Summarize pseudotime by condition.
  if ("condition" %in% colnames(cell_metadata)) {
    condition_summary <- aggregate(
      pseudotime ~ condition,
      data = cell_metadata,
      FUN = function(x) {
        c(
          n_cells = length(x),
          median = median(x, na.rm = TRUE),
          mean = mean(x, na.rm = TRUE)
        )
      }
    )

    condition_summary <- do.call(data.frame, condition_summary)

    write.csv(
      condition_summary,
      file.path(RES_DIR, paste0(analysis_id, "_pseudotime_summary_by_condition.csv")),
      row.names = FALSE
    )

    condition_values <- unique(as.character(cell_metadata$condition))
    condition_values <- condition_values[!is.na(condition_values)]

    if (all(c("HA", "other") %in% condition_values)) {
      test_df <- cell_metadata[
        cell_metadata$condition %in% c("HA", "other") &
          !is.na(cell_metadata$pseudotime),
        ,
        drop = FALSE
      ]

      wilcox_res <- wilcox.test(pseudotime ~ condition, data = test_df)

      wilcox_table <- data.frame(
        analysis_id = analysis_id,
        n_HA = sum(test_df$condition == "HA"),
        n_other = sum(test_df$condition == "other"),
        median_HA = median(test_df$pseudotime[test_df$condition == "HA"]),
        median_other = median(test_df$pseudotime[test_df$condition == "other"]),
        wilcox_pvalue = wilcox_res$p.value
      )

      write.csv(
        wilcox_table,
        file.path(RES_DIR, paste0(analysis_id, "_pseudotime_HA_vs_other_wilcox.csv")),
        row.names = FALSE
      )
    }
  }

  # 14. Summarize pseudotime by sample.
  if (all(c("sample_id", "condition") %in% colnames(cell_metadata))) {
    sample_summary <- aggregate(
      pseudotime ~ sample_id + condition,
      data = cell_metadata,
      FUN = function(x) {
        c(
          n_cells = length(x),
          median = median(x, na.rm = TRUE),
          mean = mean(x, na.rm = TRUE)
        )
      }
    )

    sample_summary <- do.call(data.frame, sample_summary)

    write.csv(
      sample_summary,
      file.path(RES_DIR, paste0(analysis_id, "_pseudotime_summary_by_sample.csv")),
      row.names = FALSE
    )
  }

  # 15. Plot trajectory colored by pseudotime, state, condition and sample.
  #     We do not save Monocle or Seurat objects here. The CSV files contain
  #     the pseudotime values needed for downstream interpretation.
  state_color_map <- build_label_colors(unique(cell_metadata[[params$state_col]]), params$state_colors)

  p_pseudotime <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_roots = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 1.5,
    cell_size = PT_SIZE
  ) +
    ggtitle(paste0(analysis_id, ": pseudotime")) +
    coord_equal() +
    theme_classic(base_size = 14)

  p_state <- plot_cells(
    cds,
    color_cells_by = params$state_col,
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_roots = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 1.5,
    cell_size = PT_SIZE
  ) +
    ggtitle(paste0(analysis_id, ": states")) +
    scale_color_manual(values = state_color_map, drop = FALSE) +
    coord_equal() +
    theme_classic(base_size = 14)

  plot_list <- list(p_pseudotime, p_state)

  if ("condition" %in% colnames(cell_metadata)) {
    condition_values <- unique(as.character(cell_metadata$condition))
    condition_color_map <- build_label_colors(condition_values, condition_colors)

    p_condition <- plot_cells(
      cds,
      color_cells_by = "condition",
      show_trajectory_graph = TRUE,
      label_cell_groups = FALSE,
      label_groups_by_cluster = FALSE,
      label_roots = FALSE,
      label_leaves = FALSE,
      label_branch_points = FALSE,
      graph_label_size = 1.5,
      cell_size = PT_SIZE
    ) +
      ggtitle(paste0(analysis_id, ": condition")) +
      scale_color_manual(values = condition_color_map, drop = FALSE) +
      coord_equal() +
      theme_classic(base_size = 14)

    plot_list[[length(plot_list) + 1]] <- p_condition
  }

  if ("sample_id" %in% colnames(cell_metadata)) {
    sample_values <- sort(unique(as.character(cell_metadata$sample_id)))
    sample_values <- sample_values[!is.na(sample_values)]
    sample_color_map <- grDevices::hcl.colors(length(sample_values), palette = "Dark 3")
    names(sample_color_map) <- sample_values

    p_sample <- plot_cells(
      cds,
      color_cells_by = "sample_id",
      show_trajectory_graph = TRUE,
      label_cell_groups = FALSE,
      label_groups_by_cluster = FALSE,
      label_roots = FALSE,
      label_leaves = FALSE,
      label_branch_points = FALSE,
      graph_label_size = 1.5,
      cell_size = PT_SIZE
    ) +
      ggtitle(paste0(analysis_id, ": sample_id")) +
      scale_color_manual(values = sample_color_map, drop = FALSE) +
      coord_equal() +
      theme_classic(base_size = 14)

    plot_list[[length(plot_list) + 1]] <- p_sample
  }

  p_review <- wrap_plots(plot_list, ncol = 2)

  ggsave(
    file.path(FIG_DIR, paste0(analysis_id, "_trajectory_review.png")),
    p_review,
    width = 16,
    height = 12,
    dpi = 600
  )

  # 16. Plot where each cell state falls along pseudotime.
  #     This is easier to read than marker plots for checking which states are
  #     early or late in the Monocle ordering.
  boxplot_df <- cell_metadata[!is.na(cell_metadata$pseudotime), ]
  boxplot_df$state_for_plot <- boxplot_df[[params$state_col]]
  boxplot_df <- boxplot_df[!is.na(boxplot_df$state_for_plot), ]

  if (nrow(boxplot_df) > 0) {
    n_states <- length(unique(boxplot_df$state_for_plot))

    p_state_boxplot <- ggplot(
      boxplot_df,
      aes(
        x = pseudotime,
        y = reorder(state_for_plot, pseudotime, median),
        fill = state_for_plot
      )
    ) +
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
      scale_fill_manual(values = state_color_map, drop = FALSE) +
      labs(
        title = paste0(analysis_id, ": pseudotime by state"),
        x = "pseudotime",
        y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")

    ggsave(
      file.path(FIG_DIR, paste0(analysis_id, "_pseudotime_by_state_boxplot.png")),
      p_state_boxplot,
      width = 10,
      height = max(4, n_states * 0.6),
      dpi = 600
    )
  }

  # 17. Run graph_test to find genes associated with the trajectory.
  marker_genes <- params$marker_genes[params$marker_genes %in% rownames(cds)]

  # 18. Optional graph_test.
  #     Use the variable genes already stored in the Seurat subset object.
  if (RUN_GRAPH_TEST) {
    graph_genes <- VariableFeatures(obj)
    graph_genes <- intersect(graph_genes, rownames(cds))
    graph_genes <- head(graph_genes, N_GRAPH_TEST_GENES)
    graph_genes <- unique(c(marker_genes, graph_genes))

    if (length(graph_genes) > 0) {
      graph_res <- graph_test(
        cds[graph_genes, ],
        neighbor_graph = "principal_graph",
        cores = GRAPH_TEST_CORES
      )

      graph_res$gene <- rownames(graph_res)
      graph_res <- graph_res[order(graph_res$q_value, -graph_res$morans_I), ]

      write.csv(
        graph_res,
        file.path(RES_DIR, paste0(analysis_id, "_graph_test_pseudotime_genes.csv")),
        row.names = FALSE
      )
    }
  }

  run_summary <- rbind(
    run_summary,
    data.frame(
      analysis_id = analysis_id,
      input_object = params$input_object,
      reduction_name = params$reduction_name,
      cluster_col = params$cluster_col,
      state_col = params$state_col,
      n_cells = ncol(obj),
      n_root_cells = length(root_cells),
      n_genes = nrow(cds),
      stringsAsFactors = FALSE
    )
  )
}

# Save a small table summarizing which analyses were completed.
write.csv(
  run_summary,
  file.path(RES_DIR, "monocle3_trajectory_run_summary.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Monocle3 trajectory analysis complete.\n")
cat("Results dir  : ", RES_DIR, "\n")
cat("Figures dir  : ", FIG_DIR, "\n")
cat("Objects dir  : ", OUT_DIR, "\n")
cat("============================================================\n")
