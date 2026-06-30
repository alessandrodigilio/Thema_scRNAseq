###################################################
### Monocle3 trajectory and pseudotime analysis ###
###################################################

# run Monocle3 on reclustered macrophage, endothelial and
# destructive lining fibroblast compartments, using their existing Seurat UMAPs

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(ggplot2)
  library(patchwork)
})

# wd
setwd("~/Thema_R")
source("src/global_config.R")
source("src/trajectory_analysis/utils.R")

# directories
fig_dir <- file.path(figures_dir, "monocle3_trajectory_analysis")
res_dir <- file.path(results_dir, "monocle3_trajectory_analysis")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# settings
assay_use <- "RNA"
min_cells <- 200
run_graph_test <- TRUE
graph_test_cores <- 4
n_graph_test_genes <- 2000
pt_size <- 0.65

condition_colors <- c(
  "other" = "#B65A5A",
  "HA" = "#5B8DB8"
)

set.seed(1234)

# define the trajectories to run
trajectory_list <- list(
  destructive_lining_fibroblast_states = list(
    input_object = file.path(data_dir, "integrated_object", "destructive_lining_fibroblasts_subclustered.rds"),
    reduction_name = "umap.harmony.destructive.lining.fibroblast",
    cluster_col = "destructive_lining_fibroblast_subcluster",
    state_col = "destructive_lining_fibroblast_subtype",
    state_colors = destructive_lining_fibroblast_subtype_colors,
    cluster_to_state = destructive_lining_fibroblast_subcluster_labels,
    states_keep = unname(destructive_lining_fibroblast_subcluster_labels),
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
    input_object = file.path(data_dir, "integrated_object", "macrophages_subclustered.rds"),
    reduction_name = "umap.harmony.macrophage",
    cluster_col = "macrophage_subcluster",
    state_col = "macrophage_subtype",
    state_colors = macrophage_subtype_colors,
    cluster_to_state = macrophage_subcluster_labels,
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
    input_object = file.path(data_dir, "integrated_object", "activated_endothelial_subclustered.rds"),
    reduction_name = "umap.harmony.endothelial",
    cluster_col = "endothelial_subcluster",
    state_col = "endothelial_subtype",
    state_colors = endothelial_subtype_colors,
    cluster_to_state = endothelial_subcluster_labels,
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

# run each trajectory for the compartments
run_summary <- data.frame()

for (analysis_id in names(trajectory_list)) {
  cat("Running trajectory:", analysis_id, "\n")

  summary_row <- run_one_monocle_trajectory(
    analysis_id = analysis_id,
    params = trajectory_list[[analysis_id]],
    assay_use = assay_use,
    min_cells = min_cells,
    pt_size = pt_size,
    run_graph_test = run_graph_test,
    n_graph_test_genes = n_graph_test_genes,
    graph_test_cores = graph_test_cores,
    condition_colors = condition_colors,
    res_dir = res_dir,
    fig_dir = fig_dir
  )

  if (!is.null(summary_row)) {
    run_summary <- rbind(run_summary, summary_row)
  }
}

write.csv(run_summary, file.path(res_dir, "monocle3_trajectory_run_summary.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("Monocle3 trajectory analysis complete.\n")
cat("Results dir: ", res_dir, "\n", sep = "")
cat("Figures dir: ", fig_dir, "\n", sep = "")
cat("============================================================\n")
