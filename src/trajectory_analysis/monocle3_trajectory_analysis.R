###################################################
### Monocle3 trajectory and pseudotime analysis ###
###################################################

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
  DefaultAssay(obj) <- assay_use

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

  if (length(cells_keep) < min_cells) {
    warning("Too few cells for ", analysis_id, ": ", length(cells_keep))
    next
  }

  obj <- subset(obj, cells = cells_keep)

  # 5. Make sure the RNA data layer is usable.
  obj <- tryCatch(
    JoinLayers(obj, assay = assay_use),
    error = function(e) obj
  )

  rna_data <- tryCatch(
    LayerData(obj[[assay_use]], layer = "data"),
    error = function(e) NULL
  )

  if (is.null(rna_data) || length(rna_data) == 0) {
    obj <- NormalizeData(obj, assay = assay_use, verbose = FALSE)
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

  write.csv(cell_metadata, file.path(res_dir, paste0(analysis_id, "_cell_pseudotime_metadata.csv")), row.names = FALSE)

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

    write.csv(condition_summary, file.path(res_dir, paste0(analysis_id, "_pseudotime_summary_by_condition.csv")), row.names = FALSE)

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

      write.csv(wilcox_table, file.path(res_dir, paste0(analysis_id, "_pseudotime_HA_vs_other_wilcox.csv")), row.names = FALSE)
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

    write.csv(sample_summary, file.path(res_dir, paste0(analysis_id, "_pseudotime_summary_by_sample.csv")), row.names = FALSE)
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
    cell_size = pt_size
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
    cell_size = pt_size
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
      cell_size = pt_size
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
      cell_size = pt_size
    ) +
      ggtitle(paste0(analysis_id, ": sample_id")) +
      scale_color_manual(values = sample_color_map, drop = FALSE) +
      coord_equal() +
      theme_classic(base_size = 14)

    plot_list[[length(plot_list) + 1]] <- p_sample
  }

  p_review <- wrap_plots(plot_list, ncol = 2)

  ggsave(
    file.path(fig_dir, paste0(analysis_id, "_trajectory_review.png")),
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
      file.path(fig_dir, paste0(analysis_id, "_pseudotime_by_state_boxplot.png")),
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
  if (run_graph_test) {
    graph_genes <- VariableFeatures(obj)
    graph_genes <- intersect(graph_genes, rownames(cds))
    graph_genes <- head(graph_genes, n_graph_test_genes)
    graph_genes <- unique(c(marker_genes, graph_genes))

    if (length(graph_genes) > 0) {
      graph_res <- graph_test(
        cds[graph_genes, ],
        neighbor_graph = "principal_graph",
        cores = graph_test_cores
      )

      graph_res$gene <- rownames(graph_res)
      graph_res <- graph_res[order(graph_res$q_value, -graph_res$morans_I), ]

      write.csv(graph_res, file.path(res_dir, paste0(analysis_id, "_graph_test_pseudotime_genes.csv")), row.names = FALSE)
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

# save a small table summarizing which analyses were completed.
write.csv(run_summary, file.path(res_dir, "monocle3_trajectory_run_summary.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("Monocle3 trajectory analysis complete.\n")
cat("Results dir  : ", res_dir, "\n")
cat("Figures dir  : ", fig_dir, "\n")
cat("Objects dir  : ", out_dir, "\n")
cat("============================================================\n")
