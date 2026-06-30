###################################
### trajectory helper functions ###
###################################

# assign colors and fill missing labels with a fallback palette
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

# add final state labels to the subset object
add_trajectory_state_labels <- function(obj, cluster_col, state_col, cluster_to_state) {
  cluster_ids <- as.character(obj@meta.data[[cluster_col]])
  state_labels <- unname(cluster_to_state[cluster_ids])
  state_labels[is.na(state_labels)] <- "Unknown"
  obj@meta.data[[state_col]] <- state_labels
  obj
}

# prepare a Seurat subset object for Monocle using the existing UMAP
make_monocle_cds <- function(obj, params, assay_use) {
  DefaultAssay(obj) <- assay_use
  obj <- JoinLayers(obj, assay = assay_use)
  obj <- NormalizeData(obj, assay = assay_use, verbose = FALSE)

  cds <- as.cell_data_set(obj)
  rowData(cds)$gene_short_name <- rownames(cds)

  seurat_umap <- Embeddings(obj, reduction = params$reduction_name)
  seurat_umap <- seurat_umap[colnames(cds), , drop = FALSE]
  colnames(seurat_umap) <- c("UMAP_1", "UMAP_2")
  cds@int_colData@listData$reducedDims$UMAP <- seurat_umap

  monocle_partition <- factor(rep(1, ncol(cds)))
  names(monocle_partition) <- colnames(cds)
  cds@clusters$UMAP$partitions <- monocle_partition

  monocle_clusters <- factor(as.character(colData(cds)[[params$state_col]]))
  names(monocle_clusters) <- colnames(cds)
  cds@clusters$UMAP$clusters <- monocle_clusters

  list(obj = obj, cds = cds)
}

# learn graph, order cells and return pseudotime metadata
order_monocle_cells <- function(obj, cds, params) {
  cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)

  root_cells <- rownames(obj@meta.data)[
    as.character(obj@meta.data[[params$state_col]]) %in% params$root_states
  ]
  root_cells <- intersect(root_cells, colnames(cds))

  if (length(root_cells) == 0) {
    stop("No root cells found for: ", params$analysis_id)
  }

  cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

  pseudotime_values <- pseudotime(cds)
  pseudotime_values[!is.finite(pseudotime_values)] <- NA_real_
  colData(cds)$pseudotime <- pseudotime_values[colnames(cds)]

  cell_metadata <- as.data.frame(colData(cds))
  cell_metadata$cell_id <- rownames(cell_metadata)
  cell_metadata$pseudotime <- as.numeric(cell_metadata$pseudotime)
  cell_metadata$is_root_cell <- cell_metadata$cell_id %in% root_cells

  list(cds = cds, cell_metadata = cell_metadata, root_cells = root_cells)
}

# save pseudotime summaries by condition and sample
save_pseudotime_tables <- function(cell_metadata, analysis_id, res_dir) {
  write.csv(cell_metadata, file.path(res_dir, paste0(analysis_id, "_cell_pseudotime_metadata.csv")), row.names = FALSE)

  condition_summary <- aggregate(
    pseudotime ~ condition,
    data = cell_metadata,
    FUN = function(x) c(n_cells = length(x), median = median(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE))
  )
  condition_summary <- do.call(data.frame, condition_summary)
  write.csv(condition_summary, file.path(res_dir, paste0(analysis_id, "_pseudotime_summary_by_condition.csv")), row.names = FALSE)

  test_df <- cell_metadata[cell_metadata$condition %in% c("HA", "other") & !is.na(cell_metadata$pseudotime), , drop = FALSE]

  if (all(c("HA", "other") %in% unique(as.character(test_df$condition)))) {
    wilcox_res <- wilcox.test(pseudotime ~ condition, data = test_df)
    wilcox_table <- data.frame(
      analysis_id = analysis_id,
      n_HA = sum(test_df$condition == "HA"),
      n_other = sum(test_df$condition == "other"),
      median_HA = median(test_df$pseudotime[test_df$condition == "HA"]),
      median_other = median(test_df$pseudotime[test_df$condition == "other"]),
      wilcox_pvalue = wilcox_res$p.value,
      stringsAsFactors = FALSE
    )
    write.csv(wilcox_table, file.path(res_dir, paste0(analysis_id, "_pseudotime_HA_vs_other_wilcox.csv")), row.names = FALSE)
  }

  sample_summary <- aggregate(
    pseudotime ~ sample_id + condition,
    data = cell_metadata,
    FUN = function(x) c(n_cells = length(x), median = median(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE))
  )
  sample_summary <- do.call(data.frame, sample_summary)
  write.csv(sample_summary, file.path(res_dir, paste0(analysis_id, "_pseudotime_summary_by_sample.csv")), row.names = FALSE)
}

# draw one Monocle trajectory plot
make_trajectory_plot <- function(cds, color_by, title, pt_size, color_map = NULL) {
  p <- plot_cells(
    cds,
    color_cells_by = color_by,
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_roots = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 1.5,
    cell_size = pt_size
  ) +
    ggtitle(title) +
    coord_equal() +
    theme_classic(base_size = 14)

  if (!is.null(color_map)) {
    p <- p + scale_color_manual(values = color_map, drop = FALSE)
  }

  p
}

# save the four-panel trajectory review plot
save_trajectory_review_plot <- function(cds, cell_metadata, params, analysis_id, fig_dir, pt_size, condition_colors) {
  state_color_map <- build_label_colors(unique(cell_metadata[[params$state_col]]), params$state_colors)
  condition_color_map <- build_label_colors(unique(cell_metadata$condition), condition_colors)

  sample_values <- sort(unique(as.character(cell_metadata$sample_id)))
  sample_values <- sample_values[!is.na(sample_values)]
  sample_color_map <- grDevices::hcl.colors(length(sample_values), palette = "Dark 3")
  names(sample_color_map) <- sample_values

  p_review <- wrap_plots(
    list(
      make_trajectory_plot(cds, "pseudotime", paste0(analysis_id, ": pseudotime"), pt_size),
      make_trajectory_plot(cds, params$state_col, paste0(analysis_id, ": states"), pt_size, state_color_map),
      make_trajectory_plot(cds, "condition", paste0(analysis_id, ": condition"), pt_size, condition_color_map),
      make_trajectory_plot(cds, "sample_id", paste0(analysis_id, ": sample_id"), pt_size, sample_color_map)
    ),
    ncol = 2
  )

  ggsave(file.path(fig_dir, paste0(analysis_id, "_trajectory_review.png")), p_review, width = 16, height = 12, dpi = 600)
}

# save the pseudotime distribution across final states
save_pseudotime_state_boxplot <- function(cell_metadata, params, analysis_id, fig_dir) {
  boxplot_df <- cell_metadata[!is.na(cell_metadata$pseudotime), ]
  boxplot_df$state_for_plot <- boxplot_df[[params$state_col]]
  boxplot_df <- boxplot_df[!is.na(boxplot_df$state_for_plot), ]
  state_color_map <- build_label_colors(unique(boxplot_df$state_for_plot), params$state_colors)

  p_state_boxplot <- ggplot(
    boxplot_df,
    aes(x = pseudotime, y = reorder(state_for_plot, pseudotime, median), fill = state_for_plot)
  ) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    scale_fill_manual(values = state_color_map, drop = FALSE) +
    labs(title = paste0(analysis_id, ": pseudotime by state"), x = "pseudotime", y = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  n_states <- length(unique(boxplot_df$state_for_plot))
  ggsave(file.path(fig_dir, paste0(analysis_id, "_pseudotime_by_state_boxplot.png")), p_state_boxplot, width = 10, height = max(4, n_states * 0.6), dpi = 600)
}

# optionally test genes associated with the principal graph
run_monocle_graph_test <- function(cds, obj, params, analysis_id, res_dir, n_graph_test_genes, graph_test_cores) {
  marker_genes <- params$marker_genes[params$marker_genes %in% rownames(cds)]
  graph_genes <- VariableFeatures(obj)
  graph_genes <- intersect(graph_genes, rownames(cds))
  graph_genes <- head(graph_genes, n_graph_test_genes)
  graph_genes <- unique(c(marker_genes, graph_genes))

  if (length(graph_genes) == 0) return(invisible(NULL))

  graph_res <- graph_test(cds[graph_genes, ], neighbor_graph = "principal_graph", cores = graph_test_cores)
  graph_res$gene <- rownames(graph_res)
  graph_res <- graph_res[order(graph_res$q_value, -graph_res$morans_I), ]
  write.csv(graph_res, file.path(res_dir, paste0(analysis_id, "_graph_test_pseudotime_genes.csv")), row.names = FALSE)

  invisible(graph_res)
}

# run one complete trajectory analysis
run_one_monocle_trajectory <- function(analysis_id,
                                       params,
                                       assay_use,
                                       min_cells,
                                       pt_size,
                                       run_graph_test,
                                       n_graph_test_genes,
                                       graph_test_cores,
                                       condition_colors,
                                       res_dir,
                                       fig_dir) {
  params$analysis_id <- analysis_id
  obj <- readRDS(params$input_object)
  obj <- add_trajectory_state_labels(obj, params$cluster_col, params$state_col, params$cluster_to_state)

  cells_keep <- rownames(obj@meta.data)[as.character(obj@meta.data[[params$state_col]]) %in% params$states_keep]
  if (length(cells_keep) < min_cells) {
    warning("Too few cells for ", analysis_id, ": ", length(cells_keep))
    return(NULL)
  }

  obj <- subset(obj, cells = cells_keep)

  monocle_input <- make_monocle_cds(obj, params, assay_use)
  obj <- monocle_input$obj
  cds <- monocle_input$cds

  ordered <- order_monocle_cells(obj, cds, params)
  cds <- ordered$cds
  cell_metadata <- ordered$cell_metadata

  save_pseudotime_tables(cell_metadata, analysis_id, res_dir)
  save_trajectory_review_plot(cds, cell_metadata, params, analysis_id, fig_dir, pt_size, condition_colors)
  save_pseudotime_state_boxplot(cell_metadata, params, analysis_id, fig_dir)

  if (run_graph_test) {
    run_monocle_graph_test(cds, obj, params, analysis_id, res_dir, n_graph_test_genes, graph_test_cores)
  }

  data.frame(
    analysis_id = analysis_id,
    input_object = params$input_object,
    reduction_name = params$reduction_name,
    cluster_col = params$cluster_col,
    state_col = params$state_col,
    n_cells = ncol(obj),
    n_root_cells = length(ordered$root_cells),
    n_genes = nrow(cds),
    stringsAsFactors = FALSE
  )
}
