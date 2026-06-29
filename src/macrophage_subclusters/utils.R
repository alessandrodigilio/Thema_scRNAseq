##############################################
### macrophage subcluster helper functions ###
##############################################

# apply final macrophage subtype labels to a subclustered object
add_macrophage_subtype_labels <- function(obj, cluster_col = "macrophage_subcluster", label_col = "macrophage_subtype") {
  cluster_ids <- as.character(obj@meta.data[[cluster_col]])
  labels <- unname(macrophage_subcluster_labels[cluster_ids])
  labels[is.na(labels)] <- "Unknown"
  obj@meta.data[[label_col]] <- factor(labels, levels = unname(macrophage_subcluster_labels))
  obj@meta.data[[label_col]] <- droplevels(obj@meta.data[[label_col]])
  obj
}

# get colors for the macrophage subtypes present in the object
get_macrophage_colors <- function(labels) {
  labels <- unique(as.character(labels))
  macrophage_subtype_colors[labels]
}

# summarize top metadata labels inside each subcluster for manual review
summarize_label_composition <- function(meta_df, cluster_col, label_col, prefix, n_top = 3) {
  if (!label_col %in% colnames(meta_df)) return(NULL)

  df <- meta_df %>%
    filter(!is.na(.data[[label_col]]) & trimws(.data[[label_col]]) != "") %>%
    group_by(.data[[cluster_col]], .data[[label_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(freq = n / sum(n)) %>%
    arrange(.data[[cluster_col]], desc(freq), desc(n)) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= n_top) %>%
    ungroup()

  if (nrow(df) == 0) return(NULL)

  label_wide <- data.frame(cluster_id = as.character(df[[cluster_col]]), rank = df$rank, value = as.character(df[[label_col]]), stringsAsFactors = FALSE) %>%
    pivot_wider(names_from = rank, values_from = value, names_glue = paste0(prefix, "_top{rank}_label"))

  freq_wide <- data.frame(cluster_id = as.character(df[[cluster_col]]), rank = df$rank, value = df$freq, stringsAsFactors = FALSE) %>%
    pivot_wider(names_from = rank, values_from = value, names_glue = paste0(prefix, "_top{rank}_freq"))

  left_join(label_wide, freq_wide, by = "cluster_id")
}

# make the marker table used by the macrophage annotation dotplot
make_macrophage_marker_table <- function(obj) {
  marker_rows <- list()

  for (subtype in names(marker_genes_macrophage)) {
    genes <- marker_genes_macrophage[[subtype]]
    genes <- genes[genes %in% rownames(obj)]
    if (length(genes) == 0) next
    marker_rows[[subtype]] <- data.frame(macrophage_subtype = subtype, gene = genes, stringsAsFactors = FALSE)
  }

  do.call(rbind, marker_rows)
}

# compute average expression and percent expression for a marker dotplot
build_macrophage_marker_dotplot_data <- function(object, features, group_col) {
  plot_df <- FetchData(object, vars = c(group_col, features))
  colnames(plot_df)[1] <- "group"
  plot_df$group <- as.character(plot_df$group)
  group_levels <- levels(object@meta.data[[group_col]])

  res_list <- list()

  for (gene in features) {
    avg_expr <- tapply(plot_df[[gene]], plot_df$group, mean, na.rm = TRUE)
    pct_expr <- tapply(plot_df[[gene]] > 0, plot_df$group, mean, na.rm = TRUE) * 100
    avg_expr <- avg_expr[group_levels]
    pct_expr <- pct_expr[group_levels]
    avg_expr[is.na(avg_expr)] <- 0
    pct_expr[is.na(pct_expr)] <- 0

    if (stats::sd(avg_expr) > 0) {
      avg_scaled <- as.numeric(scale(avg_expr))
    } else {
      avg_scaled <- rep(0, length(avg_expr))
    }

    res_list[[gene]] <- data.frame(group = group_levels, gene = gene, avg_scaled = avg_scaled, pct_expr = as.numeric(pct_expr), stringsAsFactors = FALSE)
  }

  res_df <- do.call(rbind, res_list)
  res_df <- res_df[res_df$pct_expr > 0, , drop = FALSE]
  res_df$group <- factor(res_df$group, levels = group_levels)
  res_df$gene <- factor(res_df$gene, levels = features)
  res_df
}

# plot a macrophage marker bubble plot
plot_macrophage_marker_dotplot <- function(dotplot_df, low_col, mid_col, high_col) {
  ggplot(dotplot_df, aes(x = gene, y = group, size = pct_expr, color = avg_scaled)) +
    geom_point() +
    scale_size_continuous(name = "Percent Expressed", range = c(0.6, 6)) +
    scale_color_gradient2(name = "Average Expression", low = low_col, mid = mid_col, high = high_col, midpoint = 0) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11),
      legend.position = "top",
      legend.key.width = grid::unit(0.45, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      plot.margin = margin(20, 35, 35, 25)
    ) +
    guides(
      color = guide_colorbar(title.position = "top", barwidth = grid::unit(2.4, "cm"), barheight = grid::unit(0.35, "cm")),
      size = guide_legend(title.position = "top")
    )
}

# plot sample-level mean expression for one gene in one macrophage subtype
plot_gene_by_sample_in_macrophage_subtype <- function(obj,
                                                      gene,
                                                      macrophage_subtype,
                                                      group_col,
                                                      sample_col,
                                                      subtype_col,
                                                      group_levels,
                                                      group_colors) {
  cells_use <- rownames(obj@meta.data)[obj@meta.data[[subtype_col]] == macrophage_subtype]
  obj_sub <- subset(obj, cells = cells_use)

  plot_df <- FetchData(obj_sub, vars = c(gene, group_col, sample_col))
  colnames(plot_df) <- c("expression", "condition", "sample_id")
  plot_df$condition <- ifelse(as.character(plot_df$condition) == "HA", "HA", "other")
  plot_df <- plot_df[!is.na(plot_df$condition) & !is.na(plot_df$sample_id) & as.character(plot_df$sample_id) != "", , drop = FALSE]

  sample_df <- aggregate(expression ~ sample_id + condition, data = plot_df, FUN = mean)
  sample_df$condition <- factor(sample_df$condition, levels = group_levels)

  p <- ggplot(sample_df, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot(width = 0.44, alpha = 0.9, color = "black", linewidth = 0.45, outlier.shape = NA) +
    scale_fill_manual(values = group_colors) +
    labs(title = paste0(gene, " in ", macrophage_subtype), x = NULL, y = "Mean norm. expression per sample") +
    theme_classic(base_size = 20) +
    theme(
      axis.text = element_text(size = 20, color = "black"),
      axis.title = element_text(size = 20, color = "black"),
      plot.title = element_text(size = 17, hjust = 0.5, color = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 1.1, color = "black"),
      axis.ticks = element_line(linewidth = 1.1, color = "black"),
      axis.ticks.length = grid::unit(0.2, "cm"),
      plot.margin = margin(18, 18, 18, 18)
    )

  list(plot = p, data = sample_df)
}

# draw average ferroptosis gene expression by macrophage subtype and condition
make_macrophage_ferroptosis_heatmap <- function(obj, genes_use, label_col) {
  dot_features <- genes_use[genes_use %in% rownames(obj)]
  subtype_levels <- levels(obj@meta.data[[label_col]])

  obj$subtype_group <- paste(obj@meta.data[[label_col]], obj$HA_vs_other, sep = " | ")
  group_levels <- as.vector(rbind(
    paste(subtype_levels, "HA", sep = " | "),
    paste(subtype_levels, "other", sep = " | ")
  ))

  obj$subtype_group <- factor(obj$subtype_group, levels = group_levels)

  avg_expr <- AverageExpression(
    object = obj,
    assays = "RNA",
    features = dot_features,
    group.by = "subtype_group",
    slot = "data",
    verbose = FALSE
  )$RNA

  avg_expr <- avg_expr[, group_levels[group_levels %in% colnames(avg_expr)], drop = FALSE]
  scaled_expr <- t(scale(t(as.matrix(avg_expr))))
  scaled_expr[is.na(scaled_expr)] <- 0

  heatmap_df <- as.data.frame(as.table(scaled_expr), stringsAsFactors = FALSE)
  colnames(heatmap_df) <- c("gene", "group", "z_score")
  heatmap_df$gene <- factor(heatmap_df$gene, levels = rev(dot_features))
  heatmap_df$group <- factor(heatmap_df$group, levels = colnames(scaled_expr))

  heatmap_groups <- strsplit(as.character(heatmap_df$group), " | ", fixed = TRUE)
  heatmap_df$macrophage_subtype <- factor(vapply(heatmap_groups, `[`, character(1), 1), levels = subtype_levels)
  heatmap_df$HA_vs_other <- factor(vapply(heatmap_groups, `[`, character(1), 2), levels = c("HA", "other"))

  p <- ggplot(heatmap_df, aes(x = HA_vs_other, y = gene, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(cols = vars(macrophage_subtype), scales = "free_x", space = "free_x", labeller = labeller(macrophage_subtype = wrap_each_word)) +
    scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#7A1F2B", midpoint = 0, name = "Scaled\nexpression") +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.spacing.x = grid::unit(0.22, "cm"),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 8.5, color = "black", lineheight = 0.95),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11),
      legend.position = "top",
      legend.key.width = grid::unit(0.5, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.08, "cm"),
      plot.margin = margin(20, 35, 35, 25)
    )

  list(plot = p, n_subtypes = length(unique(heatmap_df$macrophage_subtype)), n_genes = nrow(scaled_expr))
}
