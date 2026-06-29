######################################################
### destructive lining fibroblast helper functions ###
######################################################

# sort subcluster labels numerically when possible
sort_cluster_levels <- function(x) {
  x <- unique(as.character(x))
  suppressWarnings(x_num <- as.integer(x))
  if (all(!is.na(x_num))) return(as.character(sort(x_num)))
  sort(x)
}

# apply final destructive lining fibroblast subtype labels
add_destructive_lining_fibroblast_subtype_labels <- function(obj,
                                                             cluster_col = "destructive_lining_fibroblast_subcluster",
                                                             label_col = "destructive_lining_fibroblast_subtype") {
  cluster_ids <- as.character(obj@meta.data[[cluster_col]])
  labels <- unname(destructive_lining_fibroblast_subcluster_labels[cluster_ids])
  labels[is.na(labels)] <- "Unknown"
  obj@meta.data[[label_col]] <- factor(labels, levels = unname(destructive_lining_fibroblast_subcluster_labels))
  obj@meta.data[[label_col]] <- droplevels(obj@meta.data[[label_col]])
  obj
}

# get colors for the fibroblast subtypes present in the object
get_destructive_lining_fibroblast_colors <- function(labels) {
  labels <- unique(as.character(labels))
  cols <- destructive_lining_fibroblast_subtype_colors[labels]
  missing_labels <- labels[is.na(cols)]

  if (length(missing_labels) > 0) {
    fallback <- grDevices::hcl.colors(length(missing_labels), palette = "Set 3")
    names(fallback) <- missing_labels
    cols[names(fallback)] <- fallback
  }

  cols[labels]
}

# convert the fibroblast marker list into a simple subtype/gene table
make_destructive_lining_fibroblast_marker_table <- function(obj) {
  marker_rows <- list()

  for (subtype in names(marker_genes_destructive_lining_fibroblast)) {
    genes <- marker_genes_destructive_lining_fibroblast[[subtype]]
    genes <- genes[genes %in% rownames(obj)]
    if (length(genes) == 0) next
    marker_rows[[subtype]] <- data.frame(destructive_lining_fibroblast_subtype = subtype, gene = genes, stringsAsFactors = FALSE)
  }

  do.call(rbind, marker_rows)
}

# build cell ratios for subtype composition plots
build_fibroblast_ratio_plot_data <- function(meta_df, x_col, fill_col) {
  valid_rows <- !is.na(meta_df[[x_col]]) & trimws(as.character(meta_df[[x_col]])) != ""
  valid_rows <- valid_rows & !is.na(meta_df[[fill_col]]) & trimws(as.character(meta_df[[fill_col]])) != ""

  plot_df <- meta_df[valid_rows, c(x_col, fill_col), drop = FALSE]
  if (nrow(plot_df) == 0) return(NULL)

  plot_df[[x_col]] <- as.character(plot_df[[x_col]])
  plot_df[[fill_col]] <- as.character(plot_df[[fill_col]])

  count_df <- as.data.frame(table(plot_df[[x_col]], plot_df[[fill_col]]), stringsAsFactors = FALSE)
  colnames(count_df) <- c("group", "cell_type", "n_cells")
  count_df <- count_df[count_df$n_cells > 0, , drop = FALSE]
  if (nrow(count_df) == 0) return(NULL)

  totals <- aggregate(n_cells ~ group, data = count_df, FUN = sum)
  count_df <- merge(count_df, totals, by = "group", suffixes = c("", "_total"), sort = FALSE)
  count_df$ratio <- count_df$n_cells / count_df$n_cells_total
  count_df$group <- factor(count_df$group, levels = unique(count_df$group))
  count_df$cell_type <- factor(count_df$cell_type, levels = unique(as.character(plot_df[[fill_col]])))
  count_df
}

# build per-sample cell ratios and mark very small samples
build_fibroblast_sample_ratio_plot_data <- function(meta_df,
                                                    sample_col,
                                                    condition_col,
                                                    fill_col,
                                                    min_cells = 5,
                                                    low_cell_label = "Low cell count (<5 cells)") {
  valid_rows <- !is.na(meta_df[[sample_col]]) & trimws(as.character(meta_df[[sample_col]])) != ""
  valid_rows <- valid_rows & !is.na(meta_df[[condition_col]]) & trimws(as.character(meta_df[[condition_col]])) != ""
  valid_rows <- valid_rows & !is.na(meta_df[[fill_col]]) & trimws(as.character(meta_df[[fill_col]])) != ""

  plot_df <- meta_df[valid_rows, c(sample_col, condition_col, fill_col), drop = FALSE]
  if (nrow(plot_df) == 0) return(NULL)

  colnames(plot_df) <- c("sample_id", "condition", "cell_type")
  plot_df$sample_id <- as.character(plot_df$sample_id)
  plot_df$condition <- ifelse(as.character(plot_df$condition) == "HA", "HA", "other")
  plot_df$condition <- factor(plot_df$condition, levels = c("HA", "other"))
  plot_df$cell_type <- as.character(plot_df$cell_type)

  count_df <- as.data.frame(table(plot_df$sample_id, plot_df$condition, plot_df$cell_type), stringsAsFactors = FALSE)
  colnames(count_df) <- c("sample_id", "condition", "cell_type", "n_cells")
  count_df <- count_df[count_df$n_cells > 0, , drop = FALSE]
  if (nrow(count_df) == 0) return(NULL)

  totals <- aggregate(n_cells ~ sample_id, data = count_df, FUN = sum)
  colnames(totals)[2] <- "sample_total"
  count_df <- merge(count_df, totals, by = "sample_id", sort = FALSE)
  count_df$ratio <- count_df$n_cells / count_df$sample_total

  sample_info <- unique(plot_df[, c("sample_id", "condition"), drop = FALSE])
  sample_info <- sample_info[order(sample_info$condition, sample_info$sample_id), , drop = FALSE]

  low_count_samples <- totals$sample_id[totals$sample_total < min_cells]
  if (length(low_count_samples) > 0) {
    low_count_rows <- merge(
      data.frame(sample_id = low_count_samples, stringsAsFactors = FALSE),
      unique(count_df[, c("sample_id", "condition", "sample_total"), drop = FALSE]),
      by = "sample_id",
      sort = FALSE
    )
    low_count_rows$cell_type <- low_cell_label
    low_count_rows$n_cells <- low_count_rows$sample_total
    low_count_rows$ratio <- 1

    count_df <- count_df[!count_df$sample_id %in% low_count_samples, , drop = FALSE]
    count_df <- rbind(
      count_df[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total", "ratio"), drop = FALSE],
      low_count_rows[, c("sample_id", "condition", "cell_type", "n_cells", "sample_total", "ratio"), drop = FALSE]
    )
  }

  count_df$is_low_cell_count <- count_df$sample_id %in% low_count_samples
  count_df$sample_id <- factor(count_df$sample_id, levels = sample_info$sample_id)
  count_df$condition <- factor(as.character(count_df$condition), levels = c("HA", "other"))
  count_df$cell_type <- factor(count_df$cell_type, levels = c(levels(meta_df[[fill_col]]), low_cell_label))
  count_df
}

# plot a stacked ratio barplot
plot_fibroblast_ratio <- function(plot_df, fill_colors) {
  ggplot(plot_df, aes(x = group, y = ratio, fill = cell_type)) +
    geom_col(width = 0.92, color = NA) +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Ratio", fill = NULL) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.title.y = element_text(size = 18, color = "black"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length = grid::unit(0.22, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.4, "cm"),
      legend.box = "vertical",
      plot.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))
}

# plot stacked ratios for each sample
plot_fibroblast_sample_ratio <- function(plot_df, fill_colors, low_cell_label = "Low cell count (<5 cells)") {
  sample_fill_colors <- fill_colors
  sample_fill_colors[low_cell_label] <- "white"

  ggplot(plot_df, aes(x = sample_id, y = ratio, fill = cell_type)) +
    geom_col(width = 0.92, color = "black", linewidth = 0.25) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_manual(values = sample_fill_colors, breaks = setdiff(names(sample_fill_colors), low_cell_label), drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Ratio", fill = NULL) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.title.y = element_text(size = 18, color = "black"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length = grid::unit(0.22, "cm"),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 17, face = "bold", color = "black", margin = margin(t = 8)),
      panel.spacing.x = grid::unit(0.35, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.4, "cm"),
      legend.box = "vertical",
      plot.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 30, 20)
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))
}

# read one or more one-column gene-set files
load_gene_set_files <- function(files) {
  genes <- character(0)

  for (f in files) {
    if (!file.exists(f)) next
    x <- readxl::read_excel(f)
    if (ncol(x) == 0) next

    genes_here <- trimws(as.character(x[[1]]))
    genes_here <- toupper(genes_here[!is.na(genes_here) & genes_here != ""])
    genes <- c(genes, genes_here)
  }

  genes <- sort(unique(genes))
  genes <- genes[genes != "GENE"]
  genes[genes == "TRFC"] <- "TFRC"
  genes
}

# avoid infinite values in volcano plots
safe_neg_log10 <- function(x) {
  y <- -log10(x)
  if (any(is.infinite(y), na.rm = TRUE)) {
    max_finite <- max(y[is.finite(y)], na.rm = TRUE)
    if (!is.finite(max_finite)) max_finite <- 0
    y[is.infinite(y)] <- max_finite + 1
  }
  y
}
