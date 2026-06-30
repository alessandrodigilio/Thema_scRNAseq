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

# -------------------------------------------------------------------------
# functional analysis helpers used by 03_functional_analysis_ha_enriched_states.R
# -------------------------------------------------------------------------

# make file-safe pathway names
sanitize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# load one MSigDB collection with support for old and new msigdbr arguments
load_msigdb_table <- function(species, collection, subcollection) {
  out <- tryCatch(
    msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) > 0) return(as.data.frame(out))

  out <- tryCatch(
    msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) > 0) return(as.data.frame(out))

  NULL
}

# load GO BP, KEGG and Reactome pathways for fgsea
load_fibroblast_msigdb_pathways <- function(species = "Homo sapiens") {
  pathway_tables <- list()

  gobp_df <- load_msigdb_table(species = species, collection = "C5", subcollection = "GO:BP")
  if (!is.null(gobp_df) && nrow(gobp_df) > 0) {
    gobp_df$database <- "GO_BP"
    pathway_tables[["GO_BP"]] <- gobp_df
  }

  kegg_tables <- list()
  for (subcollection in c("CP:KEGG", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS")) {
    df <- load_msigdb_table(species = species, collection = "C2", subcollection = subcollection)
    if (!is.null(df) && nrow(df) > 0) {
      df$database <- "KEGG"
      kegg_tables[[subcollection]] <- df
    }
  }
  if (length(kegg_tables) > 0) pathway_tables[["KEGG"]] <- do.call(rbind, kegg_tables)

  reactome_df <- load_msigdb_table(species = species, collection = "C2", subcollection = "CP:REACTOME")
  if (!is.null(reactome_df) && nrow(reactome_df) > 0) {
    reactome_df$database <- "Reactome"
    pathway_tables[["Reactome"]] <- reactome_df
  }

  pathway_df <- do.call(rbind, pathway_tables)
  pathway_df$gene_symbol <- toupper(as.character(pathway_df$gene_symbol))
  pathway_df$pathway_id <- paste(pathway_df$database, pathway_df$gs_name, sep = "__")
  pathway_df <- pathway_df[
    !is.na(pathway_df$gene_symbol) & pathway_df$gene_symbol != "" &
      !is.na(pathway_df$gs_name) & pathway_df$gs_name != "",
    ,
    drop = FALSE
  ]

  pathways <- split(pathway_df$gene_symbol, pathway_df$pathway_id)
  pathways <- lapply(pathways, function(x) sort(unique(x)))
  pathways <- pathways[vapply(pathways, length, integer(1)) > 0]

  pathway_info <- unique(pathway_df[, c("pathway_id", "database", "gs_name"), drop = FALSE])
  pathway_info <- as.data.frame(pathway_info, stringsAsFactors = FALSE)
  rownames(pathway_info) <- pathway_info$pathway_id

  list(pathways = pathways, pathway_info = pathway_info)
}

# flag pathways by name using selected biological patterns
pathway_matches_patterns <- function(pathway_name, patterns) {
  grepl(paste(patterns, collapse = "|"), tolower(pathway_name), perl = TRUE)
}

# build ranked vector for fgsea from FindMarkers output
build_findmarkers_rank_vector <- function(de_df, rank_method = "avg_log2FC") {
  if (rank_method == "avg_log2FC") {
    rank_df <- de_df[!is.na(de_df$avg_log2FC), c("gene", "avg_log2FC"), drop = FALSE]
    colnames(rank_df) <- c("gene", "rank_value")
  } else {
    rank_df <- de_df[
      !is.na(de_df$avg_log2FC) & !is.na(de_df$p_val),
      c("gene", "avg_log2FC", "p_val"),
      drop = FALSE
    ]

    positive_p <- rank_df$p_val[rank_df$p_val > 0]
    min_positive_p <- min(positive_p, na.rm = TRUE)
    if (!is.finite(min_positive_p)) min_positive_p <- .Machine$double.xmin

    rank_df$p_val_rank <- rank_df$p_val
    rank_df$p_val_rank[rank_df$p_val_rank <= 0] <- min_positive_p * 0.1
    rank_df$rank_value <- sign(rank_df$avg_log2FC) * -log10(rank_df$p_val_rank)
    rank_df <- rank_df[, c("gene", "rank_value"), drop = FALSE]
  }

  rank_df$gene <- toupper(as.character(rank_df$gene))
  rank_df <- rank_df[!is.na(rank_df$gene) & rank_df$gene != "" & !is.na(rank_df$rank_value), , drop = FALSE]
  rank_df <- rank_df[rank_df$rank_value != 0, , drop = FALSE]
  rank_df <- rank_df[order(abs(rank_df$rank_value), decreasing = TRUE), , drop = FALSE]
  rank_df <- rank_df[!duplicated(rank_df$gene), , drop = FALSE]

  ranks <- rank_df$rank_value
  names(ranks) <- rank_df$gene
  sort(ranks, decreasing = TRUE)
}

# add readable pathway metadata to fgsea output
add_fibroblast_pathway_metadata <- function(fgsea_res, pathway_info, iron_patterns) {
  fgsea_df <- as.data.frame(fgsea_res)
  if (nrow(fgsea_df) == 0) return(fgsea_df)

  fgsea_df$pathway_id <- as.character(fgsea_df$pathway)
  fgsea_df$database <- pathway_info[fgsea_df$pathway_id, "database"]
  fgsea_df$pathway_name <- pathway_info[fgsea_df$pathway_id, "gs_name"]
  fgsea_df$is_iron_related <- pathway_matches_patterns(fgsea_df$pathway_name, iron_patterns)
  fgsea_df$leadingEdge <- vapply(fgsea_df$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))

  fgsea_df <- fgsea_df[
    order(fgsea_df$padj, -abs(fgsea_df$NES)),
    c("database", "pathway_name", "pathway_id", "pval", "padj", "ES", "NES", "size", "leadingEdge", "is_iron_related"),
    drop = FALSE
  ]

  fgsea_df
}

# save GSEA curves for selected significant pathways
save_fibroblast_gsea_curves <- function(curve_rows, ranks, pathways, fgsea_curve_dir) {
  for (i in seq_len(nrow(curve_rows))) {
    pathway_id <- curve_rows$pathway_id[i]
    pathway_name <- curve_rows$pathway_name[i]
    direction <- ifelse(curve_rows$NES[i] > 0, "HA_enriched_states", "other_states")

    p <- plotEnrichment(pathways[[pathway_id]], ranks) +
      labs(title = pathway_name, x = "Ranked genes", y = "Enrichment score") +
      theme_classic(base_size = 16) +
      theme(
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
        panel.grid = element_blank(),
        plot.margin = margin(18, 18, 18, 18)
      )

    output_file <- file.path(fgsea_curve_dir, paste0("gsea_curve_", direction, "_", sanitize_label(pathway_name), ".png"))
    ggsave(output_file, plot = p, width = 8.5, height = 5.5, dpi = 600)
  }
}

# run fgsea and save the main result tables and selected curves
run_fibroblast_fgsea_analysis <- function(de_df,
                                          func_res_dir,
                                          fgsea_table_dir,
                                          fgsea_iron_table_dir,
                                          fgsea_curve_dir,
                                          selected_gsea_curves,
                                          iron_patterns,
                                          rank_method = "avg_log2FC",
                                          padj_thr = 0.05,
                                          min_size = 5,
                                          max_size = 500,
                                          fgsea_seed = 1234,
                                          use_multilevel = TRUE,
                                          n_permutations = 10000,
                                          msigdb_species = "Homo sapiens") {
  cat("Ranking method:", rank_method, "\n")
  ranks <- build_findmarkers_rank_vector(de_df, rank_method = rank_method)
  cat("Ranked genes:", length(ranks), "\n")
  cat("Top ranked genes toward HA-enriched DLF states:\n")
  print(head(ranks, 10))
  cat("Bottom ranked genes toward other DLF states:\n")
  print(tail(ranks, 10))

  pathway_data <- load_fibroblast_msigdb_pathways(species = msigdb_species)
  pathways <- pathway_data$pathways
  pathway_info <- pathway_data$pathway_info

  cat("Pathways loaded:\n")
  print(table(pathway_info$database))

  set.seed(fgsea_seed)
  if (use_multilevel) {
    fgsea_res <- fgseaMultilevel(pathways = pathways, stats = ranks, minSize = min_size, maxSize = max_size)
  } else {
    fgsea_res <- fgsea(pathways = pathways, stats = ranks, minSize = min_size, maxSize = max_size, nperm = n_permutations)
  }

  fgsea_df <- add_fibroblast_pathway_metadata(fgsea_res, pathway_info, iron_patterns)
  sig_df <- fgsea_df[!is.na(fgsea_df$padj) & fgsea_df$padj < padj_thr, , drop = FALSE]
  iron_sig_df <- sig_df[sig_df$is_iron_related, , drop = FALSE]

  all_file <- file.path(func_res_dir, "fgsea_all_results_HA_enriched_DLF_states_vs_other_DLF_states.csv")
  sig_file <- file.path(fgsea_table_dir, "fgsea_significant_HA_enriched_DLF_states_vs_other_DLF_states.csv")
  iron_file <- file.path(fgsea_iron_table_dir, "fgsea_iron_related_significant_HA_enriched_DLF_states_vs_other_DLF_states.csv")

  write.csv(fgsea_df, all_file, row.names = FALSE)
  write.csv(sig_df, sig_file, row.names = FALSE)
  write.csv(iron_sig_df, iron_file, row.names = FALSE)

  cat("Significant pathways at FDR < ", padj_thr, ": ", nrow(sig_df), "\n", sep = "")
  cat("Significant iron/heme/ferroptosis/oxidative-stress pathways: ", nrow(iron_sig_df), "\n", sep = "")
  cat("\nTop significant pathways:\n")
  print(head(sig_df[, c("database", "pathway_name", "padj", "NES", "size", "is_iron_related"), drop = FALSE], 20))

  if (nrow(iron_sig_df) > 0) {
    cat("\nIron-related significant pathways:\n")
    print(iron_sig_df[, c("database", "pathway_name", "padj", "NES", "size"), drop = FALSE])
  }

  curve_rows <- sig_df[sig_df$pathway_name %in% selected_gsea_curves, , drop = FALSE]
  missing_curves <- setdiff(selected_gsea_curves, curve_rows$pathway_name)

  if (length(missing_curves) > 0) {
    cat("Selected GSEA curves not significant/found and therefore not plotted:\n")
    print(missing_curves)
  }

  curve_rows <- curve_rows[!duplicated(curve_rows$pathway_id), , drop = FALSE]
  write.csv(curve_rows, file.path(func_res_dir, "fgsea_pathways_selected_for_curve_plotting.csv"), row.names = FALSE)

  if (nrow(curve_rows) > 0) {
    save_fibroblast_gsea_curves(curve_rows, ranks, pathways, fgsea_curve_dir)
  }

  list(
    ranks = ranks,
    sig_df = sig_df,
    iron_sig_df = iron_sig_df,
    all_file = all_file,
    sig_file = sig_file,
    iron_file = iron_file
  )
}

# clean Enrichr pathway names for plotting
clean_enrichr_terms <- function(df) {
  df$Term <- gsub("\\s+R-HSA-[0-9]+$", "", df$Term)
  df$Term <- gsub("\\s+KEGG_[0-9]+$", "", df$Term)
  df$Term <- gsub("\\s*\\([^\\)]+\\)", "", df$Term)
  df$Term <- trimws(df$Term)
  df
}

# wrap long Enrichr labels over multiple lines
format_selected_term_label <- function(x) {
  x <- gsub(" \\(", "\n(", x, fixed = TRUE)
  vapply(strsplit(x, "\n", fixed = TRUE), function(parts) {
    first_line_words <- strsplit(parts[1], " ", fixed = TRUE)[[1]]
    wrapped_first <- paste(tapply(first_line_words, ceiling(seq_along(first_line_words) / 5), paste, collapse = " "), collapse = "\n")
    if (length(parts) > 1) paste(c(wrapped_first, parts[-1]), collapse = "\n") else wrapped_first
  }, character(1))
}

# prepare iron-related Enrichr terms for plotting
prepare_enrichr_iron_plot_df <- function(input_df, input_gene_count, db_label, iron_patterns, padj_thr = 0.05, top_n = 10) {
  df <- clean_enrichr_terms(as.data.frame(input_df))
  if (nrow(df) == 0) return(NULL)

  df$num_genes <- sapply(strsplit(df$Genes, ";"), length)
  df$gene_ratio <- df$num_genes / max(1, input_gene_count)
  df$is_significant <- !is.na(df$Adjusted.P.value) & df$Adjusted.P.value < padj_thr
  df <- df[pathway_matches_patterns(df$Term, iron_patterns) & df$is_significant, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  df <- df[order(df$gene_ratio, df$Adjusted.P.value), , drop = FALSE]
  df <- head(df, top_n)
  df$Database <- db_label
  df$Term <- paste0(df$Term, " (", df$Database, ")")
  df$Term <- format_selected_term_label(df$Term)
  df$Term <- factor(df$Term, levels = unique(df$Term))
  df
}

# dotplot for selected Enrichr pathways
build_enrichr_dotplot <- function(plot_df, color_palette = c("#8B3E2F", "coral")) {
  size_breaks <- unique(round(pretty(plot_df$num_genes)))
  size_breaks <- size_breaks[size_breaks >= min(plot_df$num_genes) & size_breaks <= max(plot_df$num_genes)]
  color_breaks <- unique(signif(seq(min(plot_df$Adjusted.P.value), max(plot_df$Adjusted.P.value), length.out = 3), 3))

  ggplot(plot_df, aes(x = gene_ratio, y = Term)) +
    geom_point(aes(size = num_genes, color = Adjusted.P.value), alpha = 0.9) +
    scale_color_continuous(low = color_palette[1], high = color_palette[2], name = "adj p-value", breaks = color_breaks) +
    scale_size_continuous(name = "Genes", range = c(3, 9), breaks = size_breaks) +
    labs(title = NULL, x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    guides(
      size = guide_legend(title.position = "top", nrow = 1, byrow = TRUE, order = 1, keywidth = grid::unit(0.18, "cm"), keyheight = grid::unit(0.32, "cm")),
      color = guide_colorbar(title.position = "top", barwidth = grid::unit(2.8, "cm"), barheight = grid::unit(0.35, "cm"), label.position = "bottom", label = scales::label_number(accuracy = 0.01), title.theme = element_text(margin = margin(b = 4)), label.theme = element_text(size = 8, margin = margin(t = 4)), order = 2)
    )
}

# tile plot showing pathway genes and their fold changes
build_selected_pathway_gene_tile_plot <- function(pathway_df, deg_df) {
  sig_df <- pathway_df[pathway_df$is_significant, , drop = FALSE]
  if (nrow(sig_df) == 0) return(NULL)

  gene_tiles <- do.call(rbind, lapply(seq_len(nrow(sig_df)), function(i) {
    genes <- unique(strsplit(sig_df$Genes[i], ";", fixed = TRUE)[[1]])
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    data.frame(Term = rep(as.character(sig_df$Term[i]), length(genes)), gene = genes, gene_index = seq_along(genes), stringsAsFactors = FALSE)
  }))

  gene_tiles <- merge(gene_tiles, unique(deg_df[, c("gene", "log2FoldChange")]), by = "gene", all.x = TRUE, sort = FALSE)
  gene_tiles$Term <- factor(gene_tiles$Term, levels = levels(pathway_df$Term))
  gene_tiles$segment_width <- 1.2
  gene_tiles <- gene_tiles[order(gene_tiles$Term, gene_tiles$gene_index), , drop = FALSE]
  gene_tiles$x_center <- ave(gene_tiles$segment_width, gene_tiles$Term, FUN = function(x) cumsum(x) - (x / 2))

  max_genes <- max(as.numeric(table(gene_tiles$Term)))
  max_abs_lfc <- max(abs(gene_tiles$log2FoldChange), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) max_abs_lfc <- 1

  ggplot(gene_tiles, aes(y = Term, x = segment_width, fill = log2FoldChange)) +
    geom_col(width = 0.78, position = "stack", color = "black", linewidth = 0.35) +
    geom_text(data = gene_tiles, aes(x = x_center, y = Term, label = gene), size = 3.2, color = "black", inherit.aes = FALSE) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, limits = c(-max_abs_lfc, max_abs_lfc), oob = scales::squish, name = "log2FC") +
    scale_x_continuous(limits = c(0, max_genes), breaks = 0:max_genes, labels = rep("", max_genes + 1), expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 14, color = "black"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank(),
      legend.position = "top"
    ) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = grid::unit(2.6, "cm"), barheight = grid::unit(0.35, "cm"), label.theme = element_text(size = 8, margin = margin(t = 4))))
}

# retry Enrichr calls because it is an online service
run_enrichr_with_retry <- function(genes, databases, max_tries = 4, retry_wait_sec = 20) {
  old_timeout <- getOption("timeout")
  options(timeout = max(300, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  last_error <- NULL

  for (attempt in seq_len(max_tries)) {
    cat("Enrichr attempt ", attempt, "/", max_tries, " with ", length(genes), " genes...\n", sep = "")

    res <- tryCatch(
      enrichr(genes, databases),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(res) && length(res) > 0) return(res)

    cat("Enrichr attempt failed: ", last_error, "\n", sep = "")
    if (attempt < max_tries) Sys.sleep(retry_wait_sec)
  }

  cat("Enrichr failed after all retries - skipping this direction.\n")
  NULL
}

# run Enrichr ORA on significant HA-enriched and other-enriched genes
run_fibroblast_enrichr_ora <- function(de_df,
                                       enrichr_table_dir,
                                       enrichr_selected_table_dir,
                                       enrichr_fig_dir,
                                       enrichr_dbs,
                                       enrichr_db_labels,
                                       iron_patterns,
                                       padj_thr = 0.05,
                                       top_n_terms = 10,
                                       max_tries = 4,
                                       retry_wait_sec = 20) {
  deg_df <- de_df[
    de_df$is_significant & !is.na(de_df$avg_log2FC) & !is.na(de_df$gene) & de_df$gene != "",
    ,
    drop = FALSE
  ]

  if (nrow(deg_df) == 0) {
    cat("No significant genes available for Enrichr ORA - skipping.\n")
    return(invisible(NULL))
  }

  deg_df$log2FoldChange <- deg_df$avg_log2FC
  deg_df$padj <- deg_df$p_val_adj

  direction_sets <- list(
    "HA_enriched_DLF_states" = deg_df[deg_df$avg_log2FC > 0, , drop = FALSE],
    "Other_DLF_states" = deg_df[deg_df$avg_log2FC < 0, , drop = FALSE]
  )

  for (direction_name in names(direction_sets)) {
    deg_here <- direction_sets[[direction_name]]
    genes_here <- unique(deg_here$gene)
    genes_here <- genes_here[!is.na(genes_here) & genes_here != ""]

    cat("\nRunning Enrichr ORA for ", direction_name, ": ", length(genes_here), " significant genes\n", sep = "")
    if (length(genes_here) < 5) next

    enrich_res <- run_enrichr_with_retry(genes_here, enrichr_dbs, max_tries = max_tries, retry_wait_sec = retry_wait_sec)
    if (is.null(enrich_res)) next

    plot_rows <- list()
    selected_rows <- list()

    for (db in names(enrich_res)) {
      write.csv(as.data.frame(enrich_res[[db]]), file.path(enrichr_table_dir, paste0("enrichr_", direction_name, "_", db, ".csv")), row.names = FALSE)

      plot_df <- prepare_enrichr_iron_plot_df(
        input_df = enrich_res[[db]],
        input_gene_count = length(genes_here),
        db_label = enrichr_db_labels[db],
        iron_patterns = iron_patterns,
        padj_thr = padj_thr,
        top_n = top_n_terms
      )

      if (is.null(plot_df)) next
      plot_rows[[db]] <- plot_df

      for (i in seq_len(nrow(plot_df))) {
        pathway_genes <- unique(strsplit(plot_df$Genes[i], ";", fixed = TRUE)[[1]])
        pathway_genes <- trimws(pathway_genes)
        pathway_genes <- pathway_genes[pathway_genes != ""]

        deg_subset <- deg_here[deg_here$gene %in% pathway_genes, c("gene", "log2FoldChange", "padj"), drop = FALSE]
        if (nrow(deg_subset) == 0) next

        deg_subset$Term <- plot_df$Term[i]
        deg_subset$Adjusted.P.value <- plot_df$Adjusted.P.value[i]
        selected_rows[[paste(db, i, sep = "_")]] <- deg_subset
      }
    }

    if (length(plot_rows) == 0 || length(selected_rows) == 0) {
      cat("No significant iron-related Enrichr terms for ", direction_name, ".\n", sep = "")
      next
    }

    pathway_plot_df <- do.call(rbind, plot_rows)
    pathway_plot_df <- pathway_plot_df[order(pathway_plot_df$gene_ratio, pathway_plot_df$Adjusted.P.value), , drop = FALSE]
    pathway_plot_df$Term <- factor(pathway_plot_df$Term, levels = unique(pathway_plot_df$Term))

    selected_df <- unique(do.call(rbind, selected_rows))
    selected_df <- selected_df[order(selected_df$Term, selected_df$log2FoldChange), , drop = FALSE]
    write.csv(selected_df, file.path(enrichr_selected_table_dir, paste0("selected_iron_related_pathway_genes_", direction_name, ".csv")), row.names = FALSE)

    p_dot <- build_enrichr_dotplot(pathway_plot_df)
    p_genes <- build_selected_pathway_gene_tile_plot(pathway_plot_df, deg_here)

    if (!is.null(p_dot) && !is.null(p_genes)) {
      combined_plot <- (p_dot + p_genes + plot_layout(widths = c(1, 1.6), guides = "collect")) &
        theme(legend.position = "top", legend.box = "horizontal", legend.box.just = "left")

      n_terms <- length(unique(selected_df$Term))
      ggsave(
        file.path(enrichr_fig_dir, paste0("selected_iron_related_pathway_gene_panel_", direction_name, ".png")),
        combined_plot,
        width = ifelse(n_terms == 1, 12, 18),
        height = ifelse(n_terms == 1, 3.9, max(4.8, 1.5 * n_terms)),
        dpi = 300
      )
    }
  }

  invisible(NULL)
}

# plot significant iron-related genes from the FindMarkers table
plot_fibroblast_iron_bubble_heatmap <- function(iron_de_file,
                                                output_table,
                                                output_file,
                                                ha_enriched_group,
                                                other_state_group) {
  iron_de_df <- read.csv(iron_de_file, stringsAsFactors = FALSE)
  iron_sig_genes <- iron_de_df[iron_de_df$is_significant & !is.na(iron_de_df$avg_log2FC), , drop = FALSE]
  iron_sig_genes <- iron_sig_genes[order(-iron_sig_genes$avg_log2FC), , drop = FALSE]

  if (nrow(iron_sig_genes) == 0) {
    cat("No significant iron-related genes found - skipping iron gene bubble heatmap.\n")
    return(invisible(NULL))
  }

  bubble_df <- iron_sig_genes
  bubble_df$comparison <- ifelse(bubble_df$avg_log2FC > 0, ha_enriched_group, other_state_group)
  bubble_df$comparison <- factor(bubble_df$comparison, levels = c(other_state_group, ha_enriched_group))
  bubble_df$gene <- factor(bubble_df$gene, levels = rev(bubble_df$gene))
  bubble_df$neg_log10_padj <- safe_neg_log10(bubble_df$p_val_adj)
  bubble_df$pct_difference <- abs(bubble_df$pct.1 - bubble_df$pct.2)

  max_abs_lfc <- max(abs(bubble_df$avg_log2FC), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) max_abs_lfc <- 1

  write.csv(bubble_df, output_table, row.names = FALSE)

  p_bubble <- ggplot(bubble_df, aes(x = comparison, y = gene)) +
    geom_point(aes(size = neg_log10_padj, fill = avg_log2FC), shape = 21, color = "black", stroke = 0.35) +
    scale_fill_gradient2(low = "#5B8DB8", mid = "white", high = "#B65A5A", midpoint = 0, limits = c(-max_abs_lfc, max_abs_lfc), oob = scales::squish, name = "avg log2FC") +
    scale_size_continuous(name = "-log10 adj. p", range = c(3.2, 9), breaks = pretty(bubble_df$neg_log10_padj, n = 4)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 16, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.line = element_line(linewidth = 1.1, color = "black"),
      axis.ticks = element_line(linewidth = 1.1, color = "black"),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size = 13, color = "black"),
      legend.text = element_text(size = 12, color = "black"),
      legend.position = "right",
      plot.margin = margin(18, 18, 18, 18)
    )

  ggsave(output_file, plot = p_bubble, width = 8.5, height = 8.5, dpi = 600)
  invisible(output_file)
}

# plot TGF-beta genes on the destructive lining fibroblast UMAP
plot_tgf_beta_featureplots <- function(input_object, output_file, reduction_name, features) {
  obj <- readRDS(input_object)

  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj, assay = "RNA")
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)

  features <- features[features %in% rownames(obj)]
  cat("TGF-beta genes found:\n")
  print(features)

  if (length(features) == 0) return(invisible(NULL))

  p_tgf_beta <- FeaturePlot(
    object = obj,
    features = features,
    reduction = reduction_name,
    cols = c("#F3F1EC", "#7A1F2B"),
    order = TRUE,
    raster = FALSE,
    pt.size = 0.35,
    ncol = 4
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

  ggsave(output_file, plot = p_tgf_beta, width = 18, height = 4.2 * ceiling(length(features) / 4), dpi = 600)
  invisible(output_file)
}
