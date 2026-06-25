#########################################################################
### Functional analysis of HA-enriched destructive lining fibroblasts ###
#########################################################################

# functional analysis of the targeted FindMarkers comparison:
# the HA-enriched destructive/MMP3+ lining fibroblast states versus the remaining
# destructive lining fibroblast states. The script saves a bubble heatmap of
# significant iron-related genes, ranked GSEA results/curves, and Enrichr ORA
# panels for significant iron-related terms.

suppressPackageStartupMessages({
  library(enrichR)
  library(fgsea)
  library(ggplot2)
  library(patchwork)
})

setwd("~/Thema_R")
source("src/global_config.R")

# input and output folders
de_res_dir <- file.path(results_dir, "destructive_lining_fibroblast_HA_enriched_state_DE")
de_file <- file.path(de_res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states.csv")
iron_de_file <- file.path(de_res_dir, "FindMarkers_HA_enriched_DLF_states_vs_other_DLF_states_iron_related.csv")

func_res_dir <- file.path(results_dir, "functional_analysis_HA_enriched_destructive_lining_fibroblast_states")
fgsea_table_dir <- file.path(func_res_dir, "fgsea_significant_tables")
fgsea_iron_table_dir <- file.path(func_res_dir, "fgsea_iron_related_significant_tables")
enrichr_table_dir <- file.path(func_res_dir, "enrichr_ora_tables")
enrichr_selected_table_dir <- file.path(func_res_dir, "enrichr_iron_related_gene_tables")

func_fig_dir <- file.path(figures_dir, "functional_analysis_HA_enriched_destructive_lining_fibroblast_states")
bubble_fig_dir <- file.path(func_fig_dir, "bubble_heatmap")
fgsea_curve_dir <- file.path(func_fig_dir, "gsea_curves")
enrichr_fig_dir <- file.path(func_fig_dir, "enrichr_ora")

dir.create(func_res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_iron_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichr_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichr_selected_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bubble_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fgsea_curve_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichr_fig_dir, recursive = TRUE, showWarnings = FALSE)

# set parameters
padj_thr <- 0.05
min_size <- 5
max_size <- 500
fgsea_seed <- 1234
use_multilevel <- TRUE
n_permutations <- 10000
top_n_gsea_curves_per_direction <- 5
top_n_enrichr_terms <- 10
run_enrichr_ora <- TRUE
enrichr_max_tries <- 4
enrichr_retry_wait_sec <- 20

# ranking method:
# avg_log2FC is the default because FindMarkers/Wilcoxon has no DESeq2-like Wald stat.
# positive ranks mean higher in HA-enriched DLF states; negative ranks mean higher in other DLF states.
rank_method <- "avg_log2FC" # options: "avg_log2FC", "signed_log10_pvalue"

ha_enriched_group <- "HA-enriched DLF states"
other_state_group <- "Other DLF states"

# curves selected for biological interpretability rather than by padj alone.
selected_gsea_curves <- c(
  "GOBP_INFLAMMATORY_RESPONSE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "GOBP_HEMOSTASIS",
  "GOBP_REGULATION_OF_COAGULATION",
  "REACTOME_DISSOLUTION_OF_FIBRIN_CLOT",
  "GOBP_RESPONSE_TO_WOUNDING",
  "GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT",
  "GOBP_CILIUM_ORGANIZATION",
  "GOBP_CELL_PROJECTION_ASSEMBLY",
  "GOBP_GOLGI_VESICLE_TRANSPORT"
)

# ranked GSEA uses MSigDB collections
msigdb_species <- "Homo sapiens"
msigdb_collections <- list(
  "GO_BP" = c("GO:BP"),
  "KEGG" = c("CP:KEGG", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"),
  "Reactome" = c("CP:REACTOME")
)

# pathways matching these terms are also saved/plotted separately.
iron_related_patterns <- c(
  "ferropt",
  "\\biron\\b",
  "iron ion",
  "ferric",
  "ferrous",
  "ferritin",
  "transferrin",
  "heme",
  "haem",
  "hemoglobin",
  "hepcidin",
  "oxidative stress"
)

enrichr_dbs <- c(
  "GO_Biological_Process_2025",
  "Reactome_2022",
  "KEGG_2021_Human"
)

enrichr_db_labels <- c(
  "GO_Biological_Process_2025" = "GO BP 2025",
  "Reactome_2022" = "Reactome 2022",
  "KEGG_2021_Human" = "KEGG 2021"
)

sanitize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# load one MSigDB subcollection, supporting both old and new msigdbr arguments.
load_msigdb_subcollection <- function(subcollection) {
  out <- tryCatch(
    msigdbr::msigdbr(
      species = msigdb_species,
      collection = "C2",
      subcollection = subcollection
    ),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) > 0) {
    return(as.data.frame(out))
  }

  out <- tryCatch(
    msigdbr::msigdbr(
      species = msigdb_species,
      category = "C2",
      subcategory = subcollection
    ),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) > 0) {
    return(as.data.frame(out))
  }

  NULL
}

# this separate loader handles GO BP from collection C5.
load_msigdb_gobp <- function() {
  out <- tryCatch(
    msigdbr::msigdbr(
      species = msigdb_species,
      collection = "C5",
      subcollection = "GO:BP"
    ),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) > 0) {
    return(as.data.frame(out))
  }

  out <- tryCatch(
    msigdbr::msigdbr(
      species = msigdb_species,
      category = "C5",
      subcategory = "GO:BP"
    ),
    error = function(e) NULL
  )

  if (!is.null(out) && nrow(out) > 0) {
    return(as.data.frame(out))
  }

  NULL
}

# build GO BP, KEGG and Reactome pathways from msigdbr.
load_msigdb_pathways <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "Package 'msigdbr' is required for GO BP / KEGG / Reactome fgsea. ",
      "Install it in thema_r_env and rerun this script."
    )
  }

  pathway_tables <- list()

  gobp_df <- load_msigdb_gobp()
  if (!is.null(gobp_df) && nrow(gobp_df) > 0) {
    gobp_df$database <- "GO_BP"
    pathway_tables[["GO_BP"]] <- gobp_df
  }

  for (db in c("KEGG", "Reactome")) {
    db_tables <- list()
    for (subcollection in msigdb_collections[[db]]) {
      df <- load_msigdb_subcollection(subcollection)
      if (!is.null(df) && nrow(df) > 0) {
        df$database <- db
        db_tables[[subcollection]] <- df
      }
    }

    if (length(db_tables) > 0) {
      pathway_tables[[db]] <- do.call(rbind, db_tables)
    }
  }

  if (length(pathway_tables) == 0) {
    stop("No MSigDB pathways loaded. Check msigdbr installation and collection names.")
  }

  pathway_df <- do.call(rbind, pathway_tables)
  required_cols <- c("database", "gs_name", "gene_symbol")
  missing_cols <- setdiff(required_cols, colnames(pathway_df))
  if (length(missing_cols) > 0) {
    stop("Missing expected msigdbr columns: ", paste(missing_cols, collapse = ", "))
  }

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

# flag pathway names that are directly related to iron, heme, ferroptosis or oxidative stress.
is_iron_related_pathway <- function(pathway_name) {
  pathway_text <- tolower(pathway_name)
  pattern <- paste(iron_related_patterns, collapse = "|")
  grepl(pattern, pathway_text, perl = TRUE)
}

clean_enrichr_terms <- function(df) {
  df$Term <- gsub("\\s+R-HSA-[0-9]+$", "", df$Term)
  df$Term <- gsub("\\s+KEGG_[0-9]+$", "", df$Term)
  df$Term <- gsub("\\s*\\([^\\)]+\\)", "", df$Term)
  df$Term <- trimws(df$Term)
  df
}

format_selected_term_label <- function(x) {
  x <- gsub(" \\(", "\n(", x, fixed = TRUE)
  vapply(strsplit(x, "\n", fixed = TRUE), function(parts) {
    first_line_words <- strsplit(parts[1], " ", fixed = TRUE)[[1]]
    wrapped_first <- paste(
      tapply(first_line_words, ceiling(seq_along(first_line_words) / 5), paste, collapse = " "),
      collapse = "\n"
    )
    if (length(parts) > 1) {
      paste(c(wrapped_first, parts[-1]), collapse = "\n")
    } else {
      wrapped_first
    }
  }, character(1))
}

prepare_enrichr_iron_plot_df <- function(input_df,
                                         input_gene_count,
                                         db_label,
                                         top_n = top_n_enrichr_terms) {
  df <- as.data.frame(input_df)
  if (nrow(df) == 0) return(NULL)

  df <- clean_enrichr_terms(df)

  if (nrow(df) == 0) return(NULL)

  term_regex <- paste(iron_related_patterns, collapse = "|")
  df$num_genes <- sapply(strsplit(df$Genes, ";"), length)
  df$gene_ratio <- df$num_genes / max(1, input_gene_count)
  df$is_significant <- !is.na(df$Adjusted.P.value) & df$Adjusted.P.value < padj_thr
  df <- df[grepl(term_regex, df$Term, ignore.case = TRUE) & df$is_significant, , drop = FALSE]

  if (nrow(df) == 0) return(NULL)

  df <- df[order(df$gene_ratio, df$Adjusted.P.value), , drop = FALSE]
  df <- head(df, top_n)
  df$Database <- db_label
  df$Term <- paste0(df$Term, " (", df$Database, ")")
  df$Term <- format_selected_term_label(df$Term)
  df$Term <- factor(df$Term, levels = unique(df$Term))
  df
}

build_enrichr_dotplot <- function(plot_df,
                                  color_palette = c("#8B3E2F", "coral")) {
  size_breaks <- unique(round(pretty(plot_df$num_genes)))
  size_breaks <- size_breaks[size_breaks >= min(plot_df$num_genes) & size_breaks <= max(plot_df$num_genes)]
  color_breaks <- unique(signif(seq(min(plot_df$Adjusted.P.value), max(plot_df$Adjusted.P.value), length.out = 3), 3))

  ggplot(plot_df, aes(x = gene_ratio, y = Term)) +
    geom_point(aes(size = num_genes, color = Adjusted.P.value), alpha = 0.9) +
    scale_color_continuous(
      low = color_palette[1],
      high = color_palette[2],
      name = "adj p-value",
      breaks = color_breaks
    ) +
    scale_size_continuous(
      name = "Genes",
      range = c(3, 9),
      breaks = size_breaks
    ) +
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
      size = guide_legend(
        title.position = "top",
        nrow = 1,
        byrow = TRUE,
        order = 1,
        keywidth = grid::unit(0.18, "cm"),
        keyheight = grid::unit(0.32, "cm")
      ),
      color = guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(2.8, "cm"),
        barheight = grid::unit(0.35, "cm"),
        label.position = "bottom",
        label = scales::label_number(accuracy = 0.01),
        title.theme = element_text(margin = margin(b = 4)),
        label.theme = element_text(size = 8, margin = margin(t = 4)),
        order = 2
      )
    )
}

build_selected_pathway_gene_tile_plot <- function(pathway_df,
                                                  deg_df) {
  term_levels <- levels(pathway_df$Term)
  sig_df <- pathway_df[pathway_df$is_significant, , drop = FALSE]
  if (nrow(sig_df) == 0) return(NULL)

  gene_tiles <- do.call(
    rbind,
    lapply(seq_len(nrow(sig_df)), function(i) {
      genes <- unique(strsplit(sig_df$Genes[i], ";", fixed = TRUE)[[1]])
      genes <- trimws(genes)
      genes <- genes[genes != ""]

      data.frame(
        Term = rep(as.character(sig_df$Term[i]), length(genes)),
        gene = genes,
        gene_index = seq_along(genes),
        stringsAsFactors = FALSE
      )
    })
  )

  gene_tiles <- merge(
    gene_tiles,
    unique(deg_df[, c("gene", "log2FoldChange")]),
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  gene_tiles$Term <- factor(gene_tiles$Term, levels = term_levels)
  gene_tiles$segment_width <- 1.2
  gene_tiles <- gene_tiles[order(gene_tiles$Term, gene_tiles$gene_index), , drop = FALSE]
  gene_tiles$x_center <- ave(
    gene_tiles$segment_width,
    gene_tiles$Term,
    FUN = function(x) cumsum(x) - (x / 2)
  )

  max_genes <- max(as.numeric(table(gene_tiles$Term)))
  max_abs_lfc <- max(abs(gene_tiles$log2FoldChange), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) {
    max_abs_lfc <- 1
  }

  ggplot(gene_tiles, aes(y = Term, x = segment_width, fill = log2FoldChange)) +
    geom_col(width = 0.78, position = "stack", color = "black", linewidth = 0.35) +
    geom_text(
      data = gene_tiles,
      aes(x = x_center, y = Term, label = gene),
      angle = 0,
      size = 3.2,
      color = "black",
      inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      limits = c(-max_abs_lfc, max_abs_lfc),
      oob = scales::squish,
      name = "log2FC"
    ) +
    scale_x_continuous(
      limits = c(0, max_genes),
      breaks = 0:max_genes,
      labels = rep("", max_genes + 1),
      expand = c(0, 0)
    ) +
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
    guides(
      fill = guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(2.6, "cm"),
        barheight = grid::unit(0.35, "cm"),
        label.theme = element_text(size = 8, margin = margin(t = 4))
      )
    )
}

# enrichr is an online service, so retry transient HTTP/curl failures instead of
# stopping the whole functional analysis after GSEA and bubble plots are done.
run_enrichr_with_retry <- function(genes, databases) {
  old_timeout <- getOption("timeout")
  options(timeout = max(300, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  last_error <- NULL

  for (attempt in seq_len(enrichr_max_tries)) {
    cat(
      "Enrichr attempt ", attempt, "/", enrichr_max_tries,
      " with ", length(genes), " genes...\n",
      sep = ""
    )

    res <- tryCatch(
      enrichr(genes, databases),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(res) && length(res) > 0) {
      return(res)
    }

    cat("Enrichr attempt failed: ", last_error, "\n", sep = "")
    if (attempt < enrichr_max_tries) {
      cat("Waiting ", enrichr_retry_wait_sec, " seconds before retrying...\n", sep = "")
      Sys.sleep(enrichr_retry_wait_sec)
    }
  }

  cat("Enrichr failed after all retries - skipping this direction.\n")
  if (!is.null(last_error)) {
    cat("Last Enrichr error: ", last_error, "\n", sep = "")
  }

  NULL
}

# build the ranked vector used by fgsea.
build_rank_vector <- function(de_df, rank_method = rank_method) {
  if (rank_method == "avg_log2FC") {
    rank_df <- de_df[
      !is.na(de_df$avg_log2FC),
      c("gene", "avg_log2FC"),
      drop = FALSE
    ]
    colnames(rank_df) <- c("gene", "rank_value")
  } else if (rank_method == "signed_log10_pvalue") {
    rank_df <- de_df[
      !is.na(de_df$avg_log2FC) &
        !is.na(de_df$p_val),
      c("gene", "avg_log2FC", "p_val"),
      drop = FALSE
    ]

    positive_p <- rank_df$p_val[rank_df$p_val > 0]
    min_positive_p <- min(positive_p, na.rm = TRUE)
    if (!is.finite(min_positive_p)) {
      min_positive_p <- .Machine$double.xmin
    }
    rank_df$p_val_rank <- rank_df$p_val
    rank_df$p_val_rank[rank_df$p_val_rank <= 0] <- min_positive_p * 0.1
    rank_df$rank_value <- sign(rank_df$avg_log2FC) * -log10(rank_df$p_val_rank)
    rank_df <- rank_df[, c("gene", "rank_value"), drop = FALSE]
  } else {
    stop("Unsupported rank_method: ", rank_method)
  }

  rank_df$gene <- toupper(as.character(rank_df$gene))
  rank_df <- rank_df[!is.na(rank_df$gene) & rank_df$gene != "", , drop = FALSE]
  rank_df <- rank_df[!is.na(rank_df$rank_value), , drop = FALSE]
  rank_df <- rank_df[rank_df$rank_value != 0, , drop = FALSE]

  # fgsea needs one value per gene, so keep the strongest duplicate if present.
  # this absolute ordering is only for duplicate removal; the final ranking stays signed.
  rank_df <- rank_df[order(abs(rank_df$rank_value), decreasing = TRUE), , drop = FALSE]
  rank_df <- rank_df[!duplicated(rank_df$gene), , drop = FALSE]

  ranks <- rank_df$rank_value
  names(ranks) <- rank_df$gene

  # ranked GSEA should use all genes, not only significant markers.
  sort(ranks, decreasing = TRUE)
}

run_fgsea <- function(pathways, ranks) {
  set.seed(fgsea_seed)

  if (use_multilevel) {
    fgseaMultilevel(
      pathways = pathways,
      stats = ranks,
      minSize = min_size,
      maxSize = max_size
    )
  } else {
    fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = min_size,
      maxSize = max_size,
      nperm = n_permutations
    )
  }
}

add_pathway_metadata <- function(fgsea_res, pathway_info) {
  fgsea_df <- as.data.frame(fgsea_res)
  if (nrow(fgsea_df) == 0) return(fgsea_df)

  fgsea_df$pathway_id <- as.character(fgsea_df$pathway)
  fgsea_df$database <- pathway_info[fgsea_df$pathway_id, "database"]
  fgsea_df$pathway_name <- pathway_info[fgsea_df$pathway_id, "gs_name"]
  fgsea_df$is_iron_related <- is_iron_related_pathway(fgsea_df$pathway_name)
  fgsea_df$leadingEdge <- vapply(fgsea_df$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))

  fgsea_df <- fgsea_df[
    order(fgsea_df$padj, -abs(fgsea_df$NES)),
    c(
      "database", "pathway_name", "pathway_id", "pval", "padj", "ES", "NES",
      "size", "leadingEdge", "is_iron_related"
    ),
    drop = FALSE
  ]

  fgsea_df
}

plot_gsea_curve <- function(pathway_id, pathway_name, ranks, pathways, output_file) {
  p <- plotEnrichment(pathways[[pathway_id]], ranks) +
    labs(
      title = pathway_name,
      x = "Ranked genes",
      y = "Enrichment score"
    ) +
    theme_classic(base_size = 16) +
    theme(
      axis.text = element_text(size = 13, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(18, 18, 18, 18)
    )

  ggsave(output_file, plot = p, width = 8.5, height = 5.5, dpi = 600)
}

plot_significant_iron_gene_bubble_heatmap <- function() {
  if (!file.exists(iron_de_file)) {
    cat("Missing iron-related DE table - skipping iron gene bubble heatmap.\n")
    return(invisible(NULL))
  }

  iron_de_df <- read.csv(iron_de_file, stringsAsFactors = FALSE)
  required_cols <- c("gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2", "direction", "is_significant")
  missing_cols <- setdiff(required_cols, colnames(iron_de_df))
  if (length(missing_cols) > 0) {
    cat("Iron-related DE table lacks required columns - skipping iron gene bubble heatmap.\n")
    return(invisible(NULL))
  }

  iron_sig_genes <- iron_de_df[
    iron_de_df$is_significant &
      !is.na(iron_de_df$avg_log2FC),
    ,
    drop = FALSE
  ]
  iron_sig_genes <- iron_sig_genes[order(-iron_sig_genes$avg_log2FC), , drop = FALSE]

  if (nrow(iron_sig_genes) == 0) {
    cat("No significant iron-related genes found - skipping iron gene bubble heatmap.\n")
    return(invisible(NULL))
  }

  cat("Plotting significant iron-related genes from FindMarkers table as bubble heatmap...\n")

  bubble_df <- iron_sig_genes
  bubble_df$comparison <- ifelse(
    bubble_df$avg_log2FC > 0,
    ha_enriched_group,
    other_state_group
  )
  bubble_df$comparison <- factor(bubble_df$comparison, levels = c(other_state_group, ha_enriched_group))
  bubble_df$gene <- factor(bubble_df$gene, levels = rev(bubble_df$gene))
  bubble_df$neg_log10_padj <- -log10(bubble_df$p_val_adj)

  if (any(is.infinite(bubble_df$neg_log10_padj), na.rm = TRUE)) {
    max_finite <- max(bubble_df$neg_log10_padj[is.finite(bubble_df$neg_log10_padj)], na.rm = TRUE)
    if (!is.finite(max_finite)) max_finite <- 0
    bubble_df$neg_log10_padj[is.infinite(bubble_df$neg_log10_padj)] <- max_finite + 1
  }

  bubble_df$pct_difference <- abs(bubble_df$pct.1 - bubble_df$pct.2)
  max_abs_lfc <- max(abs(bubble_df$avg_log2FC), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) {
    max_abs_lfc <- 1
  }

  output_table <- file.path(
    func_res_dir,
    "significant_iron_related_genes_bubble_heatmap_input_HA_enriched_DLF_states_vs_other_states.csv"
  )
  write.csv(bubble_df, output_table, row.names = FALSE)

  p_bubble <- ggplot(bubble_df, aes(x = comparison, y = gene)) +
    geom_point(aes(size = neg_log10_padj, fill = avg_log2FC), shape = 21, color = "black", stroke = 0.35) +
    scale_fill_gradient2(
      low = "#5B8DB8",
      mid = "white",
      high = "#B65A5A",
      midpoint = 0,
      limits = c(-max_abs_lfc, max_abs_lfc),
      oob = scales::squish,
      name = "avg log2FC"
    ) +
    scale_size_continuous(
      name = "-log10 adj. p",
      range = c(3.2, 9),
      breaks = pretty(bubble_df$neg_log10_padj, n = 4)
    ) +
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

  output_file <- file.path(
    bubble_fig_dir,
    "bubble_heatmap_significant_iron_related_genes_HA_enriched_DLF_states_vs_other_states.png"
  )

  ggsave(output_file, plot = p_bubble, width = 8.5, height = 8.5, dpi = 600)
  cat("Saved significant iron gene bubble heatmap: ", output_file, "\n", sep = "")
  invisible(output_file)
}

run_enrichr_iron_related_ora <- function(de_df) {
  if (!run_enrichr_ora) {
    cat("run_enrichr_ora is FALSE - skipping Enrichr ORA.\n")
    return(invisible(NULL))
  }

  required_cols <- c("gene", "avg_log2FC", "p_val_adj", "is_significant")
  missing_cols <- setdiff(required_cols, colnames(de_df))
  if (length(missing_cols) > 0) {
    cat("FindMarkers table lacks columns needed for Enrichr ORA - skipping.\n")
    return(invisible(NULL))
  }

  deg_df <- de_df[
    de_df$is_significant &
      !is.na(de_df$avg_log2FC) &
      !is.na(de_df$gene) &
      de_df$gene != "",
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

    if (length(genes_here) < 5) {
      cat("Skipping Enrichr ORA for ", direction_name, ": fewer than 5 genes.\n", sep = "")
      next
    }

    enrich_res <- run_enrichr_with_retry(genes_here, enrichr_dbs)
    if (is.null(enrich_res)) {
      next
    }

    plot_rows <- list()
    selected_rows <- list()

    for (db in names(enrich_res)) {
      raw_table_file <- file.path(
        enrichr_table_dir,
        paste0("enrichr_", direction_name, "_", db, ".csv")
      )
      write.csv(as.data.frame(enrich_res[[db]]), raw_table_file, row.names = FALSE)

      plot_df <- prepare_enrichr_iron_plot_df(
        input_df = enrich_res[[db]],
        input_gene_count = length(genes_here),
        db_label = enrichr_db_labels[db],
        top_n = top_n_enrichr_terms
      )

      if (is.null(plot_df)) next

      plot_rows[[db]] <- plot_df

      for (i in seq_len(nrow(plot_df))) {
        pathway_genes <- unique(strsplit(plot_df$Genes[i], ";", fixed = TRUE)[[1]])
        pathway_genes <- trimws(pathway_genes)
        pathway_genes <- pathway_genes[pathway_genes != ""]
        if (length(pathway_genes) == 0) next

        deg_subset <- deg_here[
          deg_here$gene %in% pathway_genes,
          c("gene", "log2FoldChange", "padj"),
          drop = FALSE
        ]
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

    selected_df <- do.call(rbind, selected_rows)
    selected_df <- unique(selected_df)
    selected_df <- selected_df[order(selected_df$Term, selected_df$log2FoldChange), , drop = FALSE]

    selected_table_file <- file.path(
      enrichr_selected_table_dir,
      paste0("selected_iron_related_pathway_genes_", direction_name, ".csv")
    )
    write.csv(selected_df, selected_table_file, row.names = FALSE)

    p_dot <- build_enrichr_dotplot(pathway_plot_df)
    p_genes <- build_selected_pathway_gene_tile_plot(
      pathway_df = pathway_plot_df,
      deg_df = deg_here
    )

    if (!is.null(p_dot) && !is.null(p_genes)) {
      combined_plot <- (p_dot + p_genes + plot_layout(widths = c(1, 1.6), guides = "collect")) &
        theme(
          legend.position = "top",
          legend.box = "horizontal",
          legend.box.just = "left",
          legend.spacing.x = grid::unit(0.45, "cm"),
          legend.spacing.y = grid::unit(0.08, "cm"),
          legend.margin = margin(0, 0, 0, 0)
        )

      plot_file <- file.path(
        enrichr_fig_dir,
        paste0("selected_iron_related_pathway_gene_panel_", direction_name, ".png")
      )

      n_terms <- length(unique(selected_df$Term))
      panel_width <- if (n_terms == 1) 12 else 18
      panel_height <- if (n_terms == 1) 3.9 else max(4.8, 1.5 * n_terms)

      ggsave(plot_file, plot = combined_plot, width = panel_width, height = panel_height, dpi = 300)
      cat("Saved Enrichr iron-related panel: ", plot_file, "\n", sep = "")
    }
  }

  invisible(NULL)
}

# load inputs
if (!file.exists(de_file)) {
  stop("Missing FindMarkers result file: ", de_file)
}

cat("Loading targeted FindMarkers table:\n")
cat(de_file, "\n")
de_df <- read.csv(de_file, stringsAsFactors = FALSE)

required_cols <- c("gene", "avg_log2FC", "p_val")
missing_cols <- setdiff(required_cols, colnames(de_df))
if (length(missing_cols) > 0) {
  stop("Missing required FindMarkers columns: ", paste(missing_cols, collapse = ", "))
}

cat("Ranking method:", rank_method, "\n")
ranks <- build_rank_vector(de_df, rank_method = rank_method)
cat("Ranked genes:", length(ranks), "\n")
cat("Top ranked genes toward HA-enriched DLF states:\n")
print(head(ranks, 10))
cat("Bottom ranked genes toward other DLF states:\n")
print(tail(ranks, 10))

pathway_data <- load_msigdb_pathways()
pathways <- pathway_data$pathways
pathway_info <- pathway_data$pathway_info

cat("Pathways loaded:\n")
print(table(pathway_info$database))

fgsea_res <- run_fgsea(pathways, ranks)
fgsea_df <- add_pathway_metadata(fgsea_res, pathway_info)

all_file <- file.path(func_res_dir, "fgsea_all_results_HA_enriched_DLF_states_vs_other_DLF_states.csv")
write.csv(fgsea_df, all_file, row.names = FALSE)

sig_df <- fgsea_df[
  !is.na(fgsea_df$padj) &
    fgsea_df$padj < padj_thr,
  ,
  drop = FALSE
]

sig_file <- file.path(fgsea_table_dir, "fgsea_significant_HA_enriched_DLF_states_vs_other_DLF_states.csv")
write.csv(sig_df, sig_file, row.names = FALSE)

iron_sig_df <- sig_df[sig_df$is_iron_related, , drop = FALSE]
iron_file <- file.path(fgsea_iron_table_dir, "fgsea_iron_related_significant_HA_enriched_DLF_states_vs_other_DLF_states.csv")
write.csv(iron_sig_df, iron_file, row.names = FALSE)

cat("Significant pathways at FDR < ", padj_thr, ": ", nrow(sig_df), "\n", sep = "")
cat("Significant iron/heme/ferroptosis/oxidative-stress pathways: ", nrow(iron_sig_df), "\n", sep = "")

cat("\nTop significant pathways:\n")
print(head(sig_df[, c("database", "pathway_name", "padj", "NES", "size", "is_iron_related"), drop = FALSE], 20))

if (nrow(iron_sig_df) > 0) {
  cat("\nIron-related significant pathways:\n")
  print(iron_sig_df[, c("database", "pathway_name", "padj", "NES", "size"), drop = FALSE])
}

# plot GSEA curves selected for biological interpretability.
curve_rows <- sig_df[sig_df$pathway_name %in% selected_gsea_curves, , drop = FALSE]

missing_selected_curves <- setdiff(selected_gsea_curves, curve_rows$pathway_name)
if (length(missing_selected_curves) > 0) {
  cat("Selected GSEA curves not significant/found and therefore not plotted:\n")
  print(missing_selected_curves)
}

curve_rows <- curve_rows[!duplicated(curve_rows$pathway_id), , drop = FALSE]

write.csv(curve_rows, file.path(func_res_dir, "fgsea_pathways_selected_for_curve_plotting.csv"), row.names = FALSE)

if (nrow(curve_rows) > 0) {
  for (i in seq_len(nrow(curve_rows))) {
    pathway_id <- curve_rows$pathway_id[i]
    pathway_name <- curve_rows$pathway_name[i]
    safe_pathway <- sanitize_label(pathway_name)
    direction <- ifelse(curve_rows$NES[i] > 0, "HA_enriched_states", "other_states")

    plot_gsea_curve(
      pathway_id = pathway_id,
      pathway_name = pathway_name,
      ranks = ranks,
      pathways = pathways,
      output_file = file.path(
        fgsea_curve_dir,
        paste0("gsea_curve_", direction, "_", safe_pathway, ".png")
      )
    )
  }
}

# add a companion bubble heatmap using significant iron-related DE genes from script 23.
plot_significant_iron_gene_bubble_heatmap()

# add Enrichr over-representation panels for significant iron-related terms.
run_enrichr_iron_related_ora(de_df)

# save a compact summary table for quick inspection.
summary_df <- data.frame(
  comparison = "HA-enriched DLF states vs other DLF states",
  rank_method = rank_method,
  n_ranked_genes = length(ranks),
  n_significant_pathways = nrow(sig_df),
  n_iron_related_significant_pathways = nrow(iron_sig_df),
  stringsAsFactors = FALSE
)

write.csv(summary_df, file.path(func_res_dir, "fgsea_summary_HA_enriched_DLF_states_vs_other_DLF_states.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("Functional analysis for HA-enriched destructive lining fibroblast states complete.\n")
cat("All results             : ", all_file, "\n")
cat("Significant table       : ", sig_file, "\n")
cat("Iron-related significant: ", iron_file, "\n")
cat("Bubble heatmap figures  : ", bubble_fig_dir, "\n")
cat("GSEA curves             : ", fgsea_curve_dir, "\n")
cat("Enrichr ORA panels      : ", enrichr_fig_dir, "\n")
cat("Interpretation          : positive NES = enriched toward HA-enriched DLF states; negative NES = enriched toward other DLF states.\n")
cat("============================================================\n")
