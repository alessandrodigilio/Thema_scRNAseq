# -----------------------------
# Enrichr over-representation analysis on pseudobulk DEGs
# -----------------------------

suppressPackageStartupMessages({
  library(enrichR)
  library(ggplot2)
  library(patchwork)
})

setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Input and output folders
PSEUDOBULK_RES_DIR <- file.path(RESULTS_DIR, "pseudobulk_deseq2")
ENRICHR_RES_DIR <- file.path(RESULTS_DIR, "enrichr_pseudobulk_HA_vs_other")
ENRICHR_TABLE_DIR <- file.path(ENRICHR_RES_DIR, "tables")
ENRICHR_FIG_DIR <- file.path(FIGURES_DIR, "enrichr_pseudobulk_HA_vs_other")
ENRICHR_PATHWAY_FIG_DIR <- file.path(ENRICHR_FIG_DIR, "selected_pathway_gene_barplots")
ENRICHR_PATHWAY_TABLE_DIR <- file.path(ENRICHR_RES_DIR, "selected_pathway_gene_tables")

dir.create(ENRICHR_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ENRICHR_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ENRICHR_PATHWAY_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ENRICHR_PATHWAY_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Enrichment parameters
PADJ_CUTOFF <- 0.05
MIN_DEGS_PER_CELLTYPE <- 30
TOP_N_TERMS <- 10
RUN_STANDARD_ENRICHR_DOTPLOTS <- FALSE
RUN_SELECTED_PATHWAY_GENE_BARPLOTS <- TRUE
SELECTED_PATHWAY_CELLTYPES <- c(
  "Lining fibroblasts (PRG4+)",
  "Destructive lining fibroblasts (MMP3+)"
)
IRON_RELATED_TERM_PATTERNS <- c(
  "Ferroptosis",
  "Regulation of Ferroptosis",
  "Negative Regulation of Ferroptosis",
  "Positive Regulation of Ferroptosis",
  "Iron Ion Transmembrane Transport",
  "Iron Ion Transport",
  "Iron Uptake And Transport",
  "Iron uptake and transport",
  "Detoxification Of Reactive Oxygen Species",
  "Reactive Oxygen Species",
  "Heme",
  "Haem",
  "Iron"
)

TGF_BETA_TERM_PATTERNS <- c(
  "TGF",
  "TGFB",
  "TGF-beta",
  "Transforming Growth Factor",
  "Transforming Growth Factor Beta",
  "SMAD",
  "SMAD2",
  "SMAD3",
  "SMAD4"
)

# Exclude very broad or low-information terms from the plots.
# Patterns are matched after the Reactome / KEGG suffix cleanup.
EXCLUDED_TERM_PATTERNS <- c(
  "Immune System",
  "Cytokine Signaling In Immune System",
  "Cellular Responses To Stimuli",
  "Cellular Responses To Stress",
  "Cellular Response to Heat",
  "Response to Unfolded Protein",
  "HSF1 Activation",
  "HSF1-dependent Transactivation",
  "Regulation Of HSF1-mediated Heat Shock Response",
  "Attenuation Phase",
  "Translation",
  "Mitochondrial Translation",
  "Mitochondrial Translation Initiation",
  "Mitochondrial Translation Elongation",
  "Mitochondrial Translation Termination",
  "Mitochondrial Gene Expression",
  "DNA-templated Transcription",
  "Transcription by RNA Polymerase II",
  "Positive Regulation of Transcription by RNA Polymerase II",
  "RNA polymerase",
  "Ribosome Biogenesis",
  "Ribosome biogenesis in eukaryotes",
  "Ribonucleoprotein Complex Biogenesis",
  "rRNA Modification In Nucleus And Cytosol",
  "Maturation of SSU-rRNA",
  "Protein processing in endoplasmic reticulum",
  "Antigen processing and presentation",
  "Lipid and atherosclerosis",
  "Pathways in cancer",
  "MAPK signaling pathway",
  "Estrogen signaling pathway",
  "Toxoplasmosis",
  "Legionellosis",
  "Measles",
  "Neutrophil Degranulation"
)

DBS <- c(
  "GO_Biological_Process_2025",
  "Reactome_2022",
  "KEGG_2021_Human"
)

DB_LABELS <- c(
  "GO_Biological_Process_2025" = "GO BP 2025",
  "Reactome_2022" = "Reactome 2022",
  "KEGG_2021_Human" = "KEGG 2021"
)

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

# Clean Enrichr output and keep the top terms used by the dotplot.
prepare_enrichr_plot_df <- function(input_df,
                                    input_gene_count,
                                    top_n = TOP_N_TERMS) {
  df <- as.data.frame(input_df)

  if (nrow(df) == 0) {
    return(NULL)
  }

  df <- clean_enrichr_terms(df)
  if (length(EXCLUDED_TERM_PATTERNS) > 0) {
    exclude_regex <- paste(EXCLUDED_TERM_PATTERNS, collapse = "|")
    df <- df[!grepl(exclude_regex, df$Term, ignore.case = FALSE), , drop = FALSE]
  }

  if (nrow(df) == 0) {
    return(NULL)
  }

  df$num_genes <- sapply(strsplit(df$Genes, ";"), length)
  df$gene_ratio <- df$num_genes / input_gene_count
  df$is_significant <- !is.na(df$Adjusted.P.value) & df$Adjusted.P.value < PADJ_CUTOFF

  if (sum(df$is_significant) == 0) {
    return(NULL)
  }

  df <- df[df$is_significant, , drop = FALSE]
  df <- df[order(df$Adjusted.P.value), , drop = FALSE]
  df <- head(df, top_n)
  df$Term <- sapply(strsplit(df$Term, " "), function(x) {
    paste(tapply(x, ceiling(seq_along(x) / 5), paste, collapse = " "), collapse = "\n")
  })
  df <- df[order(df$gene_ratio), , drop = FALSE]
  df$Term <- factor(df$Term, levels = df$Term)
  df
}

# Build the standard Enrichr dotplot.
build_enrichr_dotplot <- function(plot_df,
                                  title,
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
    labs(
      title = NULL,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
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

# Create and save the standard dotplot from raw Enrichr results.
plot_enrichr_dotplot <- function(input_df,
                                 title,
                                 output_file,
                                 input_gene_count,
                                 top_n = TOP_N_TERMS,
                                 color_palette = c("#8B3E2F", "coral")) {
  plot_df <- prepare_enrichr_plot_df(
    input_df = input_df,
    input_gene_count = input_gene_count,
    top_n = top_n
  )

  if (is.null(plot_df)) {
    message("No significant enrichment for ", title, " - skipping plot.")
    return(invisible(NULL))
  }

  p <- build_enrichr_dotplot(
    plot_df = plot_df,
    title = title,
    color_palette = color_palette
  )

  print(p)
  ggsave(output_file, plot = p, width = 9, height = 6, dpi = 300)
  message("Saved enrichment dotplot: ", output_file)
  invisible(p)
}

# Build the gene-tile panel using the same style as PEGASO.
build_selected_pathway_gene_tile_plot <- function(pathway_df,
                                                  deg_df,
                                                  title) {
  term_levels <- levels(pathway_df$Term)
  sig_df <- pathway_df[pathway_df$is_significant, , drop = FALSE]

  if (nrow(sig_df) == 0) {
    return(NULL)
  }

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
        gene_row = ifelse(seq_along(genes) %% 2 == 0, 2, 1),
        gene_col = ceiling(seq_along(genes) / 2),
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
  gene_tiles$Term_index <- as.numeric(gene_tiles$Term)
  gene_tiles$y_center <- gene_tiles$Term_index + ifelse(gene_tiles$gene_row == 1, -0.22, 0.22)
  gene_tiles$x_center <- gene_tiles$gene_col

  max_cols <- max(gene_tiles$gene_col, na.rm = TRUE)
  max_abs_lfc <- max(abs(gene_tiles$log2FoldChange), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) {
    max_abs_lfc <- 1
  }

  ggplot(gene_tiles, aes(x = x_center, y = y_center, fill = log2FoldChange)) +
    geom_tile(width = 0.96, height = 0.38, color = "black", linewidth = 0.35) +
    geom_text(
      data = gene_tiles,
      aes(x = x_center, y = y_center, label = gene),
      angle = 0,
      size = 3.4,
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
      limits = c(0.5, max_cols + 0.5),
      breaks = seq_len(max_cols),
      labels = rep("", max_cols),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(term_levels),
      labels = rep("", length(term_levels)),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 14, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
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

# Extract genes from selected Enrichr terms and save a combined Enrichr panel.
run_selected_pathway_gene_barplots <- function(enrich_res,
                                               deg_df,
                                               cell_type_label,
                                               safe_ct,
                                               term_patterns,
                                               panel_name,
                                               gene_label) {
  if (!(cell_type_label %in% SELECTED_PATHWAY_CELLTYPES)) {
    return(invisible(NULL))
  }

  plot_rows <- list()
  selected_rows <- list()

  for (db in names(enrich_res)) {
    df <- as.data.frame(enrich_res[[db]])
    if (nrow(df) == 0) next

    df <- clean_enrichr_terms(df)
    term_regex <- paste(term_patterns, collapse = "|")
    df$num_genes <- sapply(strsplit(df$Genes, ";"), length)
    df$gene_ratio <- df$num_genes / max(1, nrow(unique(deg_df["gene"])))
    df$is_significant <- !is.na(df$Adjusted.P.value) & df$Adjusted.P.value < PADJ_CUTOFF
    df <- df[grepl(term_regex, df$Term, ignore.case = TRUE) & df$is_significant, , drop = FALSE]
    if (nrow(df) == 0) next

    df <- df[order(df$gene_ratio, df$Adjusted.P.value), , drop = FALSE]
    df$Database <- DB_LABELS[db]
    df$Term <- paste0(df$Term, " (", df$Database, ")")
    df$Term <- format_selected_term_label(df$Term)
    df$Term <- factor(df$Term, levels = unique(df$Term))
    plot_rows[[db]] <- df

    for (i in seq_len(nrow(df))) {
      genes_here <- unique(strsplit(df$Genes[i], ";", fixed = TRUE)[[1]])
      genes_here <- trimws(genes_here)
      genes_here <- genes_here[genes_here != ""]
      if (length(genes_here) == 0) next

      deg_subset <- deg_df[deg_df$gene %in% genes_here, c("gene", "log2FoldChange", "padj"), drop = FALSE]
      if (nrow(deg_subset) == 0) next

      deg_subset$Term <- df$Term[i]
      deg_subset$Adjusted.P.value <- df$Adjusted.P.value[i]
      selected_rows[[paste(db, i, sep = "_")]] <- deg_subset
    }
  }

  if (length(selected_rows) == 0 || length(plot_rows) == 0) {
    message("No ", panel_name, " pathway hits for ", cell_type_label)
    return(invisible(NULL))
  }

  pathway_plot_df <- do.call(rbind, plot_rows)
  pathway_plot_df <- pathway_plot_df[order(pathway_plot_df$gene_ratio, pathway_plot_df$Adjusted.P.value), , drop = FALSE]
  pathway_plot_df$Term <- factor(pathway_plot_df$Term, levels = unique(pathway_plot_df$Term))

  selected_df <- do.call(rbind, selected_rows)
  selected_df <- unique(selected_df)
  selected_df <- selected_df[order(selected_df$Term, selected_df$log2FoldChange), , drop = FALSE]

  table_file <- file.path(
    ENRICHR_PATHWAY_TABLE_DIR,
    paste0("selected_", panel_name, "_pathway_genes_", safe_ct, ".csv")
  )
  write.csv(selected_df, table_file, row.names = FALSE)
  message("Saved selected pathway gene table: ", table_file)

  p_dot <- build_enrichr_dotplot(
    plot_df = pathway_plot_df,
    title = paste0(cell_type_label, ": ", panel_name, " pathways")
  )
  p_genes <- build_selected_pathway_gene_tile_plot(
    pathway_df = pathway_plot_df,
    deg_df = deg_df,
    title = gene_label
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
      ENRICHR_PATHWAY_FIG_DIR,
      paste0("selected_", panel_name, "_pathway_gene_panel_", safe_ct, ".png")
    )
    n_terms <- length(unique(selected_df$Term))
    print(combined_plot)
    panel_width <- if (n_terms == 1) 12 else 18
    panel_height <- if (n_terms == 1) 3.9 else max(4.8, 1.5 * n_terms)
    ggsave(plot_file, plot = combined_plot, width = panel_width, height = panel_height, dpi = 300)
    message("Saved selected pathway gene panel: ", plot_file)
  }
}

summary_file <- file.path(PSEUDOBULK_RES_DIR, "deseq2_pseudobulk_summary_by_celltype.csv")
if (!file.exists(summary_file)) {
  stop("Missing pseudobulk summary: ", summary_file)
}

summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
summary_df <- summary_df[
  summary_df$status == "tested" &
    summary_df$n_sig_padj >= MIN_DEGS_PER_CELLTYPE,
  ,
  drop = FALSE
]

cat("Cell types selected for Enrichr:\n")
print(summary_df[, c("cell_type", "n_sig_padj")])

for (i in seq_len(nrow(summary_df))) {
  cell_type_label <- summary_df$cell_type[i]
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", cell_type_label)
  safe_ct <- gsub("^_+|_+$", "", safe_ct)
  deg_file <- file.path(PSEUDOBULK_RES_DIR, paste0("deseq2_HA_vs_other_", safe_ct, "_significant.csv"))

  if (!file.exists(deg_file)) {
    message("Missing DEG file for ", cell_type_label, ": ", deg_file)
    next
  }

  df <- read.csv(deg_file, stringsAsFactors = FALSE)
  df <- df[!is.na(df$padj) & df$padj < PADJ_CUTOFF, , drop = FALSE]

  deg_genes <- unique(df$gene)
  deg_genes <- deg_genes[!is.na(deg_genes) & deg_genes != ""]

  cat("\n", cell_type_label, ": ", length(deg_genes), " significant genes\n", sep = "")

  if (length(deg_genes) == 0) {
    message("No genes for ", cell_type_label, " - skipping enrichment.")
    next
  }

  enrich_res <- enrichr(deg_genes, DBS)

  for (db in names(enrich_res)) {
    table_file <- file.path(
      ENRICHR_TABLE_DIR,
      paste0("enrichr_", safe_ct, "_", db, ".csv")
    )

    write.csv(enrich_res[[db]], table_file, row.names = FALSE)
    message("Saved enrichment table: ", table_file)

    if (RUN_STANDARD_ENRICHR_DOTPLOTS) {
      plot_file <- file.path(
        ENRICHR_FIG_DIR,
        paste0("enrichr_dotplot_", safe_ct, "_", db, ".png")
      )

      plot_enrichr_dotplot(
        enrich_res[[db]],
        title = paste0(cell_type_label, " (", DB_LABELS[db], ")"),
        output_file = plot_file,
        input_gene_count = length(deg_genes),
        top_n = TOP_N_TERMS
      )
    }
  }

  if (RUN_SELECTED_PATHWAY_GENE_BARPLOTS) {
    run_selected_pathway_gene_barplots(
      enrich_res = enrich_res,
      deg_df = df,
      cell_type_label = cell_type_label,
      safe_ct = safe_ct,
      term_patterns = IRON_RELATED_TERM_PATTERNS,
      panel_name = "iron_related",
      gene_label = "Iron-related genes"
    )

    run_selected_pathway_gene_barplots(
      enrich_res = enrich_res,
      deg_df = df,
      cell_type_label = cell_type_label,
      safe_ct = safe_ct,
      term_patterns = TGF_BETA_TERM_PATTERNS,
      panel_name = "TGF_beta_related",
      gene_label = "TGF-beta-related genes"
    )
  }
}

cat("\nEnrichr analysis complete.\n")
cat("Tables : ", ENRICHR_TABLE_DIR, "\n")
cat("Plots  : ", ENRICHR_FIG_DIR, "\n")
