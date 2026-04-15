# ---------------------------------------
# 10_GEX_iron_genes_in_pseudobulk_DEGs_HA_vs_other.R
# Iron-related genes among pseudobulk DEGs
# ---------------------------------------

# Collect iron-related / ferroptosis-related genes from the metadata gene
# lists, intersect them with the pseudobulk DESeq2 results, print the hits

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
})

sanitize_celltype_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+", "", x)
  x <- gsub("_+$", "", x)
  x
}

# Work from the project root
setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Create output directories used by this step
PSEUDOBULK_RES_DIR <- file.path(RESULTS_DIR, "pseudobulk_deseq2")
FERROPTOSIS_RES_DIR <- file.path(RESULTS_DIR, "ferroptosis")
FERROPTOSIS_FIG_DIR <- file.path(FIGURES_DIR, "ferroptosis")
dir.create(FERROPTOSIS_RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FERROPTOSIS_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_XLSX <- file.path(
  FERROPTOSIS_RES_DIR,
  "iron_related_genes_in_pseudobulk_DEGs_HA_vs_other.xlsx"
)
OUTPUT_CSV <- file.path(
  FERROPTOSIS_RES_DIR,
  "iron_related_genes_in_pseudobulk_DEGs_HA_vs_other.csv"
)
OUTPUT_BUBBLE_PLOT <- file.path(
  FERROPTOSIS_FIG_DIR,
  "iron_related_genes_in_pseudobulk_DEGs_HA_vs_other_bubble_heatmap.png"
)

IRON_GENE_FILES <- c(
  file.path(METADATA_DIR, "iron_genes", "ferroptosis_genes_curated.xlsx"),
  file.path(METADATA_DIR, "iron_genes", "iron_uptake_transport_genes.xlsx")
)

# Set significance threshold
PADJ_THR <- 0.05

# Check required input files
summary_file <- file.path(PSEUDOBULK_RES_DIR, "deseq2_pseudobulk_summary_by_celltype.csv")
if (!file.exists(summary_file)) {
  stop("Missing pseudobulk summary file: ", summary_file)
}

for (f in IRON_GENE_FILES) {
  if (!file.exists(f)) {
    stop("Missing iron-related gene file: ", f)
  }
}

# Read and combine the iron-related gene lists
iron_genes <- character(0)

for (f in IRON_GENE_FILES) {
  x <- openxlsx::read.xlsx(f)
  if (ncol(x) == 0) next

  first_col <- x[[1]]
  genes_here <- c(colnames(x)[1], first_col)
  genes_here <- toupper(trimws(as.character(genes_here)))
  genes_here <- genes_here[!is.na(genes_here) & genes_here != ""]
  iron_genes <- c(iron_genes, genes_here)
}

iron_genes <- sort(unique(iron_genes))

cat("Iron-related genes loaded:", length(iron_genes), "\n")
cat(paste(iron_genes, collapse = ", "), "\n\n")

# Keep only tested cell types from the pseudobulk summary
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
summary_df <- summary_df[summary_df$status == "tested", , drop = FALSE]

all_hits <- list()
summary_rows <- list()

# Scan each cell type result and keep only iron-related genes
for (i in seq_len(nrow(summary_df))) {
  cell_type_here <- summary_df$cell_type[i]
  safe_ct <- sanitize_celltype_name(cell_type_here)
  deg_file <- file.path(PSEUDOBULK_RES_DIR, paste0("deseq2_HA_vs_other_", safe_ct, ".csv"))

  if (!file.exists(deg_file)) {
    cat("Missing result file for", cell_type_here, "\n")
    next
  }

  df <- read.csv(deg_file, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(df)) next

  df$gene_upper <- toupper(as.character(df$gene))
  df <- df[df$gene_upper %in% iron_genes, , drop = FALSE]

  if (nrow(df) == 0) {
    cat("\n", cell_type_here, ": no iron-related genes found among tested genes\n", sep = "")
    next
  }

  df$cell_type <- cell_type_here
  df$direction_in_HA <- ifelse(
    is.na(df$log2FoldChange),
    NA_character_,
    ifelse(df$log2FoldChange > 0, "Up in HA",
           ifelse(df$log2FoldChange < 0, "Down in HA", "No change"))
  )
  df$is_significant <- !is.na(df$padj) & df$padj < PADJ_THR

  df_out <- df[, c("cell_type", "gene", "log2FoldChange", "padj", "direction_in_HA", "is_significant")]
  df_out <- df_out[order(df_out$padj, -abs(df_out$log2FoldChange)), , drop = FALSE]

  all_hits[[cell_type_here]] <- df_out
  summary_rows[[cell_type_here]] <- data.frame(
    cell_type = cell_type_here,
    n_iron_genes_found = nrow(df_out),
    n_significant_padj_lt_0.05 = sum(df_out$is_significant, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  cat("\n============================================================\n")
  cat(cell_type_here, "\n")
  cat("Iron-related genes found:", nrow(df_out), "\n")
  cat("Significant at padj <", PADJ_THR, ":", sum(df_out$is_significant, na.rm = TRUE), "\n")
  print(df_out)
}

# Stop if no iron-related genes were found in any tested cell type
if (length(all_hits) == 0) {
  stop("No iron-related genes found in the pseudobulk DESeq2 results.")
}

# Save the full iron-related hit table as CSV
all_hits_df <- do.call(rbind, all_hits)
rownames(all_hits_df) <- NULL
write.csv(all_hits_df, OUTPUT_CSV, row.names = FALSE)
summary_df_out <- do.call(rbind, summary_rows)
rownames(summary_df_out) <- NULL

# Save only significant hits
wb <- createWorkbook()

sig_hits_df <- all_hits_df[all_hits_df$is_significant, , drop = FALSE]
addWorksheet(wb, "significant_padj_lt_0.05")
writeData(wb, "significant_padj_lt_0.05", sig_hits_df)

addWorksheet(wb, "summary")
writeData(wb, "summary", summary_df_out)

for (nm in names(all_hits)) {
  df_here <- all_hits[[nm]]
  df_here <- df_here[df_here$is_significant, , drop = FALSE]
  if (nrow(df_here) == 0) next

  sheet_name <- sanitize_celltype_name(nm)
  sheet_name <- substr(sheet_name, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df_here)
}

# Write the final workbook
saveWorkbook(wb, OUTPUT_XLSX, overwrite = TRUE)

# Plot a bubble heatmap of significant iron-related DEGs
if (nrow(sig_hits_df) > 0) {
  plot_df <- sig_hits_df
  plot_df$neg_log10_padj <- -log10(plot_df$padj)

  celltype_order <- unique(summary_df_out$cell_type)
  celltype_order <- celltype_order[celltype_order %in% unique(plot_df$cell_type)]
  plot_df$cell_type <- factor(plot_df$cell_type, levels = celltype_order)

  gene_order_df <- aggregate(
    abs(log2FoldChange) ~ gene,
    data = plot_df,
    FUN = max
  )
  gene_order_df <- gene_order_df[order(gene_order_df$`abs(log2FoldChange)`, decreasing = TRUE), , drop = FALSE]
  plot_df$gene <- factor(plot_df$gene, levels = rev(gene_order_df$gene))

  max_abs_lfc <- max(abs(plot_df$log2FoldChange), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) {
    max_abs_lfc <- 1
  }

  p_bubble <- ggplot(plot_df, aes(x = cell_type, y = gene)) +
    geom_point(aes(size = neg_log10_padj, fill = log2FoldChange), shape = 21, color = "black", stroke = 0.25) +
    scale_fill_gradient2(
      low = "#4C78A8",
      mid = "white",
      high = "#7A1F2B",
      midpoint = 0,
      limits = c(-max_abs_lfc, max_abs_lfc),
      name = "log2FC"
    ) +
    scale_size_continuous(
      name = "-log10 adj p",
      range = c(2.5, 9)
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.position = "top",
      legend.key.width = grid::unit(0.45, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(2.6, "cm"),
        barheight = grid::unit(0.35, "cm")
      ),
      size = guide_legend(title.position = "top")
    )

  ggsave(
    filename = OUTPUT_BUBBLE_PLOT,
    plot = p_bubble,
    width = 12,
    height = max(6, 0.28 * length(levels(plot_df$gene)) + 2),
    dpi = 600
  )
}

cat("\n============================================================\n")
cat("Iron-related DEG summary complete.\n")
cat("\nSummary by cell type:\n")
print(summary_df_out)
cat("\nTotal iron-related genes found:", nrow(all_hits_df), "\n")
cat("Total significant iron-related genes:", nrow(sig_hits_df), "\n")
cat("CSV  : ", OUTPUT_CSV, "\n", sep = "")
cat("XLSX : ", OUTPUT_XLSX, "\n", sep = "")
if (nrow(sig_hits_df) > 0) {
  cat("PLOT : ", OUTPUT_BUBBLE_PLOT, "\n", sep = "")
}
cat("============================================================\n")
