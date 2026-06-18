# ===============================================================
#   15c_GEX_fgsea_macrophage_subtypes_HA_vs_other.R
# ===============================================================

# Run ranked GSEA with fgsea on macrophage-subtype pseudobulk DESeq2 results.
# All tested genes are ranked by the DESeq2 Wald statistic when available.

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
})

setwd("/data/home/alessandro.digilio/Thema_R")
source("src/global_config.R")

# Input and output folders
PSEUDOBULK_RES_DIR <- file.path(RESULTS_DIR, "pseudobulk_deseq2_macrophage_subtypes")
FGSEA_RES_DIR <- file.path(RESULTS_DIR, "fgsea_pseudobulk_macrophage_subtypes_HA_vs_other")
FGSEA_TABLE_DIR <- file.path(FGSEA_RES_DIR, "significant_tables")
FGSEA_IRON_TABLE_DIR <- file.path(FGSEA_RES_DIR, "iron_related_significant_tables")
FGSEA_FIG_DIR <- file.path(FIGURES_DIR, "fgsea_pseudobulk_macrophage_subtypes_HA_vs_other")
FGSEA_CURVE_DIR <- file.path(FGSEA_FIG_DIR, "gsea_curves")

dir.create(FGSEA_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FGSEA_IRON_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FGSEA_CURVE_DIR, recursive = TRUE, showWarnings = FALSE)

# Set parameters
PADJ_THR <- 0.05
MIN_SIZE <- 5
MAX_SIZE <- 500
FGSEA_SEED <- 1234
USE_MULTILEVEL <- TRUE
N_PERMUTATIONS <- 10000

# MSigDB collections used for ranked GSEA
MSIGDB_SPECIES <- "Homo sapiens"
MSIGDB_COLLECTIONS <- list(
  "GO_BP" = c("GO:BP"),
  "KEGG" = c("CP:KEGG", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"),
  "Reactome" = c("CP:REACTOME")
)

# Only significant pathways matching these terms are plotted as GSEA curves.
IRON_RELATED_PATTERNS <- c(
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
  "hepcidin"
)

sanitize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# Load one MSigDB subcollection, supporting both old and new msigdbr arguments.
load_msigdb_subcollection <- function(subcollection) {
  out <- tryCatch(
    msigdbr::msigdbr(
      species = MSIGDB_SPECIES,
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
      species = MSIGDB_SPECIES,
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

# GO BP lives in collection C5, so it needs a small separate loader.
load_msigdb_gobp <- function() {
  out <- tryCatch(
    msigdbr::msigdbr(
      species = MSIGDB_SPECIES,
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
      species = MSIGDB_SPECIES,
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

# Build GO BP, KEGG and Reactome pathways from msigdbr.
load_msigdb_pathways <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "Package 'msigdbr' is required for GO BP / KEGG / Reactome fgsea. ",
      "Install it in scatac_gex_env and rerun this script."
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
    for (subcollection in MSIGDB_COLLECTIONS[[db]]) {
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

# Flag pathway names that are directly related to iron, heme or ferroptosis.
is_iron_related_pathway <- function(pathway_name) {
  pathway_text <- tolower(pathway_name)
  pattern <- paste(IRON_RELATED_PATTERNS, collapse = "|")
  grepl(pattern, pathway_text, perl = TRUE)
}

# Build the ranked vector used by fgsea.
build_rank_vector <- function(deg_df) {
  if ("stat" %in% colnames(deg_df)) {
    rank_df <- deg_df[!is.na(deg_df$stat), c("gene", "stat"), drop = FALSE]
    colnames(rank_df) <- c("gene", "rank_value")
  } else {
    rank_df <- deg_df[
      !is.na(deg_df$log2FoldChange) &
        !is.na(deg_df$pvalue) &
        deg_df$pvalue > 0,
      c("gene", "log2FoldChange", "pvalue"),
      drop = FALSE
    ]
    rank_df$rank_value <- sign(rank_df$log2FoldChange) * -log10(rank_df$pvalue)
    rank_df <- rank_df[, c("gene", "rank_value"), drop = FALSE]
  }

  rank_df$gene <- toupper(as.character(rank_df$gene))
  rank_df <- rank_df[!is.na(rank_df$gene) & rank_df$gene != "", , drop = FALSE]
  rank_df <- rank_df[!is.na(rank_df$rank_value), , drop = FALSE]

  # fgsea needs one value per gene, so keep the strongest duplicate if present.
  # This absolute ordering is only for duplicate removal; the final ranking stays signed.
  rank_df <- rank_df[order(abs(rank_df$rank_value), decreasing = TRUE), , drop = FALSE]
  rank_df <- rank_df[!duplicated(rank_df$gene), , drop = FALSE]

  ranks <- rank_df$rank_value
  names(ranks) <- rank_df$gene

  # GSEA should use all ranked genes, not only significant DEGs.
  # Positive ranks are enriched toward HA, negative ranks toward other.
  sort(ranks, decreasing = TRUE)
}

# Load GO BP, KEGG and Reactome pathways.
pathway_data <- load_msigdb_pathways()
pathways <- pathway_data$pathways
pathway_info <- pathway_data$pathway_info

cat("Pathways loaded:\n")
print(table(pathway_info$database))

# Check macrophage-subtype pseudobulk summary.
summary_file <- file.path(PSEUDOBULK_RES_DIR, "deseq2_pseudobulk_summary_by_macrophage_subtype.csv")
if (!file.exists(summary_file)) {
  stop("Missing pseudobulk summary file: ", summary_file)
}

summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
summary_df <- summary_df[summary_df$status == "tested", , drop = FALSE]

all_sig_results <- list()
all_iron_sig_results <- list()

for (i in seq_len(nrow(summary_df))) {
  subtype_here <- summary_df$macrophage_subtype[i]
  safe_subtype <- sanitize_label(subtype_here)
  deg_file <- file.path(PSEUDOBULK_RES_DIR, paste0("deseq2_HA_vs_other_", safe_subtype, ".csv"))

  if (!file.exists(deg_file)) {
    cat("Missing DESeq2 result for ", subtype_here, ": ", deg_file, "\n", sep = "")
    next
  }

  cat("\n============================================================\n")
  cat("Running fgsea for macrophage subtype:", subtype_here, "\n")

  deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)
  ranks <- build_rank_vector(deg_df)

  cat("Ranked genes:", length(ranks), "\n")
  if (length(ranks) < MIN_SIZE) {
    cat("Skipping: too few ranked genes.\n")
    next
  }

  set.seed(FGSEA_SEED)

  if (USE_MULTILEVEL) {
    fgsea_res <- fgseaMultilevel(
      pathways = pathways,
      stats = ranks,
      minSize = MIN_SIZE,
      maxSize = MAX_SIZE,
      eps = 0,
      nproc = 1
    )
  } else {
    fgsea_res <- fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = MIN_SIZE,
      maxSize = MAX_SIZE,
      nperm = N_PERMUTATIONS,
      nproc = 1
    )
  }

  fgsea_res <- as.data.frame(fgsea_res)
  fgsea_res <- fgsea_res[order(fgsea_res$padj, fgsea_res$pval), , drop = FALSE]
  fgsea_res$macrophage_subtype <- subtype_here
  fgsea_res$database <- as.character(pathway_info[fgsea_res$pathway, "database", drop = TRUE])
  fgsea_res$pathway_name <- as.character(pathway_info[fgsea_res$pathway, "gs_name", drop = TRUE])
  fgsea_res$is_iron_related <- is_iron_related_pathway(fgsea_res$pathway_name)
  fgsea_res$leadingEdge <- vapply(fgsea_res$leadingEdge, paste, character(1), collapse = ";")
  fgsea_res <- fgsea_res[
    ,
    c(
      "macrophage_subtype", "database", "pathway_name", "is_iron_related",
      "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge"
    ),
    drop = FALSE
  ]

  sig_res <- fgsea_res[!is.na(fgsea_res$padj) & fgsea_res$padj < PADJ_THR, , drop = FALSE]
  all_sig_results[[subtype_here]] <- sig_res

  cat("Significant pathways at FDR < ", PADJ_THR, ": ", nrow(sig_res), "\n", sep = "")
  print(sig_res[, c("database", "pathway_name", "NES", "pval", "padj", "leadingEdge"), drop = FALSE])

  for (db in names(MSIGDB_COLLECTIONS)) {
    sig_db <- sig_res[sig_res$database == db, , drop = FALSE]
    if (nrow(sig_db) == 0) next

    write.csv(
      sig_db,
      file.path(FGSEA_TABLE_DIR, paste0("fgsea_HA_vs_other_", safe_subtype, "_", db, "_significant.csv")),
      row.names = FALSE
    )
  }

  iron_sig_res <- sig_res[sig_res$is_iron_related, , drop = FALSE]
  all_iron_sig_results[[subtype_here]] <- iron_sig_res

  if (nrow(iron_sig_res) > 0) {
    write.csv(
      iron_sig_res,
      file.path(FGSEA_IRON_TABLE_DIR, paste0("fgsea_HA_vs_other_", safe_subtype, "_iron_related_significant.csv")),
      row.names = FALSE
    )
  }

  # Save enrichment curves only for significant iron-related pathways.
  cat("Significant iron-related pathways: ", nrow(iron_sig_res), "\n", sep = "")
  if (nrow(iron_sig_res) == 0) {
    cat("No iron-related GSEA curves saved for ", subtype_here, ".\n", sep = "")
    next
  }

  for (j in seq_len(nrow(iron_sig_res))) {
    pathway_here <- iron_sig_res$pathway[j]
    pathway_name_here <- iron_sig_res$pathway_name[j]
    pathway_safe <- paste0(
      sanitize_label(iron_sig_res$database[j]),
      "_",
      sprintf("%03d", j),
      "_",
      substr(sanitize_label(pathway_name_here), 1, 80)
    )

    p_curve <- plotEnrichment(pathways[[pathway_here]], ranks) +
      labs(
        title = pathway_name_here,
        subtitle = paste0(
          subtype_here,
          " | ", iron_sig_res$database[j],
          " | NES = ", round(iron_sig_res$NES[j], 2),
          " | padj = ", signif(iron_sig_res$padj[j], 3)
        ),
        x = "Ranked genes",
        y = "Enrichment score"
      ) +
      theme_classic(base_size = 16) +
      theme(
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
        panel.grid = element_blank(),
        plot.margin = margin(18, 18, 18, 18)
      )

    ggsave(
      filename = file.path(FGSEA_CURVE_DIR, paste0("gsea_curve_", safe_subtype, "_", pathway_safe, ".png")),
      plot = p_curve,
      width = 8,
      height = 6,
      dpi = 600
    )
  }
}

if (length(all_sig_results) > 0) {
  all_sig_df <- do.call(rbind, all_sig_results)
  rownames(all_sig_df) <- NULL
  write.csv(
    all_sig_df,
    file.path(FGSEA_RES_DIR, "fgsea_HA_vs_other_all_macrophage_subtypes_significant.csv"),
    row.names = FALSE
  )
}

if (length(all_iron_sig_results) > 0) {
  all_iron_sig_df <- do.call(rbind, all_iron_sig_results)
  rownames(all_iron_sig_df) <- NULL
  write.csv(
    all_iron_sig_df,
    file.path(FGSEA_RES_DIR, "fgsea_HA_vs_other_all_macrophage_subtypes_iron_related_significant.csv"),
    row.names = FALSE
  )
}

cat("\n============================================================\n")
cat("Macrophage-subtype fgsea ranked GSEA complete.\n")
cat("Pathways    : MSigDB GO BP, KEGG and Reactome via msigdbr.\n")
cat("Ranking     : DESeq2 Wald statistic; fallback = sign(log2FC) * -log10(pvalue).\n")
cat("Gene input  : all ranked genes from the full DESeq2 table, not only significant DEGs.\n")
cat("Curves      : only significant iron-related pathways.\n")
cat("Reproducible: set.seed(", FGSEA_SEED, "), nproc = 1.\n", sep = "")
cat("Sig tables  : ", FGSEA_TABLE_DIR, "\n")
cat("Iron tables : ", FGSEA_IRON_TABLE_DIR, "\n")
cat("Curves      : ", FGSEA_CURVE_DIR, "\n")
cat("============================================================\n")
