setwd("~/Thema_R")

gene <- "NOTCH3"
padj_cutoff <- 0.05

res_dir <- "results/pseudobulk_deseq2"
files <- list.files(
  res_dir,
  pattern = "^deseq2_HA_vs_other_.*\\.csv$",
  full.names = TRUE
)

files <- files[!grepl("_significant\\.csv$", files)]

found <- FALSE

for (f in files) {
  tab <- read.csv(f, stringsAsFactors = FALSE)

  cell_type <- basename(f)
  cell_type <- sub("^deseq2_HA_vs_other_", "", cell_type)
  cell_type <- sub("\\.csv$", "", cell_type)
  cell_type <- gsub("_", " ", cell_type)

  row <- tab[tab$gene == gene, ]

  if (nrow(row) == 0) {
    next
  }

  found <- TRUE

  logfc <- row$log2FoldChange[1]
  padj <- row$padj[1]

  if (!is.na(padj) && padj < padj_cutoff) {
    direction <- ifelse(logfc > 0, "up in HA", "up in other")

    cat(
      gene, "is differentially expressed in", cell_type,
      "| log2FC =", round(logfc, 3),
      "| padj =", signif(padj, 3),
      "|", direction, "\n"
    )
  } else {
    cat(
      gene, "is NOT significant in", cell_type,
      "| log2FC =", round(logfc, 3),
      "| padj =", signif(padj, 3), "\n"
    )
  }
}

if (!found) {
  cat(gene, "was not found in the tested DESeq2 result tables.\n")
}