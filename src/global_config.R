##################################
### Project-wide configuration ###
##################################

# project base directory
project_root <- path.expand("~/Thema_R")

# define key directories
data_dir     <- file.path(project_root, "data")
metadata_dir <- file.path(project_root, "metadata")
src_dir      <- file.path(project_root, "src")
results_dir  <- file.path(project_root, "results")
qc_dir       <- file.path(results_dir, "qc")
figures_dir  <- file.path(project_root, "figures")
logs_dir     <- file.path(project_root, "logs")

# create directories if missing
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# data directories for RNA-specific objects
raw_data_dir    <- file.path(data_dir, "raw_counts")
filtered_data_dir   <- file.path(data_dir, "filtered_data")
integrated_data_dir <- file.path(data_dir, "integrated_object")
post_filter_qc_dir  <- file.path(qc_dir, "post_filtering")
filtering_qc_dir    <- file.path(qc_dir, "filtering")
filtering_log_dir   <- file.path(results_dir, "logs")

dir.create(filtered_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(integrated_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(post_filter_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filtering_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filtering_log_dir, recursive = TRUE, showWarnings = FALSE)

###########################
### shared input files ###
###########################

# cell-level filtering settings
scdblfinder_dims <- 20

# gene set files
ferroptosis_geneset_file <- file.path(metadata_dir, "iron_genes", "ferroptosis_genes_curated.xlsx")
iron_uptake_geneset_file <- file.path(metadata_dir, "iron_genes", "iron_uptake_transport_genes.xlsx")
iron_related_geneset_files <- c(ferroptosis_geneset_file, iron_uptake_geneset_file)

# pathway name patterns iron-related used in enrichment plots
iron_related_patterns <- c(
  "ferropt",
  "\\biron\\b",
  "iron ion",
  "iron uptake",
  "ferric",
  "ferrous",
  "ferritin",
  "transferrin",
  "heme",
  "haem",
  "hemoglobin",
  "hepcidin",
  "reactive oxygen species"
)

# pathway name patterns fibrosis-related used in enrichment plots
tgf_beta_patterns <- c(
  "TGF",
  "TGFB",
  "TGF-beta",
  "Transforming Growth Factor",
  "Transforming Growth Factor Beta",
  "SMAD"
)

######################
### atlas settings ###
######################

# integration settings
batch_var <- "sample_id"
n_pcs <- 25
harmony_dims <- 1:25
umap_neighbors <- 50
umap_min_dist <- 0.4
umap_spread <- 1.0
graph_k_param <- 30
cluster_res_grid <- seq(0.2, 1.0, by = 0.1)
selected_cluster_res <- 0.4

# final manual cluster labels
cluster_celltype <- c(
  "1"  = "Sublining fibroblasts (SFRP2+)",
  "2"  = "Lining fibroblasts (PRG4+)",
  "3"  = "Inflammatory fibroblasts (ADAM12+)",
  "4"  = "Activated / cytotoxic T cells",
  "5"  = "Remodeling fibroblasts (HTRA1+)",
  "6"  = "Inflammatory macrophages (IL1B+)",
  "7"  = "Activated endothelial cells",
  "8"  = "Destructive lining fibroblasts (MMP3+)",
  "9"  = "Resident macrophages (C1QC+)",
  "10" = "Pericytes / vascular smooth muscle cells",
  "11" = "Stress-response cells",
  "12" = "B cells",
  "13" = "cDC2",
  "14" = "Mixed-lineage cells"
)

# final cell type colors for annotated UMAPs
cluster_name_colors <- c(
  "Sublining fibroblasts (SFRP2+)" = "#1F77B4",
  "Lining fibroblasts (PRG4+)" = "#D62728",
  "Inflammatory fibroblasts (ADAM12+)" = "#FF7F0E",
  "Activated / cytotoxic T cells" = "#FFBB78",
  "Remodeling fibroblasts (HTRA1+)" = "#2CA02C",
  "Inflammatory macrophages (IL1B+)" = "#98DF8A",
  "Activated endothelial cells" = "#AEC7E8",
  "Destructive lining fibroblasts (MMP3+)" = "#FF9896",
  "Resident macrophages (C1QC+)" = "#9467BD",
  "Pericytes / vascular smooth muscle cells" = "#C5B0D5",
  "Stress-response cells" = "#8C564B",
  "B cells" = "#C49C94",
  "cDC2" = "#9EDAE5",
  "Mixed-lineage cells" = "#F7B6D2"
)

# three canonical markers per final cell type for the annotation dotplot
marker_genes <- list(
  "Sublining fibroblasts (SFRP2+)" = c("SFRP2", "COL1A1", "MFAP5"),
  "Lining fibroblasts (PRG4+)" = c("PRG4", "CLIC5", "DEFB1"),
  "Inflammatory fibroblasts (ADAM12+)" = c("ADAM12", "FGF2", "LIF"),
  "Activated / cytotoxic T cells" = c("GZMA", "BCL11B", "CD3D"),
  "Remodeling fibroblasts (HTRA1+)" = c("HTRA1", "CRTAC1", "NDUFA4L2"),
  "Inflammatory macrophages (IL1B+)" = c("IL1B", "IL1RN", "NLRP3"),
  "Activated endothelial cells" = c("SELE", "EMCN", "CLEC14A"),
  "Destructive lining fibroblasts (MMP3+)" = c("MMP3", "MMP1", "CXCL1"),
  "Resident macrophages (C1QC+)" = c("C1QC", "FOLR2", "TIMD4"),
  "Pericytes / vascular smooth muscle cells" = c("RGS5", "ACTA2", "MYH11"),
  "B cells" = c("IGKC", "CD79A", "IKZF3"),
  "Stress-response cells" = c("HSPA6", "HSPA1A", "DNAJB1"),
  "cDC2" = c("CD1C", "FCER1A", "CLEC10A"),
  "Mixed-lineage cells" = c("RNASE1", "C1QA", "CSN1S1")
)

##########################################
### macrophage annotation information ###
##########################################

# final macrophage subtype labels
macrophage_subcluster_labels <- c(
  "0" = "Inflammatory macrophages (KANK1+)",
  "1" = "Inflammatory macrophages (THBS1+)",
  "2" = "Macrophage-like state (AMTN+)",
  "3" = "Resident macrophages (HSPA6+)",
  "4" = "Red-pulp-like resident macrophages (MERTK+)",
  "5" = "Mixed macrophage-like cells (RNASE1+)",
  "6" = "Plasma-like contaminants",
  "7" = "Low-confidence cells",
  "8" = "Proliferating macrophages"
)

# final macrophage subtype colors
macrophage_subtype_colors <- c(
  "Inflammatory macrophages (KANK1+)" = "#C65A5A",
  "Inflammatory macrophages (THBS1+)" = "#b0e17b",
  "Macrophage-like state (AMTN+)" = "#D98F5C",
  "Resident macrophages (HSPA6+)" = "#8c674b",
  "Red-pulp-like resident macrophages (MERTK+)" = "#4C9F8A",
  "Mixed macrophage-like cells (RNASE1+)" = "#5f90b3",
  "Plasma-like contaminants" = "#B58ACF",
  "Low-confidence cells" = "#9FA4A9",
  "Proliferating macrophages" = "#D95FA7"
)

# canonical markers used to annotate macrophage subtypes
marker_genes_macrophage <- list(
  "Inflammatory macrophages (KANK1+)" = c("KANK1", "KHDRBS3", "BAALC"),
  "Inflammatory macrophages (THBS1+)" = c("THBS1", "IL2RA", "TNIP3"),
  "Macrophage-like state (AMTN+)" = c("AMTN", "SULF1", "ITGBL1"),
  "Resident macrophages (HSPA6+)" = c("HSPA6", "HSPE1-MOB4", "ENSG00000293472"),
  "Red-pulp-like resident macrophages (MERTK+)" = c("MERTK", "CD163", "FCGR3A"),
  "Mixed macrophage-like cells (RNASE1+)" = c("RNASE1", "CRIP1", "HCST"),
  "Plasma-like contaminants" = c("IGHA1", "MZB1", "IGKC"),
  "Low-confidence cells" = c("CCDC200", "RNU5B-1", "ENSG00000285646"),
  "Proliferating macrophages" = c("MKI67", "UBE2C", "CDCA3")
)

# iron and red-pulp-like marker panels
macrophage_iron_features <- c(
  "SPIC", "HMOX1", "SLC40A1", "TFRC", "STEAP3", "FTH1", "FTL",
  "SLC11A2", "CP", "NCOA4", "GPX4", "ACSL4", "AIFM2",
  "NQO1", "GCLC", "GCLM", "ALOX5", "ALOX15", "SAT1"
)

red_pulp_like_features <- c(
  "CD163", "CD14", "FCGR2A", "FCGR3A", "FCGR2B",
  "CYBB", "SLC48A1", "HMOX1", "SLC40A1"
)
