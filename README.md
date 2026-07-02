# THEMA: Hemophilic Arthropathy Synovial scRNA-seq

### Single-cell RNA-seq analysis of synovial membrane samples from hemophilic arthropathy (HA), osteoarthritis (OA) and rheumatoid arthritis (RA)

Hemophilic arthropathy is associated with recurrent joint bleeding, synovial iron deposition, low-grade inflammation and progressive tissue remodeling. The single-cell analysis defines the synovial cellular landscape in HA and compares it with OA/RA tissues, with focus on macrophages, fibroblast-like synoviocytes and endothelial cells.

The workflow starts from an integrated synovial atlas and focus on selected macrophage, endothelial and fibroblast compartments.

## Repository Organization

```text
src/
  global_config.R
  atlas/
  macrophage_subclusters/
  endothelial_subclusters/
  fibroblast_subcluster/
  trajectory_analysis/
  pre_processing/
  pbs/

metadata/
  filters.xlsx
  subject_info.xlsx
  iron_genes/
  samplesheets/

env/
  thema_environment.yml
```

- `global_config.R` defines paths, sample-independent settings, final annotations, marker panels and plotting palettes
- `utils.R` scripts contain helper functions used by the scripts in that folder

## Reproducibility

### Environment

The R environment (R libraries) used for the analysis is described in:

```text
env/thema_environment.yml
```

To reproduce it:

```bash
conda env create -f env/thema_environment.yml
conda activate scatac_gex_env
```

## Workflow

### FASTQ Processing

FASTQ files were processed with `nf-core/scrnaseq` using Cell Ranger and a GRCh38/GENCODE v46 reference. The PBS launcher and samplesheets used for the count-matrix generation are stored in:

```text
src/pre_processing/nfcore_scrnaseq_cellranger_all_samples.pbs
metadata/samplesheets/
```

### Atlas-Level Analysis

This block inspects disease-associated transcriptional programs across the atlas.

```text
src/atlas/00_qc_prefiltering.R
src/atlas/01_filter_samples.R
src/atlas/02_integration.R
src/atlas/03_annotation.R
src/atlas/04_celltype_composition_ha_vs_other.R
src/atlas/05_pseudobulk_deseq2_ha_vs_other.R
src/atlas/06_enrichment_ha_vs_other.R
src/atlas/07_ferroptosis_ha_vs_other.R
```

### Focused Subcluster Analyses

Macrophages:

```text
src/macrophage_subclusters/00_subclustering.R
src/macrophage_subclusters/01_annotation.R
src/macrophage_subclusters/02_subtype_composition_ha_vs_other.R
src/macrophage_subclusters/03_pseudobulk_deseq2_ha_vs_other.R
src/macrophage_subclusters/04_gsea_ha_vs_other.R
src/macrophage_subclusters/05_ferroptosis_ha_vs_other.R
```

Endothelial cells:

```text
src/endothelial_subclusters/00_subclustering.R
src/endothelial_subclusters/01_annotation.R
src/endothelial_subclusters/02_subtype_composition_ha_vs_other.R
src/endothelial_subclusters/03_pseudobulk_deseq2_ha_vs_other.R
src/endothelial_subclusters/04_gsea_ha_vs_other.R
src/endothelial_subclusters/05_ferroptosis_ha_vs_other.R
```

Destructive lining fibroblasts:

```text
src/fibroblast_subcluster/00_subclustering.R
src/fibroblast_subcluster/01_annotation.R
src/fibroblast_subcluster/02_de_ha_enriched_states.R
src/fibroblast_subcluster/03_functional_analysis_ha_enriched_states.R
```

These analyses resolve internal heterogeneity within the main compartments. The macrophage workflow characterizes resident, inflammatory and red-pulp-like states; the endothelial workflow focuses on activated and stress-response vascular states; the fibroblast workflow tests the HA-enriched destructive lining fibroblast states and their inflammatory, matrix-remodeling and iron-related programs.

### Trajectory Analysis

```text
src/trajectory_analysis/monocle3_trajectory_analysis.R
```

Monocle3 is applied to the already reclustered macrophage, endothelial and destructive lining fibroblast compartments

## Local Outputs

Large and generated files are kept outside version control:

```text
data/
results/
  atlas/
  macrophages/
  endothelial/
  fibroblasts/
  trajectory_analysis/
figures/
  atlas/
  macrophages/
  endothelial/
  fibroblasts/
  trajectory_analysis/
  paper/
logs/
```
