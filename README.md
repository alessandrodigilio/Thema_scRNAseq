# THEMA single-cell RNA-seq analysis

This repository contains the R and PBS scripts used for the single-cell RNA-seq analysis of synovial membrane samples from patients with hemophilic arthropathy and comparator osteoarthritis/rheumatoid arthritis samples.

The analysis focuses on the cellular and transcriptional architecture of hemophilic arthropathy synovium, with particular attention to iron metabolism, ferroptosis-related programs, macrophage states, fibroblast states, endothelial states, and trajectory analysis.

## Repository Structure

- `src/`: R scripts for quality control, filtering, integration, annotation, differential expression, enrichment, subclustering, ferroptosis scoring, and trajectory analysis.
- `src/pbs/`: PBS launcher scripts for the R analysis steps.
- `src/pre_processing/`: upstream nf-core/scrnaseq Cell Ranger launcher and raw-data utility scripts.
- `metadata/`: sample metadata, filtering thresholds, curated gene sets, and samplesheet templates.
- `env/`: conda environment file used for the R analysis.

Raw data, intermediate Seurat objects, results, logs, and generated figures are not tracked in this repository.

## Data Layout

The scripts expect the following local folders, which are ignored by git:

```text
data/raw_counts/
data/filtered_data/
data/integrated_object/
results/
figures/
logs/
```

Before running the workflow on a new system, update local paths in `src/global_config.R` and adapt the FASTQ paths in `metadata/samplesheets/scRNA_samplesheet_all_samples.csv`.

## Environment

The R environment can be recreated from:

```bash
conda env create -f env/environment.yml
conda activate thema_r_env
```

The upstream FASTQ-to-count-matrix processing was performed with nf-core/scrnaseq using Cell Ranger. The corresponding PBS launcher is provided in:

```text
src/pre_processing/nfcore_scrnaseq_cellranger_all_samples.pbs
```

## Analysis Workflow

Run the scripts in numerical order:

```text
00_qc_prefiltering.R
01_filter_samples.R
02_integration.R
04_annotation.R
05_celltype_composition_ha_vs_other.R
06_pseudobulk_deseq2_ha_vs_other.R
07_enrichment_ha_vs_other.R
08_ferroptosis_ha_vs_other.R: ferroptosis scoring and iron-related DEG summary
10-14: macrophage subclustering and downstream analyses
15-19: endothelial subclustering and downstream analyses
20-23: destructive lining fibroblast subclustering and downstream analyses
24_quick_fibroblast_tgfbeta_featureplots.R
25_monocle3_trajectory_analysis.R
26_paper_subcluster_umaps.R
```

PBS scripts matching the main analysis steps are available in `src/pbs/`.

## Notes

This repository is intended to store reproducible analysis code and lightweight metadata only. Large generated files should remain outside version control.
