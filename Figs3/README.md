# Figs3: Somatic-cell QC and downstream integration

## Overview

This directory contains the analysis workflow used for **Supplementary Fig. 3 (Figs3)**, focusing on somatic-cell quality control, integration, clustering, and cell-type annotation.

The workflow is organized into two sequential modules:

1. **Somatic QC**  
   Performs per-sample quality control on Seurat objects, including mitochondrial-content filtering, UMI outlier filtering, and doublet detection with `scDblFinder`.

2. **Somatic processing**  
   Uses the QC-processed object for two rounds of SCTransform/Harmony-based integration, cluster removal, dimensionality reduction, and cell-type annotation.

In the revised workflow, **multiple input `.rds` files can be provided to the QC module**, each with an explicit `sample` label. Each input object is QC-processed independently, then all QC-passed objects are merged into a single Seurat object. The merged object retains the user-specified sample label in `meta.data$sample`, which is then used in downstream processing.

---

## Repository contents

- `figs3_somatic_qc.txt`  
  R script for somatic-cell quality control.

- `figs3_somatic_processing.txt`  
  R script for downstream integration and annotation.

- `run_figs3_somatic_processing_combined.sh`  
  Combined runner for QC only, downstream processing only, or the full workflow.

- `figs3_combined.yaml`  
  YAML configuration file for both modules.

> In this repository, some R scripts may carry a `.txt` suffix for packaging consistency. They should still be executed with `Rscript`.

---

## Requirements

### R packages

The workflow requires at least the following R packages:

- `Seurat`
- `sctransform`
- `harmony`
- `yaml`
- `dplyr`
- `ggplot2`
- `Matrix`
- `SingleCellExperiment`
- `scDblFinder`

Depending on the exact local environment, additional visualization packages may also be needed.

### Input format

Each input file must be a valid **Seurat object** saved as `.rds`.

For the multi-input QC mode, each `.rds` should correspond to one sample/replicate, and the sample label is supplied explicitly at run time or in the YAML file.

---

## Workflow logic

### 1. Somatic QC

For each input Seurat object:

- compute mitochondrial fraction (`percent.mt`)
- filter cells above the configured mitochondrial threshold
- filter high-UMI outliers using an IQR-based rule
- run `scDblFinder`
- retain singlets only
- assign the user-provided `sample` label to all retained cells

After all inputs are QC-processed independently, the resulting Seurat objects are merged into a single QC-filtered object.

### 2. Somatic processing

The merged QC-filtered object is then processed through:

- sample-aware preprocessing
- feature selection excluding ribosomal/mitochondrial genes
- cell-cycle scoring
- **Round 1** Harmony integration, clustering, and UMAP
- removal of selected clusters after Round 1
- **Round 2** Harmony integration, clustering, and UMAP
- removal of selected clusters after Round 2
- annotation of final clusters using `celltype_map`

If the input object already contains a usable `sample` column, the downstream module keeps it. A fallback `sample_map` can still be used for legacy inputs that do not yet contain `sample`.

---

## Recommended usage

### A. Full workflow: multiple raw `.rds` inputs -> QC -> merge -> downstream processing

```bash
bash Figs3/run_figs3_somatic_processing_combined.sh \
  -o result/Figs3 \
  --set figs3_qc.qc_inputs.0.path=/path/PD14_1.rds \
  --set figs3_qc.qc_inputs.0.sample=PD14-1 \
  --set figs3_qc.qc_inputs.1.path=/path/PD14_2.rds \
  --set figs3_qc.qc_inputs.1.sample=PD14-2 \
  --set figs3_qc.qc_inputs.2.path=/path/PD49_1.rds \
  --set figs3_qc.qc_inputs.2.sample=PD49-1 \
  --set figs3_qc.qc_inputs.3.path=/path/PD49_2.rds \
  --set figs3_qc.qc_inputs.3.sample=PD49-2
```

This mode:

- QC-processes each input independently
- merges all QC-passed objects
- stores the provided labels in `meta.data$sample`
- runs downstream integration and annotation on the merged object

### B. QC only

```bash
bash Figs3/run_figs3_somatic_processing_combined.sh \
  --only qc \
  -o result/Figs3 \
  --set figs3_qc.qc_inputs.0.path=/path/PD14_1.rds \
  --set figs3_qc.qc_inputs.0.sample=PD14-1 \
  --set figs3_qc.qc_inputs.1.path=/path/PD14_2.rds \
  --set figs3_qc.qc_inputs.1.sample=PD14-2
```

### C. Downstream processing only

Use this mode when a merged QC-processed Seurat object is already available:

```bash
bash Figs3/run_figs3_somatic_processing_combined.sh \
  --only somatic_processing \
  -o result/Figs3 \
  --set figs3_somatic_processing.input_rds=/path/obj_merge_qc.rds
```

---

## Configuration

### QC module

The QC module accepts either a legacy single-input field:

```yaml
qc_input_rds: data/Figs3/qc/somatic_seurat_object.rds
```

or the preferred multi-input structure:

```yaml
qc_inputs:
  - path: data/Figs3/qc/PD14_1.rds
    sample: PD14-1
  - path: data/Figs3/qc/PD14_2.rds
    sample: PD14-2
  - path: data/Figs3/qc/PD49_1.rds
    sample: PD49-1
```

Important QC parameters include:

- `qc_percent_mt_max`
- `qc_umi_iqr_multiplier`
- `qc_doublet_rate`
- `qc_species`

### Downstream module

Important downstream parameters include:

- `input_rds`
- `exclude_regex`
- `nfeatures_round1`
- `nfeatures_round2`
- `pipeline_round1.*`
- `pipeline_round2.*`
- `exclude_clusters_round1`
- `exclude_clusters_round2`
- `celltype_map`

---

## Expected outputs

### QC outputs

Under `result/Figs3/qc/` (or the configured output root):

- `somatic.qc_processed.rds`
- QC violin plots at different filtering stages
- `somatic.doublet_scores.tsv`
- `somatic.qc_stats.yaml`
- parameter snapshot and `sessionInfo`

### Downstream outputs

Under `result/Figs3/somatic_processing/`:

- `somatic.round1.rds`
- `somatic.final.rds`
- UMAP plots for intermediate/final objects
- cell count summary tables
- parameter snapshot and `sessionInfo`
- processing log

---

## Notes for reproducibility

- Always prefer absolute paths for input `.rds` files in command-line overrides.
- For multi-sample analyses, provide unique and stable `sample` labels at the QC stage; these labels are propagated into the merged object.
- If running downstream processing only, ensure that the input object already contains a valid `sample` column, or provide a fallback `sample_map` in the YAML configuration.
- Cluster exclusion and cell-type annotation are dataset-specific and may require adjustment for independent datasets.

---

## Citation

If you use this code or adapt this workflow in your own work, please cite the associated study and the software packages used in the analysis, including Seurat, Harmony, and scDblFinder.
