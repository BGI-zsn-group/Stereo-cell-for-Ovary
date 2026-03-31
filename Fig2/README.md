# Fig. 2 | Cell cropping and two-pass SCTransform + Harmony integration

## Overview

This directory contains the code used to reproduce the analyses presented in **Fig. 2**. The workflow consists of two conceptually distinct components:

1. **QuPath-based cell cropping**, used to export standardized image crops centered on annotated cells for image-level inspection and documentation.
2. **Transcriptomic integration**, which combines multiple per-sample Seurat objects using a **two-pass SCTransform + Harmony** strategy and generates a single integrated Seurat object for downstream analyses and figure generation.

The image-export step and the transcriptomic integration step are complementary but independent. The former operates on microscopy annotations in QuPath, whereas the latter operates on per-sample Seurat objects in R.

Detailed step-by-step operating procedures, including how the QuPath annotations were created and how the image crops were used in the study, are described in the **article protocol**. This repository README focuses on the code structure, required inputs, and reproducible execution of the computational workflow.

---

## Repository contents

- `cell_cropping.groovy`  
  QuPath Groovy script for exporting fixed-size TIFF crops centered on annotated cells.

- `run_fig2_combined.sh`  
  Recommended entry point for the integration workflow. This wrapper extracts the Fig. 2 module from the combined YAML configuration and runs the Harmony-based integration pipeline.

- `fig2_harmony_integration.R`  
  Core R implementation of the transcriptomic integration workflow.

- `configs/fig2_combined.yaml`  
  Combined configuration file. The active module for this figure is `modules.fig2_harmony`.

---

## Component 1: QuPath-based cell cropping

### Purpose

The QuPath script exports **fixed 200 × 200 pixel TIFF crops** centered on matching annotations. The script standardizes each ROI to a square crop around the annotation centroid and exports the transformed image patch using the current viewer display settings, including gamma correction.

This step is intended for image-level preparation and visualization associated with Fig. 2.

### Script

- `cell_cropping.groovy`

### What the script does

The script:

- checks for an active QuPath image viewer and image;
- prompts the user for:
  - output directory,
  - filename prefix,
  - z-slice index,
  - annotation naming pattern (`A1/A2/...` or `1/2/...`);
- identifies matching annotations;
- replaces each matched ROI with a centered **200 × 200 pixel square**;
- applies current viewer display transforms and gamma correction;
- exports each crop as a TIFF file.

### Outputs

Typical outputs are:

- one TIFF file per matched annotation;
- file names based on annotation names and the user-specified prefix.

### Note

Because the QuPath interaction includes annotation preparation, viewer settings, and manual choices that depend on the study design, the **full operational procedure is documented in the article protocol rather than repeated here**.

---

## Component 2: Transcriptomic integration

### Purpose

The transcriptomic component integrates multiple per-sample Seurat objects using a **two-pass SCTransform + Harmony** workflow and produces a single integrated Seurat object.

### Input requirements

The integration pipeline expects a directory containing **one Seurat object (`.rds`) per sample**.

Each input object should:

- be a valid Seurat object;
- contain an `RNA` assay;
- use a compatible gene naming scheme across samples.

Sample identifiers are inferred from file names (extension removed) and are written to metadata as `sample`.

---

## Integration workflow summary

### Pass 1

1. Read all per-sample Seurat objects.
2. Add sample identifiers to metadata.
3. Run `SCTransform()` independently for each sample.
4. Select integration features with `SelectIntegrationFeatures()`.
5. Exclude ribosomal/mitochondrial genes and optional user-defined genes.
6. Merge all samples.
7. Run `SCTransform()` on the merged object using `residual.features`.
8. Run PCA, Harmony, neighbor graph construction, clustering, and UMAP.
9. Remove predefined clusters interpreted as contamination or unwanted populations.

### Pass 2

1. Split the filtered merged object by sample.
2. Re-run per-sample `SCTransform()`.
3. Recompute integration features.
4. Run a second merged `SCTransform()` using filtered residual features.
5. Run PCA, Harmony, neighbor graph construction, clustering, and UMAP again.
6. Save the final integrated Seurat object.

---

## Outputs

The main integration output is a single integrated Seurat object:

- `out` — integrated Seurat object (`.rds`), typically containing:
  - Harmony reduction (`harmony`)
  - UMAP embedding (`umap`)
  - cluster labels (`seurat_clusters`)

For reproducibility, the pipeline also writes:

- `params_used_fig2_harmony.yaml` — resolved parameters used for the run
- `sessionInfo_fig2.txt` — R session and package information

These files are written to the same output directory as the integrated object.

---

## Recommended usage

Run from the repository root:

```bash
bash Fig2/run_fig2_combined.sh \
  -i /path/to/per_sample_rds_dir \
  -o /path/to/output_dir
```

The convenience behavior of `-o/--out` is:

- if the value ends with `.rds`, it is treated as the full output file path;
- otherwise, it is treated as an output directory and the pipeline writes `obj_oo.rds` inside that directory.

Example with explicit parameter overrides:

```bash
bash Fig2/run_fig2_combined.sh \
  -i /path/to/per_sample_rds_dir \
  -o /path/to/output_dir \
  --set seed=1 \
  --set res_pass2=1.2 \
  --set theta_pass2=1
```

---

## Key configuration fields

The following parameters are read from `modules.fig2_harmony` in `configs/fig2_combined.yaml`.

| Parameter | Description |
|---|---|
| `rds_dir` | Directory containing one Seurat `.rds` object per sample |
| `out` | Output path for the integrated Seurat object |
| `pattern` | Regular expression used to select input files |
| `seed` | Random seed |
| `sct_vfeatures_n` | Number of variable features used in per-sample `SCTransform()` |
| `nfeatures_integrate` | Number of integration features selected by `SelectIntegrationFeatures()` |
| `exclude_regex` | Regex pattern for genes to exclude from integration features |
| `exclude_genes` | Additional user-defined genes to exclude |
| `dims_pass1` | PCA/Harmony dimensions used in pass 1 |
| `res_pass1` | Clustering resolution in pass 1 |
| `remove_clusters` | Cluster IDs removed after pass 1 |
| `harmony_vars_pass1` | Batch variable(s) used by Harmony in pass 1 |
| `dims_pass2` | PCA/Harmony dimensions used in pass 2 |
| `res_pass2` | Clustering resolution in pass 2 |
| `theta_pass2` | Harmony `theta` used in pass 2 |
| `harmony_vars_pass2` | Batch variable(s) used by Harmony in pass 2 |

---

## Running the R script directly

The R script expects a **module-level YAML** file containing keys such as `rds_dir`, `out`, and other Fig. 2 parameters at the top level.

If your repository uses the combined YAML structure, extract the module first and then run:

```bash
python - <<'PY'
import yaml
root = yaml.safe_load(open('Fig2/configs/fig2_combined.yaml', 'r', encoding='utf-8'))
mod = root['modules']['fig2_harmony']
yaml.safe_dump(mod, open('Fig2/configs/fig2_harmony.yaml', 'w', encoding='utf-8'),
               sort_keys=False, allow_unicode=True)
print('wrote: Fig2/configs/fig2_harmony.yaml')
PY

Rscript Fig2/fig2_harmony_integration.R --config Fig2/configs/fig2_harmony.yaml
```

---

## Software requirements

### Core requirements

- Bash
- Python 3
- R (recommended: R >= 4.2)
- QuPath (for `cell_cropping.groovy`)

### Python packages

- `pyyaml`

Install with:

```bash
pip install pyyaml
```

### R packages

- `Seurat`
- `harmony`
- `dplyr`
- `patchwork`
- `ggplot2`
- `ggrepel`
- `yaml`

---

## Troubleshooting

- **No RDS files found**: check `rds_dir` and `pattern`.
- **Missing `pyyaml`**: install `pyyaml` in the Python environment used by the wrapper.
- **Missing R package `yaml`**: install it with `install.packages('yaml')`.
- **Cluster removal error after pass 1**: use the current `fig2_harmony_integration.R`, which removes clusters using Seurat identities rather than unstable expression-based subsetting.
- **Merge or SCTransform failure**: verify that each input file is a valid Seurat object and that genes are named consistently across samples.
- **QuPath export produces no files**: verify that matching annotations are present and that the chosen naming pattern matches the annotation names.

---

## Reproducibility

For each integration run, the pipeline records both resolved parameters and the R software environment. Users are encouraged to archive:

- the exact YAML used for the run,
- the generated `params_used_fig2_harmony.yaml`, and
- `sessionInfo_fig2.txt`.

For the image-cropping component, readers should refer to the **article protocol** for the full annotation and export procedure used in the manuscript.
