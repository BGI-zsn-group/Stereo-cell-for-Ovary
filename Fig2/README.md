# Fig. 2 | Two-pass SCTransform + Harmony integration

## Overview

This directory contains the code used to reproduce the integration analysis shown in **Fig. 2**. The workflow integrates multiple per-sample Seurat objects using a **two-pass SCTransform + Harmony** strategy and writes a single integrated Seurat object for downstream analyses and figure generation.

The implementation consists of:

- `run_fig2_combined.sh`: recommended entry point; extracts the Fig. 2 module from the combined YAML configuration and runs the integration pipeline.
- `fig2_harmony_integration.R`: core R implementation of the integration workflow.
- `configs/fig2_combined.yaml`: combined configuration file; the active module is `modules.fig2_harmony`.

## Input requirements

The pipeline expects a directory containing **one Seurat object (`.rds`) per sample**.

Each input object should:

- be a valid Seurat object;
- contain an `RNA` assay;
- use a compatible gene naming scheme across samples.

Sample identifiers are inferred from file names (extension removed) and added to the metadata as `sample`.

## Workflow summary

### Pass 1

1. Read all per-sample Seurat objects.
2. Add sample identifiers to metadata.
3. Run `SCTransform()` independently for each sample.
4. Select integration features with `SelectIntegrationFeatures()`.
5. Exclude ribosomal/mitochondrial and optionally user-specified genes.
6. Merge all samples.
7. Run `SCTransform()` on the merged object using `residual.features`.
8. Run PCA, Harmony, neighbor graph construction, clustering, and UMAP.
9. Remove pre-specified clusters interpreted as contamination or undesired populations.

### Pass 2

1. Split the filtered merged object by sample.
2. Re-run per-sample `SCTransform()`.
3. Recompute integration features.
4. Run a second merged `SCTransform()` using filtered residual features.
5. Run PCA, Harmony, neighbor graph construction, clustering, and UMAP again.
6. Save the final integrated Seurat object.

## Outputs

The main output is a single integrated Seurat object:

- `out` — integrated Seurat object (`.rds`), typically containing:
  - Harmony reduction (`harmony`)
  - UMAP embedding (`umap`)
  - cluster labels (`seurat_clusters`)

For reproducibility, the pipeline also writes:

- `params_used_fig2_harmony.yaml` — resolved parameters used for the run
- `sessionInfo_fig2.txt` — R session and package information

These files are written to the same output directory as the integrated object.

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

## Software requirements

### Core requirements

- Bash
- Python 3
- R (recommended: R >= 4.2)

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

## Troubleshooting

- **No RDS files found**: check `rds_dir` and `pattern`.
- **Missing `pyyaml`**: install `pyyaml` in the Python environment used by the wrapper.
- **Missing R package `yaml`**: install it with `install.packages('yaml')`.
- **Cluster removal error after pass 1**: use the current `fig2_harmony_integration.R`, which removes clusters using Seurat identities rather than unstable expression-based subsetting.
- **Merge or SCTransform failure**: verify that each input file is a valid Seurat object and that genes are named consistently across samples.

## Reproducibility

For each run, the pipeline records both resolved parameters and the R software environment. Users are encouraged to archive:

- the exact YAML used for the run,
- the generated `params_used_fig2_harmony.yaml`, and
- `sessionInfo_fig2.txt`.

These files provide sufficient provenance to reproduce the integration settings used in the manuscript.
