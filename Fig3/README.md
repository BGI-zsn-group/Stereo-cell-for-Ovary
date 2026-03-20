# Fig. 3 | Pseudotime, gene modules, regulatory programs, and trajectory validation

## Overview

This directory contains the code used to reproduce the analyses shown in **Fig. 3**. The full workflow is organized as modular steps that can be executed jointly or independently through the wrapper script `run_fig3_combined.sh`.

The Fig. 3 workflow includes:

- **Monocle3 trajectory analysis** (`monocle3`)
- **Cell diameter versus pseudotime analysis** (`length`)
- **Gene Ontology enrichment of Monocle3 modules** (`go`)
- **SCENIC regulatory network analysis** (`scenic1`, `scenic2`, `scenic_downstream`)
- **scTour trajectory validation / alternative pseudotime modeling** (`rds2h5ad`, `sctour`)

## Directory contents

```text
Fig3/
  configs/
    fig3_combined.yaml
  fig3_monocle3_modules.R
  fig3_length_pseudotime.R
  fig3_go_enrichment.R
  fig3_scenic_rds_to_csv.R
  fig3_scenic_pyscenic.py
  fig3_scenic_downstream.R
  fig3_rds_to_h5ad.R
  fig3_sctour_tnode.py
  run_fig3_combined.sh
```

## Pipeline structure

### 1. Monocle3 pseudotime and module analysis (`monocle3`)

Input: integrated Seurat object from Fig. 2.

This step:

1. converts the Seurat object to a Monocle3 `cell_data_set`;
2. preprocesses cells and computes a reduced embedding;
3. reconstructs the Monocle3 object using Seurat cluster labels and UMAP coordinates;
4. learns the trajectory graph and orders cells in pseudotime from a user-specified root cluster;
5. identifies genes varying along the trajectory with `graph_test()`;
6. performs an initial gene-module call;
7. applies an optional stage-aware prefiltering strategy using positive markers from early and late stages;
8. reruns the final gene-module call on the filtered gene set;
9. writes the Seurat object augmented with pseudotime and binned pseudotime.

Primary outputs:

- `out_obj_rds` — Seurat object with `pseudotime` and `pseudotime_bin_q`
- `out_deg_csv` — trajectory-associated genes from Monocle3 `graph_test()`
- `out_pr_deg_ids` — gene IDs retained for the final module call
- `out_gene_modules_csv` — final gene-module assignments
- `out_prefilter_modules_csv` — prefilter module assignments (optional)
- `params_used_fig3_monocle3.yaml`
- `sessionInfo_fig3.txt`

### 2. Length versus pseudotime (`length`)

Input: Seurat object containing both `pseudotime` and `Length_px`.

This step fits a LOESS-smoothed trend of cell diameter versus pseudotime and writes:

- `length_vs_pseudotime.pdf`
- `length_vs_pseudotime.png`
- `length_vs_pseudotime_correlation.txt`

### 3. Gene Ontology enrichment (`go`)

Input: Monocle3 module assignment table.

This step merges user-defined module groups and performs GO enrichment analysis, returning a simplified GO results table.

### 4. SCENIC (`scenic1`, `scenic2`, `scenic_downstream`)

#### Step 4.1 — Export expression matrix (`scenic1`)

Input: Seurat object.

This step:

- optionally removes a specified stage (default: `MII`),
- extracts an expression matrix from the selected assay/slot,
- writes a CSV matrix for pySCENIC.

#### Step 4.2 — Run pySCENIC (`scenic2`)

Input: CSV from `scenic1`.

This step:

- converts the CSV matrix to `.loom`,
- runs pySCENIC (`grn`, `ctx`, `aucell`),
- preserves cell IDs from the exported CSV so that downstream matching is consistent.

Primary outputs:

- `scenic_out_csv`
- `scenic_loom`
- `scenic_adj_tsv`
- `scenic_reg_csv`
- `scenic_result_loom`

#### Step 4.3 — SCENIC downstream RSS analysis (`scenic_downstream`)

Input: SCENIC result loom and the corresponding Seurat object.

This step:

- loads `RegulonsAUC` from the SCENIC result loom,
- reloads the Seurat object used as SCENIC input,
- applies the same stage filtering used for SCENIC export,
- aligns cells between Seurat and SCENIC outputs,
- computes regulon specificity scores (RSS),
- saves the RSS matrix and publication-ready figures.

Primary outputs:

- `rss_matrix.csv`
- `rss_plot.pdf`
- `rss_plot.png`
- regulon and threshold RDS artifacts for reproducibility

### 5. scTour (`rds2h5ad`, `sctour`)

#### Step 5.1 — Convert Seurat to AnnData (`rds2h5ad`)

Input: Seurat object.

This step writes an `.h5ad` file suitable for scTour, typically under `out_dir/sctour/`.

#### Step 5.2 — Run scTour (`sctour`)

Input: `.h5ad` from `rds2h5ad`.

This step fits the scTour TNODE model and writes:

- updated `.h5ad` with `ptime` and latent embeddings,
- per-cell pseudotime output table.

## Recommended usage

Run all desired steps from the repository root through the wrapper.

### Monocle3 only

```bash
bash Fig3/run_fig3_combined.sh \
  --only monocle3 \
  -i /path/to/obj_oo.rds \
  -o results/Fig3
```

### Monocle3 + length + GO

```bash
bash Fig3/run_fig3_combined.sh \
  --only monocle3,length,go \
  -i /path/to/obj_oo.rds \
  -o results/Fig3
```

### SCENIC

```bash
bash Fig3/run_fig3_combined.sh \
  --only scenic1,scenic2,scenic_downstream \
  -i /path/to/obj_with_pseudotime.rds \
  -o results/Fig3 \
  --set scenic_tf_list=/path/to/allTFs_mm.txt \
  --set scenic_rankings_feather=/path/to/mm10_rankings.feather \
  --set scenic_motif_annotations=/path/to/motifs.tbl
```

### scTour

```bash
bash Fig3/run_fig3_combined.sh \
  --only rds2h5ad,sctour \
  -i /path/to/obj_with_pseudotime.rds \
  -o results/Fig3
```

## Output naming conventions

The wrapper interprets `-o` as the **Fig. 3 output directory**. Individual module outputs are generated under this directory.

For SCENIC, default file names are derived from:

1. `scenic_prefix` if explicitly provided; otherwise
2. the basename of the `-i` input file.

If stage exclusion is active (default excludes `MII`), the default SCENIC basename becomes:

```text
<input_basename>_without<excluded_stage>
```

For example, with:

```bash
-i result/oocyte_edit_1216.rds
```

SCENIC outputs will be written as:

```text
results/Fig3/scenic/oocyte_edit_1216_withoutMII.csv
results/Fig3/scenic/oocyte_edit_1216_withoutMII.loom
results/Fig3/scenic/adj_oocyte_edit_1216_withoutMII.tsv
results/Fig3/scenic/reg_oocyte_edit_1216_withoutMII.csv
results/Fig3/scenic/oocyte_edit_1216_withoutMII_result.loom
```

## Key configuration fields

Configuration is stored under `modules.fig3` in `configs/fig3_combined.yaml`.

### Core trajectory analysis

| Parameter | Description |
|---|---|
| `input_rds` | Input Seurat object |
| `out_dir` | Main output directory |
| `out_obj_rds` | Monocle3-updated Seurat object with pseudotime |
| `root_cluster` | Root cluster used for `order_cells()` |
| `num_dim` | Number of dimensions used for Monocle3 preprocessing |
| `cores` | Number of CPU cores for Monocle3 operations |
| `q_value_thr` | FDR threshold for `graph_test()` |
| `morans_I_thr` | Moran's I threshold for selecting trajectory-associated genes |
| `pseudotime_bins` | Number of pseudotime bins |
| `module_resolution` | Final resolution for `find_gene_modules()` |
| `module_seed` | Random seed for module identification |

### Stage-aware Monocle3 prefiltering

| Parameter | Description |
|---|---|
| `prefilter_enabled` | Enable the stage-aware prefiltering step |
| `prefilter_module_resolution` | Resolution used for the initial module call |
| `prefilter_target_module` | Module targeted for stage-aware filtering |
| `stage_col` | Metadata column used for stage labels |
| `early_stage_levels` | Stage labels considered early |
| `late_stage_levels` | Stage labels considered late |
| `out_prefilter_modules_csv` | Output path for the prefilter module table |

### Length analysis

| Parameter | Description |
|---|---|
| `plot_out_pdf`, `plot_out_png` | Output figure paths |
| `plot_loess_span` | LOESS smoothing span |
| `plot_cor_method` | Correlation method |

### SCENIC

| Parameter | Description |
|---|---|
| `scenic_input_rds` | Seurat object used for SCENIC |
| `scenic_stage_column` | Metadata column used for stage filtering |
| `scenic_exclude_stage_value` | Stage value excluded before SCENIC export |
| `scenic_assay`, `scenic_slot` | Assay and slot exported to CSV |
| `scenic_out_csv` | Expression matrix CSV for pySCENIC |
| `scenic_tf_list` | TF list used by pySCENIC |
| `scenic_rankings_feather` | cisTarget ranking database |
| `scenic_motif_annotations` | motif annotation table |
| `scenic_result_loom` | final SCENIC result loom |
| `scenic_rss_var` | metadata variable used for RSS grouping |

### scTour

| Parameter | Description |
|---|---|
| `rds2h5ad_input_rds` | Input Seurat object for AnnData conversion |
| `rds2h5ad_output_h5ad` | Output `.h5ad` path |
| `sctour_input_h5ad` | Input `.h5ad` for scTour |
| `sctour_output_h5ad` | Output `.h5ad` containing scTour results |
| `sctour_ptime_csv` | Per-cell pseudotime table |

## Software requirements

### General

- Bash
- Python 3
- R (recommended: R >= 4.2)

### Python packages

- `pyyaml`
- `scanpy`
- `loompy`
- `numpy`
- `pandas`
- `scipy`
- `sctour` (for the scTour step)
- `pySCENIC` CLI available as `pyscenic`

### R packages

- `Seurat`
- `SeuratWrappers`
- `monocle3`
- `yaml`
- `dplyr`
- `tidyverse`
- `ggplot2`
- `patchwork`
- `ggrepel`
- `reshape2`
- `clusterProfiler` and related GO dependencies (for `go`)
- `SCopeLoomR`, `AUCell`, `SCENIC`, `ComplexHeatmap`, `data.table`, `grid` (for SCENIC downstream)

## Troubleshooting

- **SCENIC result loom cannot be matched back to Seurat cells**: ensure `fig3_scenic_pyscenic.py` preserves cell IDs when reading the CSV, and rerun `scenic2` if older result looms were generated with numeric cell IDs.
- **No overlapping cells between Seurat and loom regulonAUC**: make sure the same stage filtering is applied in both SCENIC export and SCENIC downstream analysis.
- **Missing motif annotation file**: set `scenic_motif_annotations` to the full path of the motif annotation table.
- **scTour input not found**: verify the `.h5ad` path written by `rds2h5ad` and used by `sctour`.
- **Output appears under an unexpected file name**: check whether `scenic_prefix`, `-i`, and `scenic_exclude_stage_value` are influencing the derived basename.

## Reproducibility

The Fig. 3 pipeline writes parameter snapshots and environment metadata for key steps. For manuscript-level reproducibility, archive:

- the combined YAML used for the run,
- all generated `params_used_*.yaml` files,
- all `sessionInfo*.txt` files,
- and the exact external database versions used for SCENIC.

These artifacts document both analytical settings and software environments for each module of the figure.
