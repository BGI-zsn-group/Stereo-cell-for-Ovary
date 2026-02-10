# Single-cell transcriptomic atlas of mouse oocyte development from primary to peri-ovulatory follicles

This repository contains the analysis workflows used to generate the main and supplementary figures for the manuscript:

**_Single-cell transcriptomic atlas of mouse oocyte development from primary to peri-ovulatory follicles_**

The codebase is organized **by figure**. Each figure folder is intended to be runnable and reproducible on its own, with:
- a combined YAML config (`configs/*_combined.yaml`) describing inputs/parameters,
- a bash entrypoint (`run_*_combined.sh`) that extracts the relevant module(s) from the combined config and executes the corresponding R/Python scripts,
- figure-specific scripts (`*.R`, `*.py`),
- a figure-level `README.md` describing inputs/outputs and common troubleshooting.

> **Tip:** Start from the figure-level README. The top-level README is a navigation page + quickstart.

---

## Repository layout

```
.
├── Fig1/
├── Fig2/
├── Fig3/
├── Fig4/
├── Fig5/
├── Figs3/
├── .gitattributes
└── README.md
```

---

## Quickstart

### 1) Prerequisites

This repo mixes **R** and **Python**. Exact dependencies differ slightly by figure; each figure folder documents the required packages.

Typical requirements:

- **R (>= 4.2)**  
  Common packages: `Seurat`, `harmony`, `monocle3`, `CellChat`, `GSVA`, `msigdbr`, `pheatmap`, `yaml`, `dplyr`, `ggplot2`, `patchwork`, `reshape2`

- **Python (>= 3.9)**  
  Common packages: `pyyaml` (+ optional `scanpy/anndata`, `cell2location` for spatial workflows)

> If you hit a missing-package error, install the package indicated by the error message, then rerun.

### 2) Recommended data & output convention

To keep the configs portable (GitHub-friendly), we recommend:

- `data/` for **inputs**
- `results/` for **outputs**

Example:
```
data/
  Fig2/...
  Fig3/...
  Fig4/...
  Fig5/...
results/
  Fig2/...
  Fig3/...
  Fig4/...
  Fig5/...
```

Then update the input paths in each figure’s `configs/*_combined.yaml` to match your local files.

---

## Running the analyses

All commands below are intended to be run from the **repository root**.

### Fig2
```bash
bash Fig2/run_fig2_combined.sh -i Fig2/configs/fig2_combined.yaml
```

### Fig3
```bash
bash Fig3/run_fig3_combined.sh -i Fig3/configs/fig3_combined.yaml
```

### Fig4
```bash
bash Fig4/run_fig4_combined.sh -i Fig4/configs/fig4_combined.yaml
```

### Fig5
```bash
bash Fig5/run_fig5_combined.sh -i Fig5/configs/fig5_combined.yaml
```

### Supplementary Figs3
```bash
bash Figs3/run_figs3_somatic_processing_combined.sh -i Figs3/configs/figs3_combined.yaml
```

> For Fig1, see `Fig1/` for the figure-specific entrypoint and README.

---

## Common runner options

Most `run_*_combined.sh` scripts support a consistent set of patterns (see each figure README for the exact options):

- `--only <step>`: run a single step/module (when applicable)
- `-o/--out <dir|file>`: redirect outputs
- `--set key=value`: override config values from the command line (including nested keys)

Example:
```bash
bash Fig4/run_fig4_combined.sh --only ssgsea --set parallel_sz=16 -o results/Fig4/debug_ssgsea
```

---

## Reproducibility notes

1. **Paths:** The most common failure mode is missing input files. If you see “file not found”, update the relevant paths in the figure’s YAML config.
2. **Downstream dependencies:** Some steps consume outputs from prior steps. If you change output directories, update downstream inputs accordingly.
3. **Randomness:** When algorithms involve randomness (e.g., UMAP/integration/training), prefer fixed `seed` values and record your environment versions.

---

## Citation

If you use this repository (or parts of it) in your work, please cite the manuscript:

**Single-cell transcriptomic atlas of mouse oocyte development from primary to peri-ovulatory follicles**

(Full citation details to be added upon publication.)
