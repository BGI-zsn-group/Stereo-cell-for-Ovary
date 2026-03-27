# Fig. 4 | Spatial mapping and pathway analysis of granulosa cells

## Overview

This directory contains the full pipeline used to reproduce **Fig. 4**, integrating:

- granulosa cell refinement
- pathway activity analysis with ssGSEA
- trajectory inference
- spatial transcriptomics mapping with Cell2location


The workflow is modular and is executed through a unified wrapper script.

## Directory structure

```text
Fig4/
  configs/
    fig4_combined.yaml

  data/
    B04372C211.roi_cell_ids.txt

  fig4_gc_processing.R
  fig4_ssgsea_hallmark.R
  fig4_monocle3.R
  fig4_cell2location_sc_prep.R
  fig4_cell2location_spatial_prep.py
  fig4_cell2location_train.py

  run_fig4_combined.sh
```

## Pipeline overview

The pipeline is executed with:

```bash
bash Fig4/run_fig4_combined.sh -o results/Fig4
```

Supported modules:

- `gc`
- `ssgsea`
- `monocle3`
- `c2l_sc`
- `c2l_spatial`
- `c2l_train`

## 1. Granulosa cell processing (`gc`)

Script: `fig4_gc_processing.R`

This step:

- subsets granulosa cells from the full annotated object
- performs Harmony-based multi-round clustering
- removes undesired clusters
- assigns refined final `celltype` labels

Primary output:

```text
results/Fig4/gc_processing/obj_gr_newanno.rds
```

## 2. Hallmark pathway analysis (`ssgsea`)

Script: `fig4_ssgsea_hallmark.R`

This step:

- loads the refined granulosa cell Seurat object
- computes average expression by cell type
- runs GSVA/ssGSEA using MSigDB Hallmark gene sets
- writes pathway score matrices and heatmaps

Default grouping:

```r
group_by = "celltype"
```

Typical output directory:

```text
results/Fig4/ssgsea_hallmark
```

## 3. Trajectory inference (`monocle3`)

Script: `fig4_monocle3.R`

This step:

- converts the refined Seurat object to a Monocle3 object
- learns the trajectory graph
- orders cells in pseudotime
- writes pseudotime values back to the Seurat object

Primary outputs typically include:

- Monocle3 `cds` object
- Seurat object with pseudotime
- pseudotime table
- UMAP pseudotime plot

## 4. Cell2location reference preparation (`c2l_sc`)

Script: `fig4_cell2location_sc_prep.R`

This step:

- combines the refined GC object with non-GC somatic cells from the full object
- keeps the desired sample
- exports a merged `.h5ad` file for Cell2location

Typical output:

```text
results/Fig4/cell2location_sc_ref/somatic_1226.h5ad
```

## 5. Spatial preprocessing (`c2l_spatial`)

Script: `fig4_cell2location_spatial_prep.py`

This step:

- performs QC filtering
- aggregates duplicated genes if needed
- computes neighborhood graph and clustering
- prepares the spatial AnnData used for Cell2location

### ROI filtering

To exclude non-ovarian regions, the workflow supports a manual ROI cell list stored at:

```text
Fig4/data/B04372C211.roi_cell_ids.txt
```

Format:

```text
cell_id_1
cell_id_2
cell_id_3
...
```

Behavior:

- if the ROI file is provided, only ROI cells are retained
- if the ROI path is left empty, ROI filtering is skipped

This manual ROI step is intentionally kept explicit for reproducibility and biological transparency.

## 6. Cell2location training (`c2l_train`)

Script: `fig4_cell2location_train.py`

This step has two stages.

### 6.1 Reference model training

- input: single-cell reference `.h5ad`
- model: `RegressionModel`
- output: inferred average expression signatures

### 6.2 Spatial mapping

- input: spatial `.h5ad`
- model: `Cell2location`
- output: estimated cell-type abundance for each spatial location

Typical output directory:

```text
results/Fig4/cell2location_train
```

## Usage examples

### Run the full pipeline

```bash
bash Fig4/run_fig4_combined.sh -o results/Fig4
```

### Run only ssGSEA

```bash
bash Fig4/run_fig4_combined.sh \
  --only ssgsea \
  -o results/Fig4 \
  -i ssgsea_rds=results/Fig4/gc_processing/obj_gr_newanno.rds
```

### Run Cell2location-related steps

```bash
bash Fig4/run_fig4_combined.sh \
  --only c2l_sc,c2l_spatial,c2l_train \
  -o results/Fig4
```

## Dependencies

### R packages

- Seurat
- harmony
- monocle3
- SeuratWrappers
- GSVA
- msigdbr
- pheatmap
- yaml

### Python packages

- scanpy
- anndata
- pandas
- numpy
- cell2location
- scvi-tools
- squidpy (optional)

## Reproducibility notes

- the refined GC object stores final labels directly in `meta.data$celltype`
- the ROI filtering step is manual but fixed through the provided ROI file
- downstream results are reproducible given the same input data, config, and ROI file
