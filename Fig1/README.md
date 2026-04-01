# Figure 1

## Overview

This directory contains the preprocessing workflow used for Figure 1. The workflow starts from a raw GEM file and QuPath-defined ROI annotations, extracts ROI-specific transcriptomic profiles, optionally generates GEM-based visualization / mask images, and converts the cropped GEM into a Seurat object for downstream analysis.

The workflow now consists of three components:

1. **ROI-based GEM cropping**
2. **GEM-based mask / visualization generation**
3. **Conversion of cropped GEM to Seurat object**

Detailed operational procedures, including QuPath preparation and the practical use of the GEM visualization script, are described in the **article protocol**. This README focuses on the repository structure, required inputs, and the role of each script.

---

## Input files

The workflow expects the following inputs, depending on which component is used:

- A raw GEM file (`*.gem` or `*.gem.gz`)
- A QuPath ROI annotation file in GeoJSON format
- Optional gene lists / channel assignments for GEM visualization

### Required columns

#### Raw GEM
The raw GEM should contain coordinate columns (by default `x` and `y`) and the expression fields needed for downstream conversion or visualization.

#### Cropped GEM
The cropped GEM generated in step 1 must contain at least:

- `geneID`
- `MIDCount`
- `Label`

For `gem_mask2.py`, the GEM input is expected to be interpretable as four columns corresponding to:

- `variable`
- `x`
- `y`
- `value`

The exact practical formatting and usage conventions are described in the **article protocol**.

---

## Scripts

### 1. `fig1_cutgem.py`
Crops a raw GEM file according to ROIs exported from QuPath.

**Inputs**
- GEM file
- ROI GeoJSON
- output prefix
- output directory

**Outputs**
- `<sample>.QuPath.cut.gem.gz`
- `<sample>.QuPath.roi.geojson`

This script rasterizes ROI polygons, assigns transcript coordinates to labeled ROIs, and outputs a GEM restricted to the selected regions.

---

### 2. `gem_mask2.py`
Generates GEM-based visualization images or masks.

**Purpose**
This script can be used to generate:

- merged RGB expression images,
- split per-gene channel images,
- total-count (`nCount`) images.

**Inputs**
- GEM file
- selected genes for red / green / blue channels
- optional cropping window
- optional normalization and contrast parameters
- output prefix

**Outputs**
- TIFF images written under `STimg/`

**Notes**
The script supports several modes, including:

- `all`
- `merge`
- `split`
- `nCount`

Because practical usage depends on the study design, selected markers, and QuPath/image coordination, the **detailed step-by-step operation is documented in the article protocol rather than repeated here**.

---

### 3. `fig1_gem2rds.R`
Converts a cropped GEM file into a Seurat object.

**Inputs**
- Cropped GEM from step 1
- Output `.rds` path
- Project name

**Output**
- `<sample>.rds`

This script constructs a sparse count matrix from the cropped GEM and creates a Seurat object for downstream analysis.

---

## Recommended execution order

The recommended order is:

1. `fig1_cutgem.py`
2. `gem_mask2.py` *(optional / protocol-dependent)*
3. `fig1_gem2rds.R`

---

## Example: single-sample run

Assume the following files are present:

- `A04086F2.gem.gz`
- `A04086F2_roi.geojson`

### Step 1: crop GEM by QuPath ROI

```bash
python fig1_cutgem.py   -g A04086F2.gem.gz   -l A04086F2_roi.geojson   -b 50   -C x y   -o A04086F2   --outdir gem_cut
```

### Step 2: optional GEM-based visualization / mask generation

The detailed operation of `gem_mask2.py`, including channel assignment, normalization, and visualization choices, is described in the **article protocol**.

### Step 3: convert cropped GEM to Seurat

```bash
Rscript fig1_gem2rds.R   gem_cut/A04086F2.QuPath.cut.gem.gz   A04086F2.rds   A04086F2
```

---

## Example: multiple samples

```bash
for s in A04086F2 A04087F4; do
  python fig1_cutgem.py     -g ${s}.gem.gz     -l ${s}_roi.geojson     -b 50     -C x y     -o ${s}     --outdir gem_cut

  Rscript fig1_gem2rds.R     gem_cut/${s}.QuPath.cut.gem.gz     ${s}.rds     ${s}
done
```

---

## Notes

- Although some legacy file names may end in `.txt`, the corresponding scripts should be treated as **R scripts** where appropriate.
- Sample naming must remain consistent across the cropped GEM and downstream `.rds` object.
- If the coordinate columns in the raw GEM are not `x` and `y`, update them using `-C`.
- `gem_mask2.py` is best treated as a figure-supporting visualization / masking utility rather than a mandatory preprocessing step for Seurat conversion.

---

## Output summary

Typical outputs include:

- ROI-cropped GEM files
- GEM-derived TIFF visualization / mask images
- Seurat `.rds` objects

These outputs can be used directly for downstream analyses and figure generation in subsequent steps of the study.
