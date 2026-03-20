# Figure 1

## Overview

This directory contains the preprocessing workflow used for Figure 1. The workflow starts from a raw GEM file and QuPath-defined ROI annotations, extracts ROI-specific transcriptomic profiles, converts the cropped GEM into a Seurat object, and appends per-cell morphological length measurements (`Length_px`) to the Seurat metadata for downstream analysis.

The pipeline consists of three steps:

1. **ROI-based GEM cropping**
2. **Conversion of cropped GEM to Seurat object**
3. **Injection of per-cell `Length_px` metadata into the Seurat object**

---

## Input files

The workflow expects the following inputs:

- A raw GEM file (`*.gem` or `*.gem.gz`)
- A QuPath ROI annotation file in GeoJSON format
- A per-cell length table (`line_lengths_<sample>.csv`)

### Required columns

#### Raw GEM
The raw GEM should contain coordinate columns (by default `x` and `y`) and the expression fields needed for downstream conversion.

#### Cropped GEM
The cropped GEM generated in step 1 must contain at least:

- `geneID`
- `MIDCount`
- `Label`

#### Length table
The file `line_lengths_<sample>.csv` must contain at least:

- `Name`
- `Length_px`

`Name` must match the Seurat cell names used in the corresponding `.rds` object.

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

### 2. `fig1_gem2rds.R`
Converts a cropped GEM file into a Seurat object.

**Inputs**
- Cropped GEM from step 1
- Output `.rds` path
- Project name

**Output**
- `<sample>.rds`

This script constructs a sparse count matrix from the cropped GEM and creates a Seurat object for downstream analysis.

---

### 3. `fig1_add_lengthpx.R`
Adds `Length_px` to Seurat metadata and removes unmatched cells.

**Inputs**
- One or more sample names
- Directory containing `.rds` files
- Directory containing `line_lengths_<sample>.csv`

**Output**
- Updated `.rds` file(s), either overwritten or written with a suffix

This script merges cell-level morphological length information into the Seurat object metadata and drops cells without matched `Length_px` values.

---

## Recommended execution order

The recommended order is:

1. `fig1_cutgem.py`
2. `fig1_gem2rds.R`
3. `fig1_add_lengthpx.R`

---

## Example: single-sample run

Assume the following files are present:

- `A04086F2.gem.gz`
- `A04086F2_roi.geojson`
- `line_lengths_A04086F2.csv`

### Step 1: crop GEM by QuPath ROI

```bash
python fig1_cutgem.py \
  -g A04086F2.gem.gz \
  -l A04086F2_roi.geojson \
  -b 50 \
  -C x y \
  -o A04086F2 \
  --outdir gem_cut
```

### Step 2: convert cropped GEM to Seurat

```bash
Rscript fig1_gem2rds.R \
  gem_cut/A04086F2.QuPath.cut.gem.gz \
  A04086F2.rds \
  A04086F2
```

### Step 3: add `Length_px`

```bash
Rscript fig1_add_lengthpx.R \
  --samples A04086F2 \
  --rds-dir . \
  --csv-dir . \
  --overwrite
```

---

## Example: multiple samples

```bash
for s in A04086F2 A04087F4; do
  python fig1_cutgem.py \
    -g ${s}.gem.gz \
    -l ${s}_roi.geojson \
    -b 50 \
    -C x y \
    -o ${s} \
    --outdir gem_cut

  Rscript fig1_gem2rds.R \
    gem_cut/${s}.QuPath.cut.gem.gz \
    ${s}.rds \
    ${s}
done

Rscript fig1_add_lengthpx.R \
  --samples A04086F2,A04087F4 \
  --rds-dir . \
  --csv-dir . \
  --overwrite
```

---

## Notes

- Although some legacy file names may end in `.txt`, the corresponding scripts should be treated as **R scripts**.
- Sample naming must be consistent across the cropped GEM, `.rds`, and `line_lengths_<sample>.csv`.
- If the coordinate columns in the raw GEM are not `x` and `y`, update them using `-C`.
- If you do not want to overwrite the original `.rds`, use `--out-suffix` in `fig1_add_lengthpx.R`.

---

## Output summary

Typical outputs include:

- ROI-cropped GEM files
- Seurat `.rds` objects
- Seurat objects with appended `Length_px` metadata

These outputs can be used directly for downstream integration and trajectory analyses in subsequent figures.
