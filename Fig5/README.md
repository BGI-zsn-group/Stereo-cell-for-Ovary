# Fig5: Cell–cell communication analysis (CellChat)

This module performs cell–cell communication inference between granulosa cells (GC) and oocytes across developmental stages using CellChat.

---

## Overview

This workflow:

1. Loads processed GC and oocyte Seurat objects  
2. Filters unwanted cell populations  
   - GC: removes Atretic cells  
   - Oocyte: removes MII stage  
3. Splits data by developmental stage  
4. Runs CellChat independently for each stage  
5. Outputs one CellChat object per stage  

Note:  
We do NOT merge CellChat objects across stages. Each stage is analyzed independently.

---

## Input

### Required

- gr_final.rds  
  Seurat object of granulosa cells  

- Oocyte_final.rds  
  Seurat object of oocytes  

### Metadata requirements

#### GC object must contain:

- celltype (cell type annotation)

#### Oocyte object must contain:

- stage (developmental stage)

---

## Run

### Recommended

bash Stereo-cell-oocyte/Fig5/run_fig5_combined.sh \
  --gr-rds gr_final.rds \
  --oo-rds Oocyte_final.rds \
  -o Stereo-cell-oocyte/result/Fig5

### Optional

bash Stereo-cell-oocyte/Fig5/run_fig5_combined.sh \
  -i Stereo-cell-oocyte/Fig5/configs/fig5_combined.yaml \
  --gr-rds gr_final.rds \
  --oo-rds Oocyte_final.rds \
  -o Stereo-cell-oocyte/result/Fig5

---

## Output

### Per-stage CellChat objects

cellchat_EGO_<method>.rds  
cellchat_GO1_<method>.rds  
cellchat_GO2_<method>.rds  
cellchat_GO3_<method>.rds  
cellchat_FGO_<method>.rds  

Example:

cc <- readRDS("cellchat_EGO_triMean.rds")

---

### Stage index file

fig5.cellchat_stage_index_<method>.csv

---

## Downstream usage

cc <- readRDS("cellchat_EGO_triMean.rds")

table(cc@DB$interaction$annotation)

netVisual_circle(cc@net$count)

---

## Notes

### No merged CellChat object

We intentionally avoid mergeCellChat() to preserve stage-specific biology and simplify interpretation.

### Filtering strategy

- GC: remove "Atretic"  
- Oocyte: remove "MII"  

### Label usage

- celltype (GC)  
- stage (oocyte)  

---

## Dependencies

- Seurat  
- CellChat  
- dplyr  
- yaml  

---

## Workflow summary

GC (celltype) ─┐  
               ├── merge by stage ──> CellChat (per stage)  
Oocyte (stage) ┘  

EGO → CellChat  
GO1 → CellChat  
GO2 → CellChat  
GO3 → CellChat  
FGO → CellChat  
 
