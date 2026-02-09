#!/usr/bin/env bash
set -euo pipefail

# Run Fig.3 pipeline (linked I/O) + SCENIC (Step1+2+downstream) + scTour:
#   1) Monocle3 pseudotime + gene modules  -> out_obj_rds, gene_modules.csv
#   2) Length~Pseudotime plot (+ correlation) (reads out_obj_rds)
#   3) GO enrichment (reads gene_modules.csv)
#   4) SCENIC Step1: Seurat -> CSV (filters excluded stage) (reads out_obj_rds by default)
#   5) SCENIC Step2: CSV -> loom + pyscenic grn/ctx/aucell
#   6) SCENIC downstream: RSS by pseudotime_bin_q
#   7) Seurat out_obj_rds -> h5ad (for scTour input)
#   8) scTour TNODE pseudotime + latents (reads h5ad from step 7)
#
CONFIG="Fig3/configs/fig3_monocle3.yaml"

SCRIPT1="Fig3/fig3_monocle3_modules.R"
SCRIPT2="Fig3/fig3_length_pseudotime_plot.R"
SCRIPT3="Fig3/fig3_go_enrichment.R"
SCRIPT4="Fig3/fig3_scenic_rds_to_csv.R"
SCRIPT5="Fig3/fig3_scenic_pyscenic.py"
SCRIPT6="Fig3/fig3_scenic_downstream.R"
SCRIPT7="Fig3/fig3_rds_to_h5ad.R"
SCRIPT8="Fig3/fig3_sctour_tnode.py"

PYTHON_BIN="${PYTHON_BIN:-python}"

Rscript "$SCRIPT1" --config "$CONFIG"
Rscript "$SCRIPT2" --config "$CONFIG"
Rscript "$SCRIPT3" --config "$CONFIG"
Rscript "$SCRIPT4" --config "$CONFIG"
"$PYTHON_BIN" "$SCRIPT5" --config "$CONFIG"
Rscript "$SCRIPT6" --config "$CONFIG"
Rscript "$SCRIPT7" --config "$CONFIG"
"$PYTHON_BIN" "$SCRIPT8" --config "$CONFIG"
