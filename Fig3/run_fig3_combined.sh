#!/usr/bin/env bash
set -euo pipefail

# Fig3 pipeline runner (uses *combined* YAML)
#
# Steps:
#   1) Monocle3 pseudotime + gene modules
#   2) Length ~ Pseudotime plot (+ correlation)
#   3) GO enrichment
#   4) SCENIC Step1 (Seurat -> CSV)
#   5) SCENIC Step2 (CSV -> loom + pyscenic grn/ctx/aucell)
#   6) SCENIC downstream (RSS)
#   7) Seurat -> h5ad (for scTour input)
#   8) scTour TNODE
#
# From repo root:
#   bash Fig3/run_fig3.sh
#
# This script extracts the module config from:
#   Fig3/configs/fig3_combined.yaml
# and passes it to all scripts as a flat YAML.

COMBINED="Fig3/configs/fig3_combined.yaml"
MODULE_KEY="fig3"   # module name (derived from original YAML file: fig3.yaml)

SCRIPT1="Fig3/fig3_monocle3_modules.R"
SCRIPT2="Fig3/fig3_length_pseudotime_plot.R"
SCRIPT3="Fig3/fig3_go_enrichment.R"
SCRIPT4="Fig3/fig3_scenic_rds_to_csv.R"
SCRIPT5="Fig3/fig3_scenic_pyscenic.py"
SCRIPT6="Fig3/fig3_scenic_downstream.R"
SCRIPT7="Fig3/fig3_rds_to_h5ad.R"
SCRIPT8="Fig3/fig3_sctour_tnode.py"

PYTHON_BIN="${PYTHON_BIN:-python}"
TMP_CFG="$(mktemp -t fig3_cfg_XXXXXX.yaml)"

"$PYTHON_BIN" - <<'PY' "$COMBINED" "$MODULE_KEY" "$TMP_CFG"
import sys
try:
    import yaml
except ImportError as e:
    raise SystemExit("ERROR: Python package 'pyyaml' is required to parse combined YAML. Please install PyYAML.") from e

combined_path, module_key, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
cfg = yaml.safe_load(open(combined_path, "r", encoding="utf-8"))
mod = cfg.get("modules", {}).get(module_key)
if mod is None:
    raise SystemExit(f"ERROR: module '{module_key}' not found under 'modules' in {combined_path}")
yaml.safe_dump(mod, open(out_path, "w", encoding="utf-8"), sort_keys=False, allow_unicode=True)
PY

Rscript "$SCRIPT1" --config "$TMP_CFG"
Rscript "$SCRIPT2" --config "$TMP_CFG"
Rscript "$SCRIPT3" --config "$TMP_CFG"
Rscript "$SCRIPT4" --config "$TMP_CFG"
"$PYTHON_BIN" "$SCRIPT5" --config "$TMP_CFG"
Rscript "$SCRIPT6" --config "$TMP_CFG"
Rscript "$SCRIPT7" --config "$TMP_CFG"
"$PYTHON_BIN" "$SCRIPT8" --config "$TMP_CFG"

rm -f "$TMP_CFG"
echo "[OK] Fig3 finished using combined YAML: $COMBINED (module: $MODULE_KEY)"
