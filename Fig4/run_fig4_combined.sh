#!/usr/bin/env bash
set -euo pipefail

# Fig4 runner (uses *combined* YAML)
#
# This script extracts module configs from:
#   Fig4/configs/fig4_combined.yaml
# then runs each step using a temporary flat YAML.
#
# From repo root:
#   bash Fig4/run_fig4.sh
#
# Dependencies:
# - python + pyyaml (to parse combined yaml)
# - Rscript (for Seurat/Monocle3 steps)
# - scanpy/anndata stack (for spatial AnnData prep)

COMBINED="Fig4/configs/fig4_combined.yaml"
PYTHON_BIN="${PYTHON_BIN:-python}"

# Helper: extract one module to a temp yaml
extract_module () {
  local module_key="$1"
  local tmp_cfg="$2"
  "$PYTHON_BIN" - <<'PY' "$COMBINED" "$module_key" "$tmp_cfg"
import sys
try:
    import yaml
except ImportError as e:
    raise SystemExit("ERROR: Python package 'pyyaml' is required (pip install pyyaml).") from e

combined_path, module_key, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
cfg = yaml.safe_load(open(combined_path, "r", encoding="utf-8"))
mod = cfg.get("modules", {}).get(module_key)
if mod is None:
    raise SystemExit(f"ERROR: module '{module_key}' not found under 'modules' in {combined_path}")
yaml.safe_dump(mod, open(out_path, "w", encoding="utf-8"), sort_keys=False, allow_unicode=True)
PY
}

# ---- 1) GC processing (Seurat/Harmony) ----
TMP1="$(mktemp -t fig4_gc_XXXXXX.yaml)"
extract_module "fig4_gc_processing" "$TMP1"
Rscript "Fig4/fig4_gc_processing.R" --config "$TMP1"
rm -f "$TMP1"

# ---- 2) Monocle3 on GC ----
TMP2="$(mktemp -t fig4_monocle3_XXXXXX.yaml)"
extract_module "fig4_monocle3" "$TMP2"
Rscript "Fig4/fig4_monocle3.R" --config "$TMP2"
rm -f "$TMP2"

# ---- 3) Cell2location single-cell reference prep (Seurat -> h5ad) ----
TMP3="$(mktemp -t fig4_c2l_sc_XXXXXX.yaml)"
extract_module "fig4_cell2location_sc_prep" "$TMP3"
Rscript "Fig4/fig4_cell2location_sc_prep.R" --config "$TMP3"
rm -f "$TMP3"

# ---- 4) Cell2location spatial prep (AnnData) ----
TMP4="$(mktemp -t fig4_c2l_spatial_XXXXXX.yaml)"
extract_module "fig4_cell2location_spatial_prep" "$TMP4"
"$PYTHON_BIN" "Fig4/fig4_cell2location_spatial_prep.py" --config "$TMP4"
rm -f "$TMP4"


# ---- 5) Cell2location training (reference + spatial mapping) ----
TMP5="$(mktemp -t fig4_c2l_train_XXXXXX.yaml)"
extract_module "fig4_cell2location_train" "$TMP5"
"$PYTHON_BIN" "Fig4/fig4_cell2location_train.py" --config "$TMP5"
rm -f "$TMP5"

echo "[OK] Fig4 pipeline finished using combined YAML: $COMBINED"
