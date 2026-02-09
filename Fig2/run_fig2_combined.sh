#!/usr/bin/env bash
set -euo pipefail

# Fig2 runner (uses *combined* YAML)
# From repo root:
#   bash Fig2/run_fig2.sh
#
# This script extracts the module config from:
#   Fig2/configs/fig2_combined.yaml
# and passes it to the R script as a flat YAML.

COMBINED="Fig2/configs/fig2_combined.yaml"
MODULE_KEY="fig2_harmony"   # module name (derived from original YAML file: fig2_harmony.yaml)

SCRIPT="Fig2/fig2_harmony_integration.R"
PYTHON_BIN="${PYTHON_BIN:-python}"

TMP_CFG="$(mktemp -t fig2_cfg_XXXXXX.yaml)"

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

Rscript "$SCRIPT" --config "$TMP_CFG"

rm -f "$TMP_CFG"
echo "[OK] Fig2 finished using combined YAML: $COMBINED (module: $MODULE_KEY)"
