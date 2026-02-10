#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Figs3 runner (combined YAML): somatic processing
#
# Recommended (from repo root):
#   bash Figs3/run_figs3_somatic_processing_combined.sh \
#     -i Figs3/configs/figs3_combined.yaml
#
# Optional:
#   -o/--out <dir|final.rds>      Override output directory or final RDS path
#   --set key=value              Override config values (supports nested keys with dots)
#
# Examples:
#   # write all outputs under a custom directory
#   bash Figs3/run_figs3_somatic_processing_combined.sh -o results/Figs3/somatic_processing_alt
#
#   # override a nested harmony parameter
#   bash Figs3/run_figs3_somatic_processing_combined.sh --set pipeline_round2.theta=2.0
# -----------------------------------------------------------------------------

die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat <<'USAGE'
Usage:
  bash Figs3/run_figs3_somatic_processing_combined.sh [options]

Options:
  -i, --in <combined.yaml>     Combined YAML (default: Figs3/configs/figs3_combined.yaml)
  -m, --module <name>          Module key under `modules` (default: figs3_somatic_processing)
  -s, --script <path>          R script path (default: Figs3/figs3_somatic_processing.R)
  -o, --out <dir|final.rds>    Override output directory or final RDS path
      --set key=value          Override config values (can repeat; supports dot paths)
  -h, --help                   Show help

Notes:
  - The wrapper requires Python + PyYAML to extract a module from the combined YAML.
  - The R script requires Seurat + harmony + yaml + dplyr + ggplot2 + patchwork.
USAGE
}

COMBINED_DEFAULT="Figs3/configs/figs3_combined.yaml"
MODULE_DEFAULT="figs3_somatic_processing"
SCRIPT_DEFAULT="Figs3/figs3_somatic_processing.R"

PYTHON_BIN="${PYTHON_BIN:-python}"

# Under `set -u`, arrays must be declared.
declare -a OV_SET=()

COMBINED="$COMBINED_DEFAULT"
MODULE_KEY="$MODULE_DEFAULT"
SCRIPT="$SCRIPT_DEFAULT"
OUT_OVERRIDE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--in)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      COMBINED="$2"; shift 2;;
    -m|--module)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      MODULE_KEY="$2"; shift 2;;
    -s|--script)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      SCRIPT="$2"; shift 2;;
    -o|--out)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      OUT_OVERRIDE="$2"; shift 2;;
    --set)
      [[ $# -ge 2 ]] || die "Missing value for --set"
      [[ -n "$2" ]] || die "Override cannot be empty"
      [[ "$2" == *"="* ]] || die "Override must be key=value (got '$2')"
      KEY="${2%%=*}"
      [[ -n "$KEY" ]] || die "Override key cannot be empty (got '$2')"
      OV_SET+=("$2")
      shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      die "Unknown arg: $1 (use -h for help)";;
  esac
done

[[ -f "$COMBINED" ]] || die "Combined YAML not found: $COMBINED"
[[ -f "$SCRIPT" ]] || die "R script not found: $SCRIPT"

TMP_CFG="$(mktemp -t figs3_cfg_XXXXXX.yaml)"
cleanup() { rm -f "$TMP_CFG"; }
trap cleanup EXIT

"$PYTHON_BIN" - "$COMBINED" "$MODULE_KEY" "$TMP_CFG" "$OUT_OVERRIDE" "${OV_SET[@]}" <<'PY'
import os, sys
try:
    import yaml
except ImportError as e:
    raise SystemExit("ERROR: Python package 'pyyaml' is required. Install with: pip install pyyaml") from e

combined_path, module_key, out_path, out_override = sys.argv[1:5]
overrides = sys.argv[5:]

with open(combined_path, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

mod = (cfg.get("modules") or {}).get(module_key)
if mod is None:
    raise SystemExit(f"ERROR: module '{module_key}' not found under 'modules' in {combined_path}")

def parse_value(v: str):
    # Best-effort type parsing using YAML itself
    try:
        return yaml.safe_load(v)
    except Exception:
        return v

def set_by_path(d, path, value):
    parts = path.split(".")
    cur = d
    for p in parts[:-1]:
        if p not in cur or not isinstance(cur[p], dict):
            cur[p] = {}
        cur = cur[p]
    cur[parts[-1]] = value

# Apply --out override first (so --set can still overwrite afterwards)
if out_override:
    if out_override.endswith(".rds"):
        out_dir = os.path.dirname(out_override) or "."
        mod["out_dir"] = out_dir
        mod["out_final_rds"] = out_override
        # Round1 RDS goes alongside final unless explicitly overridden later
        prefix = str(mod.get("prefix") or "somatic")
        mod["out_round1_rds"] = os.path.join(out_dir, f"{prefix}.round1.rds")
    else:
        # treat as directory
        out_dir = out_override
        mod["out_dir"] = out_dir
        prefix = str(mod.get("prefix") or "somatic")
        mod["out_round1_rds"] = os.path.join(out_dir, f"{prefix}.round1.rds")
        mod["out_final_rds"] = os.path.join(out_dir, f"{prefix}.final.rds")

# Apply --set overrides (supports nested keys via dots)
for ov in overrides:
    if not ov:
        continue
    if "=" not in ov:
        raise SystemExit(f"ERROR: override must be key=value (got '{ov}')")
    k, v = ov.split("=", 1)
    k = k.strip()
    if not k:
        raise SystemExit(f"ERROR: override key cannot be empty (got '{ov}')")
    set_by_path(mod, k, parse_value(v.strip()))

with open(out_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(mod, f, sort_keys=False, allow_unicode=True)
PY

Rscript "$SCRIPT" --config "$TMP_CFG"
echo "[OK] Figs3 somatic processing finished (combined: $COMBINED | module: $MODULE_KEY)"
