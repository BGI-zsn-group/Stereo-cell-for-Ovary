#!/usr/bin/env bash
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat <<'USAGE'
Usage:
  bash Fig5/run_fig5_combined.sh [options]

Options:
  -i, --in <combined.yaml>     Combined YAML.
                              Default: auto-resolved as <script_dir>/configs/fig5_combined.yaml
  -m, --module <name>          Module key under `modules` (default: fig5_cellchat_build)
  -s, --script <path>          R script path.
                              Default: auto-resolved as <script_dir>/fig5_cellchat_build.R

  --gr-rds <path>              Override GC input RDS (modules.<module>.obj_gr_rds)
  --oo-rds <path>              Override oocyte input RDS (modules.<module>.obj_oo_rds)

  -o, --out <dir|prefix.rds>   Override output directory; if a *.rds path is given,
                               derive {out_dir, prefix} from it
  --prefix <name>              Override output prefix only

  --set key=value              Override config (repeatable; supports dot paths inside module)
  -h, --help                   Show help

Behavior:
  - If you only specify the wrapper path, YAML and R script are auto-resolved relative to it.
  - You can still explicitly override YAML/R paths with -i / -s.
  - This wrapper extracts one module from the combined YAML and writes a temporary module YAML
    for the underlying R script.

Examples:
  # 1) Easiest: rely on auto-resolved YAML + R script relative to wrapper
  bash Stereo-cell-oocyte/Fig5/run_fig5_combined.sh \
    --gr-rds gr_final.rds \
    --oo-rds Oocyte_final.rds \
    -o Stereo-cell-oocyte/result/Fig5

  # 2) Explicit combined YAML
  bash Stereo-cell-oocyte/Fig5/run_fig5_combined.sh \
    -i Stereo-cell-oocyte/Fig5/configs/fig5_combined.yaml \
    --gr-rds gr_final.rds \
    --oo-rds Oocyte_final.rds \
    -o Stereo-cell-oocyte/result/Fig5

  # 3) Explicit R script
  bash Stereo-cell-oocyte/Fig5/run_fig5_combined.sh \
    -s Stereo-cell-oocyte/Fig5/fig5_cellchat_build.R \
    --gr-rds gr_final.rds \
    --oo-rds Oocyte_final.rds

  # 4) Direct output prefix inference from an .rds-like path
  bash Stereo-cell-oocyte/Fig5/run_fig5_combined.sh \
    --gr-rds gr_final.rds \
    --oo-rds Oocyte_final.rds \
    -o Stereo-cell-oocyte/result/Fig5/my_run.rds

  # 5) Generic config override
  bash Stereo-cell-oocyte/Fig5/run_fig5_combined.sh \
    --gr-rds gr_final.rds \
    --oo-rds Oocyte_final.rds \
    --set mode=cellcell_contact
USAGE
}

SCRIPT_PATH="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT_PATH" ]; do
  SCRIPT_DIR="$(cd -P "$(dirname "$SCRIPT_PATH")" && pwd)"
  SCRIPT_PATH="$(readlink "$SCRIPT_PATH")"
  [[ "$SCRIPT_PATH" != /* ]] && SCRIPT_PATH="$SCRIPT_DIR/$SCRIPT_PATH"
done
SCRIPT_DIR="$(cd -P "$(dirname "$SCRIPT_PATH")" && pwd)"

COMBINED_DEFAULT="$SCRIPT_DIR/configs/fig5_combined.yaml"
MODULE_DEFAULT="fig5_cellchat_build"
RSCRIPT_DEFAULT="$SCRIPT_DIR/fig5_cellchat_build.R"

PYTHON_BIN="${PYTHON_BIN:-python}"
R_BIN="${R_BIN:-Rscript}"

declare -a OV_SET=()

COMBINED="$COMBINED_DEFAULT"
MODULE_KEY="$MODULE_DEFAULT"
R_SCRIPT="$RSCRIPT_DEFAULT"
OUT_OVERRIDE=""
PREFIX_OVERRIDE=""
GR_RDS_OVERRIDE=""
OO_RDS_OVERRIDE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--in) [[ $# -ge 2 ]] || die "Missing value for $1"; COMBINED="$2"; shift 2;;
    -m|--module) [[ $# -ge 2 ]] || die "Missing value for $1"; MODULE_KEY="$2"; shift 2;;
    -s|--script) [[ $# -ge 2 ]] || die "Missing value for $1"; R_SCRIPT="$2"; shift 2;;
    --gr-rds) [[ $# -ge 2 ]] || die "Missing value for $1"; GR_RDS_OVERRIDE="$2"; shift 2;;
    --oo-rds) [[ $# -ge 2 ]] || die "Missing value for $1"; OO_RDS_OVERRIDE="$2"; shift 2;;
    -o|--out) [[ $# -ge 2 ]] || die "Missing value for $1"; OUT_OVERRIDE="$2"; shift 2;;
    --prefix) [[ $# -ge 2 ]] || die "Missing value for $1"; PREFIX_OVERRIDE="$2"; shift 2;;
    --set)
      [[ $# -ge 2 ]] || die "Missing value for --set"
      [[ -n "$2" ]] || die "Override cannot be empty"
      [[ "$2" == *"="* ]] || die "Override must be key=value (got '$2')"
      KEY="${2%%=*}"
      [[ -n "$KEY" ]] || die "Override key cannot be empty (got '$2')"
      OV_SET+=("$2")
      shift 2
      ;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1 (use -h for help)";;
  esac
done

[[ -f "$COMBINED" ]] || die "Combined YAML not found: $COMBINED"
[[ -f "$R_SCRIPT" ]] || die "R script not found: $R_SCRIPT"

TMP_CFG="$(mktemp -t fig5_cfg_XXXXXX.yaml)"
cleanup() { rm -f "$TMP_CFG"; }
trap cleanup EXIT

# Assemble python args safely under `set -u`
PY_ARGS=(
  "$COMBINED"
  "$MODULE_KEY"
  "$TMP_CFG"
  "$OUT_OVERRIDE"
  "$PREFIX_OVERRIDE"
  "$GR_RDS_OVERRIDE"
  "$OO_RDS_OVERRIDE"
)
if ((${#OV_SET[@]})); then
  PY_ARGS+=("${OV_SET[@]}")
fi

"$PYTHON_BIN" - "${PY_ARGS[@]}" <<'PY'
import os, sys
try:
    import yaml
except ImportError as e:
    raise SystemExit("ERROR: Python package 'pyyaml' is required. Install with: pip install pyyaml") from e

combined_path, module_key, out_path, out_override, prefix_override, gr_rds_override, oo_rds_override = sys.argv[1:8]
overrides = sys.argv[8:]

with open(combined_path, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

mods = cfg.get("modules") or {}
mod = mods.get(module_key)
if mod is None:
    raise SystemExit(f"ERROR: module '{module_key}' not found under 'modules' in {combined_path}")

def parse_value(v: str):
    try:
        return yaml.safe_load(v)
    except Exception:
        return v

def set_by_path(d, path, value):
    parts = [p for p in path.split(".") if p]
    if not parts:
        raise SystemExit("ERROR: empty override path")
    cur = d
    for p in parts[:-1]:
        if p not in cur or not isinstance(cur[p], dict):
            cur[p] = {}
        cur = cur[p]
    cur[parts[-1]] = value

# Explicit input overrides
if gr_rds_override:
    mod["obj_gr_rds"] = gr_rds_override
if oo_rds_override:
    mod["obj_oo_rds"] = oo_rds_override

# Output overrides
if out_override:
    if out_override.endswith(".rds"):
        out_dir = os.path.dirname(out_override) or "."
        prefix = os.path.splitext(os.path.basename(out_override))[0]
        mod["out_dir"] = out_dir
        mod["prefix"] = prefix
    else:
        mod["out_dir"] = out_override

if prefix_override:
    mod["prefix"] = prefix_override

# Generic module-local overrides
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

"$R_BIN" "$R_SCRIPT" --config "$TMP_CFG"
echo "[OK] Fig5 CellChat build finished"
echo "      wrapper:  $SCRIPT_PATH"
echo "      combined: $COMBINED"
echo "      module:   $MODULE_KEY"
echo "      rscript:  $R_SCRIPT"
echo "      tmp_cfg:  $TMP_CFG"
