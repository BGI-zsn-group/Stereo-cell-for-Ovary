#!/usr/bin/env bash
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat <<'USAGE'
Usage:
  bash Fig5/run_fig5_combined.sh [options]

Options:
  -i, --in <combined.yaml>     Combined YAML (default: Fig5/configs/fig5_combined.yaml)
  -m, --module <name>          Module key under `modules` (default: fig5_cellchat_build)
  -s, --script <path>          R script path (default: Fig5/fig5_cellchat_build.R)

  --gr-rds <path>              Override GC input RDS (modules.<module>.obj_gr_rds)
  --oo-rds <path>              Override oocyte input RDS (modules.<module>.obj_oo_rds)

  -o, --out <dir|prefix.rds>   Override output directory; if a *.rds path is given,
                               derive {out_dir, prefix} from it
  --prefix <name>              Override output prefix only

  --set key=value              Override config (repeatable; supports dot paths)
  -h, --help                   Show help

Examples:
  # default run
  bash Fig5/run_fig5_combined.sh

  # override both input RDS paths
  bash Fig5/run_fig5_combined.sh \
    --gr-rds data/Fig5/obj_gr_newanno.rds \
    --oo-rds data/Fig5/oocyte_edit.rds

  # override output directory
  bash Fig5/run_fig5_combined.sh \
    -o results/Fig5/cellchat_debug

  # override both inputs + output dir
  bash Fig5/run_fig5_combined.sh \
    -i Fig5/configs/fig5_combined.yaml \
    --gr-rds /path/to/gc.rds \
    --oo-rds /path/to/oocyte.rds \
    -o /path/to/outdir

  # derive out_dir and prefix from an rds-like output path
  bash Fig5/run_fig5_combined.sh \
    --gr-rds /path/to/gc.rds \
    --oo-rds /path/to/oocyte.rds \
    -o results/Fig5/cellchat/my_run.rds

  # keep generic override support
  bash Fig5/run_fig5_combined.sh --set mode=cellcell_contact
USAGE
}

COMBINED_DEFAULT="Fig5/configs/fig5_combined.yaml"
MODULE_DEFAULT="fig5_cellchat_build"
SCRIPT_DEFAULT="Fig5/fig5_cellchat_build.R"
PYTHON_BIN="${PYTHON_BIN:-python}"

declare -a OV_SET=()

COMBINED="$COMBINED_DEFAULT"
MODULE_KEY="$MODULE_DEFAULT"
SCRIPT="$SCRIPT_DEFAULT"
OUT_OVERRIDE=""
PREFIX_OVERRIDE=""
GR_RDS_OVERRIDE=""
OO_RDS_OVERRIDE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--in) [[ $# -ge 2 ]] || die "Missing value for $1"; COMBINED="$2"; shift 2;;
    -m|--module) [[ $# -ge 2 ]] || die "Missing value for $1"; MODULE_KEY="$2"; shift 2;;
    -s|--script) [[ $# -ge 2 ]] || die "Missing value for $1"; SCRIPT="$2"; shift 2;;
    --gr-rds) [[ $# -ge 2 ]] || die "Missing value for $1"; GR_RDS_OVERRIDE="$2"; shift 2;;
    --oo-rds) [[ $# -ge 2 ]] || die "Missing value for $1"; OO_RDS_OVERRIDE="$2"; shift 2;;
    -o|--out) [[ $# -ge 2 ]] || die "Missing value for $1"; OUT_OVERRIDE="$2"; shift 2;;
    --prefix) [[ $# -ge 2 ]] || die "Missing value for $1"; PREFIX_OVERRIDE="$2"; shift 2;;
    --set)
      [[ $# -ge 2 ]] || die "Missing value for --set"
      [[ -n "$2" ]] || die "Override cannot be empty"
      [[ "$2" == *"="* ]] || die "Override must be key=value (got '$2')"
      KEY="${2%%=*}"; [[ -n "$KEY" ]] || die "Override key cannot be empty (got '$2')"
      OV_SET+=("$2"); shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1 (use -h for help)";;
  esac
done

[[ -f "$COMBINED" ]] || die "Combined YAML not found: $COMBINED"
[[ -f "$SCRIPT" ]] || die "R script not found: $SCRIPT"

TMP_CFG="$(mktemp -t fig5_cfg_XXXXXX.yaml)"
cleanup() { rm -f "$TMP_CFG"; }
trap cleanup EXIT

"$PYTHON_BIN" - "$COMBINED" "$MODULE_KEY" "$TMP_CFG" "$OUT_OVERRIDE" "$PREFIX_OVERRIDE" "$GR_RDS_OVERRIDE" "$OO_RDS_OVERRIDE" "${OV_SET[@]}" <<'PY'
import os, sys
try:
    import yaml
except ImportError as e:
    raise SystemExit("ERROR: Python package 'pyyaml' is required. Install with: pip install pyyaml") from e

combined_path, module_key, out_path, out_override, prefix_override, gr_rds_override, oo_rds_override = sys.argv[1:8]
overrides = sys.argv[8:]

with open(combined_path, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

mod = (cfg.get("modules") or {}).get(module_key)
if mod is None:
    raise SystemExit(f"ERROR: module '{module_key}' not found under 'modules' in {combined_path}")

def parse_value(v: str):
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
echo "[OK] Fig5 CellChat build finished (combined: $COMBINED | module: $MODULE_KEY)"
