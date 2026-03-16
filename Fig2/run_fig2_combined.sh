\
#!/usr/bin/env bash
set -euo pipefail

# Fig2 runner (Seurat + Harmony) — stable, nounset-safe, supports CLI input/output overrides.
#
# Run with:
#   bash Fig2/run_fig2_combined.sh ...
#
# Requires:
#   - python + pyyaml (pip install pyyaml)
#   - Rscript available in PATH (or set R_BIN=/path/to/Rscript)

usage() {
  cat <<'EOF'
Fig2 runner

Usage:
  bash Fig2/run_fig2_combined.sh [--only harmony|all]
                                 [-i INPUTDIR] [-o OUTDIR]
                                 [--set key=value]...
                                 [--config PATH | --config-dir DIR]
                                 [--scripts-dir DIR]
                                 [--module-key KEY]
                                 [--script harmony=FILE]
                                 [--dry-run]
                                 [--print-config]
                                 [--keep-config]

Short flags:
  -i, --in   INPUTDIR   Common input directory override (repeatable; last wins).
                        Mapped to: rds_dir and io.rds_dir
  -o, --out  OUTDIR     Common output directory override (repeatable; last wins).
                        Mapped to: out_dir/output_dir/results_dir and out_rds=OUTDIR/obj_oo.rds

Overrides:
  --set key=value       Override any field inside the extracted module config. Repeatable.
                        Supports dotted keys, e.g. --set io.rds_dir=/path
                        Supports module prefix, e.g. --set fig2_harmony.io.rds_dir=/path

Examples:
  # Default run (uses Fig2/configs/fig2_combined.yaml)
  bash Fig2/run_fig2_combined.sh

  # Set input + output quickly
  bash Fig2/run_fig2_combined.sh -i /data/originalRDS -o results/Fig2

  # Use a specific config file
  bash Fig2/run_fig2_combined.sh --config Fig2/configs/fig2_combined.yaml

  # Precise override
  bash Fig2/run_fig2_combined.sh --set io.rds_dir=/data/originalRDS --set out_rds=results/Fig2/obj_oo.rds
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG_DIR="$SCRIPT_DIR/configs"

ONLY="all"
CONFIG=""
CONFIG_DIR="$DEFAULT_CONFIG_DIR"
SCRIPTS_DIR="$SCRIPT_DIR"
MODULE_KEY="fig2_harmony"
DRY_RUN=0
PRINT_CONFIG=0
KEEP_CONFIG=0

PYTHON_BIN="${PYTHON_BIN:-python}"
R_BIN="${R_BIN:-Rscript}"

# nounset-safe: always initialize
declare -A SCRIPT_OVERRIDE=()
declare -a OVERRIDES=()
declare -a IN_PATHS=()
declare -a OUT_DIRS=()

run_cmd () {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '+'
    printf ' %q' "$@"
    printf '\n'
  else
    "$@"
  fi
}

# -------------------------
# arg parsing
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY="${2:?ERROR: --only requires a value}"; shift 2;;
    --config) CONFIG="${2:?ERROR: --config requires a path}"; shift 2;;
    --config-dir) CONFIG_DIR="${2:?ERROR: --config-dir requires a dir}"; shift 2;;
    --scripts-dir) SCRIPTS_DIR="${2:?ERROR: --scripts-dir requires a dir}"; shift 2;;
    --module-key) MODULE_KEY="${2:?ERROR: --module-key requires a key}"; shift 2;;
    --script)
      kv="${2:?ERROR: --script expects step=file}"; shift 2
      [[ "$kv" == *=* ]] || { echo "ERROR: --script expects step=file"; exit 2; }
      SCRIPT_OVERRIDE["${kv%%=*}"]="${kv#*=}"
      ;;
    --set) OVERRIDES+=("${2:?ERROR: --set requires key=value}"); shift 2;;
    -i|--in) IN_PATHS+=("${2:?ERROR: -i/--in requires a dir}"); shift 2;;
    -o|--out) OUT_DIRS+=("${2:?ERROR: -o/--out requires a dir}"); shift 2;;
    --dry-run) DRY_RUN=1; shift 1;;
    --print-config) PRINT_CONFIG=1; shift 1;;
    --keep-config) KEEP_CONFIG=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

# -------------------------
# Build overrides:
#   short flags first, explicit --set last (so --set wins)
# -------------------------
declare -a EXPLICIT_OVERRIDES=("${OVERRIDES[@]}")
OVERRIDES=()

if [[ "${#OUT_DIRS[@]}" -gt 0 ]]; then
  od="${OUT_DIRS[-1]}"
  OVERRIDES+=("out_dir=$od" "output_dir=$od" "results_dir=$od" "out_rds=$od/obj_oo.rds")
fi

if [[ "${#IN_PATHS[@]}" -gt 0 ]]; then
  ip="${IN_PATHS[-1]}"
  OVERRIDES+=("rds_dir=$ip" "io.rds_dir=$ip")
fi

OVERRIDES+=("${EXPLICIT_OVERRIDES[@]}")

discover_config () {
  if [[ -n "$CONFIG" ]]; then echo "$CONFIG"; return 0; fi
  local candidate="$CONFIG_DIR/fig2_combined.yaml"
  if [[ -f "$candidate" ]]; then echo "$candidate"; return 0; fi

  # fallback: if only one yaml in config dir, use it
  mapfile -t ys < <(ls -1 "$CONFIG_DIR"/*.yaml 2>/dev/null || true)
  if [[ "${#ys[@]}" -eq 1 ]]; then echo "${ys[0]}"; return 0; fi

  echo "ERROR: cannot find fig2_combined.yaml in: $CONFIG_DIR" >&2
  echo "       tip: pass --config or --config-dir" >&2
  exit 1
}

resolve_script () {
  local step="$1"; shift
  local override="${SCRIPT_OVERRIDE[$step]-}"
  if [[ -n "$override" ]]; then
    if [[ "$override" = /* ]]; then
      [[ -f "$override" ]] || { echo "ERROR: script not found: $override"; exit 1; }
      echo "$override"; return 0
    fi
    [[ -f "$SCRIPTS_DIR/$override" ]] || { echo "ERROR: script not found: $SCRIPTS_DIR/$override"; exit 1; }
    echo "$SCRIPTS_DIR/$override"; return 0
  fi
  for cand in "$@"; do
    [[ -f "$SCRIPTS_DIR/$cand" ]] && { echo "$SCRIPTS_DIR/$cand"; return 0; }
  done
  echo "ERROR: cannot find script for step '$step' in $SCRIPTS_DIR (tried: $*)" >&2
  exit 1
}

# -------------------------
# Extract module and apply overrides (filters empty args)
# -------------------------
extract_module () {
  local combined="$1"
  local module_key="$2"
  local out_cfg="$3"
  shift 3

  # filter empty overrides
  local -a clean=()
  local ov
  for ov in "$@"; do
    [[ -n "${ov:-}" ]] && clean+=("$ov")
  done

  "$PYTHON_BIN" - <<'PY' "$combined" "$module_key" "$out_cfg" "${clean[@]}"
import sys
try:
    import yaml
except Exception:
    raise SystemExit("ERROR: missing dependency pyyaml. Install: pip install pyyaml")

combined_path, module_key, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
overrides = sys.argv[4:]

root = yaml.safe_load(open(combined_path, "r", encoding="utf-8")) or {}
mods = root.get("modules") or {}
if module_key not in mods:
    raise SystemExit(f"ERROR: module '{module_key}' not found under `modules` in {combined_path}")
mod = mods[module_key]

def set_path(d, path, value):
    keys = [k for k in path.split(".") if k]
    cur = d
    for k in keys[:-1]:
        if not isinstance(cur, dict):
            raise SystemExit(f"ERROR: cannot set '{path}': '{k}' is not a mapping")
        if k not in cur or cur[k] is None:
            cur[k] = {}
        cur = cur[k]
    if not keys:
        raise SystemExit("ERROR: empty override key")
    cur[keys[-1]] = value

for ov in overrides:
    if not ov or "=" not in ov:
        raise SystemExit(f"ERROR: override must be key=value (got {ov})")
    k, v = ov.split("=", 1)
    k = k.strip()
    if k.startswith(module_key + "."):
        k = k[len(module_key) + 1 :]
    set_path(mod, k, v)

yaml.safe_dump(mod, open(out_path, "w", encoding="utf-8"), sort_keys=False, allow_unicode=True)
PY
}

COMBINED="$(discover_config)"
TMP_CFG="$(mktemp -t fig2_cfg_XXXXXX.yaml)"
extract_module "$COMBINED" "$MODULE_KEY" "$TMP_CFG" "${OVERRIDES[@]}"

[[ "$PRINT_CONFIG" -eq 1 ]] && echo "[config] $TMP_CFG"

SCRIPT="$(resolve_script harmony fig2_harmony_integration.R fig2_harmony_integration.txt fig2_harmony.R)"
case "$ONLY" in
  all|harmony|"") run_cmd "$R_BIN" "$SCRIPT" --config "$TMP_CFG";;
  *) echo "ERROR: unsupported --only '$ONLY' (supported: harmony|all)"; exit 2;;
esac

if [[ "$KEEP_CONFIG" -eq 0 ]]; then
  rm -f "$TMP_CFG"
fi

echo "[OK] Fig2 finished. combined=$COMBINED module=$MODULE_KEY only=$ONLY"
