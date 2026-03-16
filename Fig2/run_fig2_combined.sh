#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Fig2 runner

Usage:
  bash Fig2/run_fig2_combined.sh [--only harmony|all]
                                 [-i INPUTDIR] [-o OUT]
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

  -o, --out  OUT        Common output override (repeatable; last wins).
                        If OUT ends with .rds, treat it as the full output file path
                        and map to: out=OUT
                        Otherwise, treat it as an output directory and map to:
                        out=OUT/obj_oo.rds

Overrides:
  --set key=value       Override any field inside the extracted module config. Repeatable.
                        Supports dotted keys, e.g. --set io.rds_dir=/path
                        Supports module prefix, e.g. --set fig2_harmony.io.rds_dir=/path
                        Supports legacy out_rds=... which will be normalized to out=...

Examples:
  bash Fig2/run_fig2_combined.sh
  bash Fig2/run_fig2_combined.sh -i /data/originalRDS -o results/Fig2
  bash Fig2/run_fig2_combined.sh -i /data/originalRDS -o results/Fig2/my_obj.rds
  bash Fig2/run_fig2_combined.sh --set io.rds_dir=/data/originalRDS --set out=results/Fig2/obj_oo.rds
USAGE
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

declare -A SCRIPT_OVERRIDE=()
declare -a USER_OVERRIDES=()
declare -a EFFECTIVE_OVERRIDES=()
declare -a IN_PATHS=()
declare -a OUT_PATHS=()

run_cmd() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '+'
    printf ' %q' "$@"
    printf '\n'
  else
    "$@"
  fi
}

append_if_nonempty() {
  local arr_name="$1"
  local value="$2"
  [[ -n "$value" ]] || return 0
  eval "$arr_name+=(\"\$value\")"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY="${2:?ERROR: --only requires a value}"; shift 2 ;;
    --config) CONFIG="${2:?ERROR: --config requires a path}"; shift 2 ;;
    --config-dir) CONFIG_DIR="${2:?ERROR: --config-dir requires a dir}"; shift 2 ;;
    --scripts-dir) SCRIPTS_DIR="${2:?ERROR: --scripts-dir requires a dir}"; shift 2 ;;
    --module-key) MODULE_KEY="${2:?ERROR: --module-key requires a key}"; shift 2 ;;
    --script)
      kv="${2:?ERROR: --script expects step=file}"; shift 2
      [[ "$kv" == *=* ]] || { echo "ERROR: --script expects step=file" >&2; exit 2; }
      SCRIPT_OVERRIDE["${kv%%=*}"]="${kv#*=}"
      ;;
    --set) USER_OVERRIDES+=("${2:?ERROR: --set requires key=value}"); shift 2 ;;
    -i|--in) IN_PATHS+=("${2:?ERROR: -i/--in requires a dir}"); shift 2 ;;
    -o|--out) OUT_PATHS+=("${2:?ERROR: -o/--out requires a path}"); shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --print-config) PRINT_CONFIG=1; shift ;;
    --keep-config) KEEP_CONFIG=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

# Build overrides. Short flags first; explicit --set last so --set wins.
if (( ${#OUT_PATHS[@]} > 0 )); then
  op="${OUT_PATHS[-1]}"
  if [[ "$op" == *.rds ]]; then
    EFFECTIVE_OVERRIDES+=("out=$op")
  else
    op="${op%/}"
    EFFECTIVE_OVERRIDES+=("out=$op/obj_oo.rds")
  fi
fi

if (( ${#IN_PATHS[@]} > 0 )); then
  ip="${IN_PATHS[-1]}"
  EFFECTIVE_OVERRIDES+=("rds_dir=$ip" "io.rds_dir=$ip")
fi

if (( ${#USER_OVERRIDES[@]} > 0 )); then
  EFFECTIVE_OVERRIDES+=("${USER_OVERRIDES[@]}")
fi

discover_config() {
  if [[ -n "$CONFIG" ]]; then
    echo "$CONFIG"
    return 0
  fi

  local candidate="$CONFIG_DIR/fig2_combined.yaml"
  if [[ -f "$candidate" ]]; then
    echo "$candidate"
    return 0
  fi

  mapfile -t ys < <(ls -1 "$CONFIG_DIR"/*.yaml 2>/dev/null || true)
  if (( ${#ys[@]} == 1 )); then
    echo "${ys[0]}"
    return 0
  fi

  echo "ERROR: cannot find fig2_combined.yaml in: $CONFIG_DIR" >&2
  echo "       tip: pass --config or --config-dir" >&2
  exit 1
}

resolve_script() {
  local step="$1"
  shift
  local override="${SCRIPT_OVERRIDE[$step]-}"
  if [[ -n "$override" ]]; then
    if [[ "$override" = /* ]]; then
      [[ -f "$override" ]] || { echo "ERROR: script not found: $override" >&2; exit 1; }
      echo "$override"
      return 0
    fi
    [[ -f "$SCRIPTS_DIR/$override" ]] || { echo "ERROR: script not found: $SCRIPTS_DIR/$override" >&2; exit 1; }
    echo "$SCRIPTS_DIR/$override"
    return 0
  fi

  local cand
  for cand in "$@"; do
    if [[ -f "$SCRIPTS_DIR/$cand" ]]; then
      echo "$SCRIPTS_DIR/$cand"
      return 0
    fi
  done

  echo "ERROR: cannot find script for step '$step' in $SCRIPTS_DIR (tried: $*)" >&2
  exit 1
}

extract_module() {
  local combined="$1"
  local module_key="$2"
  local out_cfg="$3"
  shift 3

  local -a clean=()
  local ov
  for ov in "$@"; do
    [[ -n "${ov:-}" ]] && clean+=("$ov")
  done

  local -a py_args=("$combined" "$module_key" "$out_cfg")
  if (( ${#clean[@]} > 0 )); then
    py_args+=("${clean[@]}")
  fi

  "$PYTHON_BIN" - "${py_args[@]}" <<'PY'
import sys
from pathlib import Path

try:
    import yaml
except Exception:
    raise SystemExit("ERROR: missing dependency pyyaml. Install: pip install pyyaml")

_, combined_path, module_key, out_path, *overrides = sys.argv

with open(combined_path, "r", encoding="utf-8") as f:
    root = yaml.safe_load(f) or {}
mods = root.get("modules") or {}
if module_key not in mods:
    raise SystemExit(f"ERROR: module '{module_key}' not found under `modules` in {combined_path}")
mod = mods[module_key]
if mod is None:
    mod = {}


def set_path(d, path, value):
    keys = [k for k in path.split(".") if k]
    if not keys:
        raise SystemExit("ERROR: empty override key")
    cur = d
    for k in keys[:-1]:
        if not isinstance(cur, dict):
            raise SystemExit(f"ERROR: cannot set '{path}': '{k}' is not a mapping")
        if k not in cur or cur[k] is None:
            cur[k] = {}
        cur = cur[k]
    cur[keys[-1]] = value

for ov in overrides:
    if not ov or "=" not in ov:
        raise SystemExit(f"ERROR: override must be key=value (got {ov})")
    k, v = ov.split("=", 1)
    k = k.strip()
    if k.startswith(module_key + "."):
        k = k[len(module_key) + 1 :]
    if k == "out_rds":
        k = "out"
    set_path(mod, k, v)

with open(out_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(mod, f, sort_keys=False, allow_unicode=True)
PY
}

COMBINED="$(discover_config)"
TMP_CFG="$(mktemp -t fig2_cfg_XXXXXX.yaml)"

if (( ${#EFFECTIVE_OVERRIDES[@]} > 0 )); then
  extract_module "$COMBINED" "$MODULE_KEY" "$TMP_CFG" "${EFFECTIVE_OVERRIDES[@]}"
else
  extract_module "$COMBINED" "$MODULE_KEY" "$TMP_CFG"
fi

if [[ "$PRINT_CONFIG" -eq 1 ]]; then
  echo "[config] $TMP_CFG"
  if command -v sed >/dev/null 2>&1; then
    sed -n '1,200p' "$TMP_CFG"
  fi
fi

SCRIPT="$(resolve_script harmony fig2_harmony_integration.R fig2_harmony_integration.txt fig2_harmony.R)"
case "$ONLY" in
  all|harmony|"")
    run_cmd "$R_BIN" "$SCRIPT" --config "$TMP_CFG"
    ;;
  *)
    echo "ERROR: unsupported --only '$ONLY' (supported: harmony|all)" >&2
    exit 2
    ;;
esac

if [[ "$KEEP_CONFIG" -eq 0 ]]; then
  rm -f "$TMP_CFG"
fi

echo "[OK] Fig2 finished. combined=$COMBINED module=$MODULE_KEY only=$ONLY"
