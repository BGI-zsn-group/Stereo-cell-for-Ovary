#!/usr/bin/env bash
set -euo pipefail

usage(){
  cat <<'USAGE'
Usage:
  bash Fig4/run_fig4_combined.sh [--only STEP[,STEP...]|all] [-o OUTDIR]
                                 [--set key=value]...
                                 [--config PATH | --config-dir DIR]
                                 [--scripts-dir DIR]
                                 [--dry-run]

Steps:
  gc,monocle3,c2l_sc,c2l_spatial,c2l_train
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG_DIR="$SCRIPT_DIR/configs"

ONLY="all"
CONFIG=""
CONFIG_DIR="$DEFAULT_CONFIG_DIR"
SCRIPTS_DIR="$SCRIPT_DIR"
DRY_RUN=0

PYTHON_BIN="${PYTHON_BIN:-python}"
R_BIN="${R_BIN:-Rscript}"

OVERRIDES=()
OUT_DIR=""


die(){ echo "ERROR: $*" >&2; exit 2; }
require_arg(){
  local opt="$1"; local val="${2-}"
  [[ -n "$val" ]] || die "$opt requires an argument"
}
add_override(){
  local ov="$1"
  [[ -n "$ov" ]] || die "--set requires key=value (got empty)"
  [[ "$ov" == *"="* ]] || die "override must be key=value (got: $ov)"
  local k="${ov%%=*}"
  [[ -n "$k" ]] || die "override key cannot be empty (got: $ov)"
  OVERRIDES+=("$ov")
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) require_arg "$1" "${2-}"; ONLY="$2"; shift 2;;
    -o|--out) require_arg "$1" "${2-}"; OUT_DIR="$2"; shift 2;;
    --set) require_arg "$1" "${2-}"; add_override "$2"; shift 2;;
    --config) require_arg "$1" "${2-}"; CONFIG="$2"; shift 2;;
    --config-dir) require_arg "$1" "${2-}"; CONFIG_DIR="$2"; shift 2;;
    --scripts-dir) require_arg "$1" "${2-}"; SCRIPTS_DIR="$2"; shift 2;;
    --dry-run) DRY_RUN=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

# -o convenience overrides
if [[ -n "$OUT_DIR" ]]; then
  OVERRIDES=("out_dir=$OUT_DIR" "output_dir=$OUT_DIR" "results_dir=$OUT_DIR" "${OVERRIDES[@]}")
fi

discover_config(){
  if [[ -n "$CONFIG" ]]; then echo "$CONFIG"; return 0; fi
  for cand in "$CONFIG_DIR/fig4_combined_v2.yaml" "$CONFIG_DIR/fig4_combined.yaml"; do
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  echo "ERROR: cannot find fig4_combined(.yaml|_v2.yaml) under: $CONFIG_DIR" >&2
  exit 1
}

run(){
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '+'; printf ' %q' "$@"; printf '\n'
  else
    "$@"
  fi
}

extract_module(){
  local combined="$1"; local module_key="$2"; local out_cfg="$3"; shift 3
  "$PYTHON_BIN" - <<'PY' "$combined" "$module_key" "$out_cfg" "$@"
import sys
try:
  import yaml
except Exception:
  raise SystemExit("ERROR: missing dependency pyyaml. Install: pip install pyyaml")
combined_path, module_key, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
overrides = sys.argv[4:]

data = yaml.safe_load(open(combined_path,'r',encoding='utf-8')) or {}
mods = data.get('modules') or {}
if module_key not in mods:
  raise SystemExit(f"ERROR: module '{module_key}' not found under modules in {combined_path}")
mod = mods[module_key]

def set_path(d, path, val):
  ks=[k for k in path.split('.') if k]
  cur=d
  for k in ks[:-1]:
    if k not in cur or cur[k] is None:
      cur[k] = {}
    cur = cur[k]
  cur[ks[-1]] = val

for ov in overrides:
  if not ov:
    continue
  if '=' not in ov:
    raise SystemExit(f"ERROR: override must be key=value (got: {ov})")
  k,v = ov.split('=',1)
  k=k.strip()
  if k.startswith(module_key+'.'):
    k = k[len(module_key)+1:]
  set_path(mod, k, v)

yaml.safe_dump(mod, open(out_path,'w',encoding='utf-8'), sort_keys=False, allow_unicode=True)
PY
}

need(){ [[ -f "$1" ]] || { echo "ERROR: missing script: $1" >&2; exit 1; }; }

run_step(){
  local step="$1"; local module="$2"; local kind="$3"; local script="$4"
  need "$script"
  local tmp_cfg
  tmp_cfg="$(mktemp -t fig4_${step}_XXXXXX.yaml)"
  extract_module "$COMBINED" "$module" "$tmp_cfg" "${OVERRIDES[@]}"
  if [[ "$kind" == "R" ]]; then
    run "$R_BIN" "$script" --config "$tmp_cfg"
  else
    run "$PYTHON_BIN" "$script" --config "$tmp_cfg"
  fi
  rm -f "$tmp_cfg"
}

COMBINED="$(discover_config)"

# Map steps -> (module_key, kind, script)
step_gc(){ run_step gc        fig4_gc_processing            R  "$SCRIPTS_DIR/fig4_gc_processing.R"; }
step_monocle3(){ run_step monocle3 fig4_monocle3                 R  "$SCRIPTS_DIR/fig4_monocle3.R"; }
step_c2l_sc(){ run_step c2l_sc   fig4_cell2location_sc_prep      R  "$SCRIPTS_DIR/fig4_cell2location_sc_prep.R"; }
step_c2l_spatial(){ run_step c2l_spatial fig4_cell2location_spatial_prep PY "$SCRIPTS_DIR/fig4_cell2location_spatial_prep.py"; }
step_c2l_train(){ run_step c2l_train fig4_cell2location_train    PY "$SCRIPTS_DIR/fig4_cell2location_train.py"; }

run_only(){
  local s="$1"
  case "$s" in
    gc) step_gc;;
    monocle3) step_monocle3;;
    c2l_sc) step_c2l_sc;;
    c2l_spatial) step_c2l_spatial;;
    c2l_train) step_c2l_train;;
    *) echo "ERROR: unknown step '$s'"; usage; exit 2;;
  esac
}

if [[ "$ONLY" == "all" || -z "$ONLY" ]]; then
  step_gc; step_monocle3; step_c2l_sc; step_c2l_spatial; step_c2l_train
else
  IFS=',' read -r -a STEPS <<< "$ONLY"
  for s in "${STEPS[@]}"; do run_only "$s"; done
fi

echo "[OK] Fig4 finished."
