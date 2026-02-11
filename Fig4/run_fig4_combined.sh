#!/usr/bin/env bash
set -euo pipefail

usage(){
  cat <<'USAGE'
Usage:
  bash Fig4/run_fig4_combined.sh [--only STEP[,STEP...]|all] [-o OUT]
                                 [--set key=value]...
                                 [--config PATH | --config-dir DIR]
                                 [--scripts-dir DIR]
                                 [--dry-run]

Steps:
  gc,ssgsea,monocle3,c2l_sc,c2l_spatial,c2l_train

Notes:
  - `--set` supports nested keys with dots (e.g., pipeline_round2.cluster_resolution=1.2).
  - Module-scoped overrides are supported: <module>.<key>=<value>
      e.g. fig4_gc_processing.out_dir=results/Fig4_alt/gc_processing
    These will apply ONLY to the matching module; others are ignored.
  - `-o/--out` is a convenience option to redirect ALL outputs under a base directory,
    while keeping per-step subfolders to avoid collisions.

Examples:
  # run all steps using default config
  bash Fig4/run_fig4_combined.sh

  # run only GC processing
  bash Fig4/run_fig4_combined.sh --only gc

  # redirect all outputs
  bash Fig4/run_fig4_combined.sh -o results/Fig4_alt

  # override a nested parameter
  bash Fig4/run_fig4_combined.sh --set fig4_gc_processing.pipeline_round2.theta=2.0
USAGE
}

die(){ echo "ERROR: $*" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG_DIR="$SCRIPT_DIR/configs"

ONLY="all"
CONFIG=""
CONFIG_DIR="$DEFAULT_CONFIG_DIR"
SCRIPTS_DIR="$SCRIPT_DIR"
DRY_RUN=0

PYTHON_BIN="${PYTHON_BIN:-python}"
R_BIN="${R_BIN:-Rscript}"

# Under `set -u`, arrays must be declared.
declare -a OVERRIDES=()

OUT_BASE=""

need(){ [[ -f "$1" ]] || die "missing script: $1"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --only)
      [[ $# -ge 2 ]] || die "Missing value for --only"
      [[ -n "${2}" ]] || die "--only cannot be empty"
      ONLY="$2"; shift 2;;
    -o|--out)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      [[ -n "${2}" ]] || die "$1 cannot be empty"
      OUT_BASE="$2"; shift 2;;
    --set)
      [[ $# -ge 2 ]] || die "Missing value for --set"
      [[ -n "${2}" ]] || die "Override cannot be empty"
      [[ "${2}" == *"="* ]] || die "Override must be key=value (got '${2}')"
      KEY="${2%%=*}"
      [[ -n "$KEY" ]] || die "Override key cannot be empty (got '${2}')"
      OVERRIDES+=("$2"); shift 2;;
    --config)
      [[ $# -ge 2 ]] || die "Missing value for --config"
      CONFIG="$2"; shift 2;;
    --config-dir)
      [[ $# -ge 2 ]] || die "Missing value for --config-dir"
      CONFIG_DIR="$2"; shift 2;;
    --scripts-dir)
      [[ $# -ge 2 ]] || die "Missing value for --scripts-dir"
      SCRIPTS_DIR="$2"; shift 2;;
    --dry-run) DRY_RUN=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1 (use -h for help)";;
  esac
done

discover_config(){
  if [[ -n "$CONFIG" ]]; then echo "$CONFIG"; return 0; fi
  for cand in "$CONFIG_DIR/fig4_combined_v2.yaml" "$CONFIG_DIR/fig4_combined.yaml"; do
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  die "cannot find fig4_combined(.yaml|_v2.yaml) under: $CONFIG_DIR"
}

run(){
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '+'; printf ' %q' "$@"; printf '
'
  else
    "$@"
  fi
}

extract_module(){
  local combined="$1"; local module_key="$2"; local out_cfg="$3"; shift 3
  "$PYTHON_BIN" - <<'PY' "$combined" "$module_key" "$out_cfg" "$@"
import sys, os
try:
  import yaml
except Exception as e:
  raise SystemExit("ERROR: missing dependency pyyaml. Install: pip install pyyaml") from e

combined_path, module_key, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
overrides = sys.argv[4:]

data = yaml.safe_load(open(combined_path,'r',encoding='utf-8')) or {}
mods = data.get('modules') or {}
if module_key not in mods:
  raise SystemExit(f"ERROR: module '{module_key}' not found under modules in {combined_path}")
mod = mods[module_key]

def parse_value(v: str):
  v = v.strip()
  if v == "":
    return ""
  try:
    return yaml.safe_load(v)
  except Exception:
    return v

def set_path(d, path, val):
  ks=[k for k in path.split('.') if k]
  if not ks:
    raise SystemExit("ERROR: override key is empty after normalization")
  cur=d
  for k in ks[:-1]:
    if k not in cur or cur[k] is None or not isinstance(cur[k], dict):
      cur[k] = {}
    cur = cur[k]
  cur[ks[-1]] = val

module_keys = set(mods.keys())

for ov in overrides:
  if not ov:
    continue
  if '=' not in ov:
    raise SystemExit(f"ERROR: override must be key=value (got '{ov}')")
  k,v = ov.split('=',1)
  k = k.strip()
  if not k:
    raise SystemExit(f"ERROR: override key cannot be empty (got '{ov}')")

  # Module-scoped overrides: <module>.<key>=...
  head = k.split('.', 1)[0]
  if head in module_keys:
    if head != module_key:
      continue
    # strip "<module>."
    k = k.split('.', 1)[1] if '.' in k else ""
    if not k:
      raise SystemExit(f"ERROR: override key missing after module prefix (got '{ov}')")

  set_path(mod, k, parse_value(v))

yaml.safe_dump(mod, open(out_path,'w',encoding='utf-8'), sort_keys=False, allow_unicode=True)
PY
}

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

# -o convenience: redirect outputs under a base directory while keeping subfolders
if [[ -n "$OUT_BASE" ]]; then
  OVERRIDES+=(
    "fig4_gc_processing.out_dir=$OUT_BASE/gc_processing"
    "fig4_gc_processing.out_round1_rds=$OUT_BASE/gc_processing/gc.round1.rds"
    "fig4_gc_processing.out_final_rds=$OUT_BASE/gc_processing/obj_gr_newanno.rds"

    "fig4_monocle3.out_dir=$OUT_BASE/monocle3"
    "fig4_monocle3.out_cds_rds=$OUT_BASE/monocle3/gc.cds.rds"
    "fig4_monocle3.out_obj_rds=$OUT_BASE/monocle3/gc.obj_with_pseudotime.rds"
    "fig4_monocle3.out_pseudotime_csv=$OUT_BASE/monocle3/gc.pseudotime.csv"
    "fig4_monocle3.out_plot_pdf=$OUT_BASE/monocle3/gc.monocle3_pseudotime_umap.pdf"

    "fig4_cell2location_sc_prep.out_dir=$OUT_BASE/cell2location_sc_ref"
    "fig4_cell2location_sc_prep.out_h5ad=$OUT_BASE/cell2location_sc_ref/somatic_1226.h5ad"
    "fig4_cell2location_sc_prep.out_merged_rds=$OUT_BASE/cell2location_sc_ref/somatic_1226.merged_for_cell2location.rds"

    "fig4_cell2location_spatial_prep.out_dir=$OUT_BASE/cell2location_spatial"

    "fig4_cell2location_train.out_dir=$OUT_BASE/cell2location_train"
    "fig4_cell2location_train.single_cell.input_h5ad=$OUT_BASE/cell2location_sc_ref/somatic_1226.h5ad"
    "fig4_cell2location_train.spatial.input_h5ad=$OUT_BASE/cell2location_spatial/B04372C211.cell2location_spatial_counts.h5ad"
  )
fi

# Map steps -> (module_key, kind, script)
step_gc(){         run_step gc         fig4_gc_processing                 R  "$SCRIPTS_DIR/fig4_gc_processing.R"; }
step_ssgsea(){ run_step ssgsea    fig4_ssgsea_hallmark          R  "$SCRIPTS_DIR/fig4_ssgsea_hallmark.R"; }
step_monocle3(){   run_step monocle3   fig4_monocle3                      R  "$SCRIPTS_DIR/fig4_monocle3.R"; }
step_c2l_sc(){     run_step c2l_sc     fig4_cell2location_sc_prep         R  "$SCRIPTS_DIR/fig4_cell2location_sc_prep.R"; }
step_c2l_spatial(){run_step c2l_spatial fig4_cell2location_spatial_prep   PY "$SCRIPTS_DIR/fig4_cell2location_spatial_prep.py"; }
step_c2l_train(){  run_step c2l_train  fig4_cell2location_train           PY "$SCRIPTS_DIR/fig4_cell2location_train.py"; }

run_only(){
  local s="$1"
  case "$s" in
    gc) step_gc;;
    monocle3) step_monocle3;;
    c2l_sc) step_c2l_sc;;
    c2l_spatial) step_c2l_spatial;;
    c2l_train) step_c2l_train;;
    *) die "unknown step '$s'";;
  esac
}

if [[ "$ONLY" == "all" || -z "$ONLY" ]]; then
  step_gc; step_ssgsea; step_monocle3; step_c2l_sc; step_c2l_spatial; step_c2l_train
else
  IFS=',' read -r -a STEPS <<< "$ONLY"
  for s in "${STEPS[@]}"; do run_only "$s"; done
fi

echo "[OK] Fig4 finished."
