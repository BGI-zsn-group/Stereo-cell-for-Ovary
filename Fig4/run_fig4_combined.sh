#!/usr/bin/env bash
set -euo pipefail

usage(){
  cat <<'USAGE'
Usage:
  bash Fig4/run_fig4_combined.sh [--only STEP[,STEP...]|all] [-o OUT]
                                 [-i PATH | -i alias=PATH]...
                                 [--set key=value]...
                                 [--config PATH | --config-dir DIR]
                                 [--scripts-dir DIR]
                                 [--dry-run]

Steps:
  gc,ssgsea,monocle3,c2l_sc,c2l_spatial,c2l_train

Input shortcuts:
  - `-i PATH`
      Fast path for a single-step run, or for the compatible middle-step chain
      `ssgsea,monocle3,c2l_sc` (all will use the same GC RDS as input).
  - `-i alias=PATH`
      Explicit input override. Supported aliases:
        gc_rds              -> fig4_gc_processing.input_rds
        gc_result_rds       -> set the same GC result RDS for ssgsea/monocle3/c2l_sc
        ssgsea_rds          -> fig4_ssgsea_hallmark.input_rds
        monocle3_rds        -> fig4_monocle3.input_rds
        c2l_sc_gc_rds       -> fig4_cell2location_sc_prep.obj_gr_rds
        c2l_sc_total_rds    -> fig4_cell2location_sc_prep.obj_total_rds
        spatial_h5ad        -> fig4_cell2location_spatial_prep.input_h5ad
        train_sc_h5ad       -> fig4_cell2location_train.single_cell.input_h5ad
        train_spatial_h5ad  -> fig4_cell2location_train.spatial.input_h5ad

Notes:
  - `--set` supports nested keys with dots (e.g., pipeline_round2.cluster_resolution=1.2).
  - Module-scoped overrides are supported: <module>.<key>=<value>
      e.g. fig4_gc_processing.out_dir=results/Fig4_alt/gc_processing
    These will apply ONLY to the matching module; others are ignored.
  - `-o/--out` redirects ALL outputs under a base directory,
    and also rewires downstream auto-inputs to those new outputs.
  - Explicit `-i/--input` overrides win over `-o` auto-wiring.

Examples:
  # run all steps using default config
  bash Fig4/run_fig4_combined.sh

  # run only GC processing with a custom input
  bash Fig4/run_fig4_combined.sh --only gc -i data/Fig4/obj_all_withanno.rds

  # run middle steps from an existing GC object
  bash Fig4/run_fig4_combined.sh --only ssgsea,monocle3,c2l_sc -i results/Fig4/gc_processing/obj_gr_newanno.rds

  # run train with explicit sc + spatial inputs
  bash Fig4/run_fig4_combined.sh --only c2l_train \
    -i train_sc_h5ad=results/Fig4/cell2location_sc_ref/somatic_1226.h5ad \
    -i train_spatial_h5ad=results/Fig4/cell2location_spatial/B04372C211.cell2location_spatial_counts.h5ad

  # redirect all outputs
  bash Fig4/run_fig4_combined.sh -o results/Fig4_alt
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

declare -a OVERRIDES=()
declare -a INPUT_SPECS=()
declare -a SELECTED_STEPS=()

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
    -i|--input)
      [[ $# -ge 2 ]] || die "Missing value for $1"
      [[ -n "${2}" ]] || die "$1 cannot be empty"
      INPUT_SPECS+=("$2"); shift 2;;
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
    printf '+'; printf ' %q' "$@"; printf '\n'
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

  head = k.split('.', 1)[0]
  if head in module_keys:
    if head != module_key:
      continue
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

append_override(){
  local key="$1"; local val="$2"
  OVERRIDES+=("${key}=${val}")
}

set_primary_input_for_step(){
  local step="$1"; local path="$2"
  case "$step" in
    gc)          append_override "fig4_gc_processing.input_rds" "$path" ;;
    ssgsea)      append_override "fig4_ssgsea_hallmark.input_rds" "$path" ;;
    monocle3)    append_override "fig4_monocle3.input_rds" "$path" ;;
    c2l_sc)      append_override "fig4_cell2location_sc_prep.obj_gr_rds" "$path" ;;
    c2l_spatial) append_override "fig4_cell2location_spatial_prep.input_h5ad" "$path" ;;
    c2l_train)   append_override "fig4_cell2location_train.single_cell.input_h5ad" "$path" ;;
    *) die "cannot infer primary input mapping for step '$step'" ;;
  esac
}

alias_input_override(){
  local alias_name="$1"; local path="$2"
  case "$alias_name" in
    gc_rds|gc) append_override "fig4_gc_processing.input_rds" "$path" ;;
    gc_result_rds|gc_result|gc_obj_rds)
      append_override "fig4_ssgsea_hallmark.input_rds" "$path"
      append_override "fig4_monocle3.input_rds" "$path"
      append_override "fig4_cell2location_sc_prep.obj_gr_rds" "$path"
      ;;
    ssgsea_rds|ssgsea) append_override "fig4_ssgsea_hallmark.input_rds" "$path" ;;
    monocle3_rds|monocle3) append_override "fig4_monocle3.input_rds" "$path" ;;
    c2l_sc_gc_rds|c2l_sc) append_override "fig4_cell2location_sc_prep.obj_gr_rds" "$path" ;;
    c2l_sc_total_rds|c2l_sc_total) append_override "fig4_cell2location_sc_prep.obj_total_rds" "$path" ;;
    spatial_h5ad|c2l_spatial|spatial) append_override "fig4_cell2location_spatial_prep.input_h5ad" "$path" ;;
    train_sc_h5ad|c2l_train_sc|train_sc) append_override "fig4_cell2location_train.single_cell.input_h5ad" "$path" ;;
    train_spatial_h5ad|c2l_train_spatial|train_spatial) append_override "fig4_cell2location_train.spatial.input_h5ad" "$path" ;;
    *)
      die "unknown input alias '$alias_name'. Use -h to see supported aliases."
      ;;
  esac
}

parse_selected_steps(){
  SELECTED_STEPS=()
  if [[ "$ONLY" == "all" || -z "$ONLY" ]]; then
    SELECTED_STEPS=(gc ssgsea monocle3 c2l_sc c2l_spatial c2l_train)
  else
    IFS=',' read -r -a SELECTED_STEPS <<< "$ONLY"
    [[ "${#SELECTED_STEPS[@]}" -gt 0 ]] || die "no valid steps parsed from --only '$ONLY'"
  fi
}

all_steps_in_middle_gc_chain(){
  local s
  for s in "${SELECTED_STEPS[@]}"; do
    case "$s" in
      ssgsea|monocle3|c2l_sc) ;;
      *) return 1 ;;
    esac
  done
  return 0
}

apply_input_specs(){
  local spec alias_name path s
  [[ "${#INPUT_SPECS[@]}" -eq 0 ]] && return 0

  for spec in "${INPUT_SPECS[@]}"; do
    if [[ "$spec" == *=* ]]; then
      alias_name="${spec%%=*}"
      path="${spec#*=}"
      [[ -n "$alias_name" && -n "$path" ]] || die "bad -i alias=path spec: '$spec'"
      alias_input_override "$alias_name" "$path"
    else
      path="$spec"
      if [[ "${#SELECTED_STEPS[@]}" -eq 1 ]]; then
        set_primary_input_for_step "${SELECTED_STEPS[0]}" "$path"
      elif all_steps_in_middle_gc_chain; then
        for s in "${SELECTED_STEPS[@]}"; do
          set_primary_input_for_step "$s" "$path"
        done
      else
        die "plain '-i PATH' is ambiguous for steps: ${SELECTED_STEPS[*]}. Use named inputs like -i gc_result_rds=... or -i train_sc_h5ad=..."
      fi
    fi
  done
}

COMBINED="$(discover_config)"
parse_selected_steps

if [[ -n "$OUT_BASE" ]]; then
  OVERRIDES+=(
    "fig4_gc_processing.out_dir=$OUT_BASE/gc_processing"
    "fig4_gc_processing.out_round1_rds=$OUT_BASE/gc_processing/gc.round1.rds"
    "fig4_gc_processing.out_final_rds=$OUT_BASE/gc_processing/obj_gr_newanno.rds"

    "fig4_ssgsea_hallmark.input_rds=$OUT_BASE/gc_processing/obj_gr_newanno.rds"
    "fig4_ssgsea_hallmark.out_dir=$OUT_BASE/ssgsea_hallmark"

    "fig4_monocle3.input_rds=$OUT_BASE/gc_processing/obj_gr_newanno.rds"
    "fig4_monocle3.out_dir=$OUT_BASE/monocle3"
    "fig4_monocle3.out_cds_rds=$OUT_BASE/monocle3/gc.cds.rds"
    "fig4_monocle3.out_obj_rds=$OUT_BASE/monocle3/gc.obj_with_pseudotime.rds"
    "fig4_monocle3.out_pseudotime_csv=$OUT_BASE/monocle3/gc.pseudotime.csv"
    "fig4_monocle3.out_plot_pdf=$OUT_BASE/monocle3/gc.monocle3_pseudotime_umap.pdf"

    "fig4_cell2location_sc_prep.obj_gr_rds=$OUT_BASE/gc_processing/obj_gr_newanno.rds"
    "fig4_cell2location_sc_prep.out_dir=$OUT_BASE/cell2location_sc_ref"
    "fig4_cell2location_sc_prep.out_h5ad=$OUT_BASE/cell2location_sc_ref/somatic_1226.h5ad"
    "fig4_cell2location_sc_prep.out_merged_rds=$OUT_BASE/cell2location_sc_ref/somatic_1226.merged_for_cell2location.rds"

    "fig4_cell2location_spatial_prep.out_dir=$OUT_BASE/cell2location_spatial"

    "fig4_cell2location_train.out_dir=$OUT_BASE/cell2location_train"
    "fig4_cell2location_train.single_cell.input_h5ad=$OUT_BASE/cell2location_sc_ref/somatic_1226.h5ad"
    "fig4_cell2location_train.spatial.input_h5ad=$OUT_BASE/cell2location_spatial/B04372C211.cell2location_spatial_counts.h5ad"
  )
fi

apply_input_specs

step_gc(){          run_step gc          fig4_gc_processing               R  "$SCRIPTS_DIR/fig4_gc_processing.R"; }
step_ssgsea(){      run_step ssgsea      fig4_ssgsea_hallmark             R  "$SCRIPTS_DIR/fig4_ssgsea_hallmark.R"; }
step_monocle3(){    run_step monocle3    fig4_monocle3                    R  "$SCRIPTS_DIR/fig4_monocle3.R"; }
step_c2l_sc(){      run_step c2l_sc      fig4_cell2location_sc_prep       R  "$SCRIPTS_DIR/fig4_cell2location_sc_prep.R"; }
step_c2l_spatial(){ run_step c2l_spatial fig4_cell2location_spatial_prep  PY "$SCRIPTS_DIR/fig4_cell2location_spatial_prep.py"; }
step_c2l_train(){   run_step c2l_train   fig4_cell2location_train         PY "$SCRIPTS_DIR/fig4_cell2location_train.py"; }

run_only(){
  local s="$1"
  case "$s" in
    gc) step_gc ;;
    ssgsea) step_ssgsea ;;
    monocle3) step_monocle3 ;;
    c2l_sc) step_c2l_sc ;;
    c2l_spatial) step_c2l_spatial ;;
    c2l_train) step_c2l_train ;;
    *) die "unknown step '$s'" ;;
  esac
}

if [[ "$ONLY" == "all" || -z "$ONLY" ]]; then
  step_gc; step_ssgsea; step_monocle3; step_c2l_sc; step_c2l_spatial; step_c2l_train
else
  for s in "${SELECTED_STEPS[@]}"; do run_only "$s"; done
fi

echo "[OK] Fig4 finished."
