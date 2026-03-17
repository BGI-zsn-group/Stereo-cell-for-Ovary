#!/usr/bin/env bash
set -euo pipefail

# Fig3 runner: run selected scripts using a module extracted from combined YAML,
# with optional CLI overrides for input/output.

usage(){
  cat <<'USAGE'
Usage:
  bash Fig3/run_fig3_combined.sh [--only STEP[,STEP...]|all] [-i INPUT_RDS] [-o OUTDIR|OUT_OBJ_RDS]
                                 [--set key=value]... [--config PATH|--config-dir DIR]
                                 [--scripts-dir DIR] [--script step=PATH]... [--module-key KEY]
                                 [--dry-run]
Steps: monocle3,length,go,scenic1,scenic2,scenic_downstream,rds2h5ad,sctour
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR_DEFAULT="$SCRIPT_DIR/configs"

ONLY="all"; CONFIG=""; CONFIG_DIR="$CONFIG_DIR_DEFAULT"; SCRIPTS_DIR="$SCRIPT_DIR"; MODULE_KEY="fig3"; DRY_RUN=0
PYTHON_BIN="${PYTHON_BIN:-python}"; R_BIN="${R_BIN:-Rscript}"

# overrides (short flags expanded to many key=val), explicit --set appended last
# IMPORTANT under `set -u`: declare arrays upfront so `${arr[@]}` never trips "unbound variable".
declare -a OV_SHORT=()
declare -a OV_EXPL=()
declare -a IN_PATHS=()
declare -a OUT_DIRS=()

# Back-compat: some older variants used these names; keep them declared so `set -u` never breaks.
declare -a USER_SET_OVERRIDES=()
declare -A SCRIPT_OVERRIDE=()

# ---- arg parse ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: --only needs a value" >&2; exit 2; }
      ONLY="$2"; shift 2;;
    --config)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: --config needs a value" >&2; exit 2; }
      CONFIG="$2"; shift 2;;
    --config-dir)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: --config-dir needs a value" >&2; exit 2; }
      CONFIG_DIR="$2"; shift 2;;
    --scripts-dir)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: --scripts-dir needs a value" >&2; exit 2; }
      SCRIPTS_DIR="$2"; shift 2;;
    --module-key)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: --module-key needs a value" >&2; exit 2; }
      MODULE_KEY="$2"; shift 2;;
    --script)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: --script needs step=PATH" >&2; exit 2; }
      kv="$2"; shift 2
      [[ "$kv" == *=* ]] || { echo "ERROR: --script needs step=PATH"; exit 2; }
      SCRIPT_OVERRIDE["${kv%%=*}"]="${kv#*=}";;
    --set)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: --set needs key=value" >&2; exit 2; }
      OV_EXPL+=("$2"); shift 2;;
    -i|--in)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: -i/--in needs a value" >&2; exit 2; }
      IN_PATHS+=("$2"); shift 2;;
    -o|--out)
      [[ $# -ge 2 && -n "${2:-}" ]] || { echo "ERROR: -o/--out needs a value" >&2; exit 2; }
      OUT_DIRS+=("$2"); shift 2;;
    --dry-run) DRY_RUN=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

# ---- helpers ----
has_step(){
  local want="$1"; [[ "$ONLY" == "all" || -z "$ONLY" ]] && return 0
  IFS=',' read -r -a _s <<< "$ONLY"; for s in "${_s[@]}"; do [[ "$s" == "$want" ]] && return 0; done
  return 1
}

discover_config(){
  if [[ -n "$CONFIG" ]]; then echo "$CONFIG"; return 0; fi
  for cand in "$CONFIG_DIR/fig3_combined_v2.yaml" "$CONFIG_DIR/fig3_combined.yaml"; do
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  echo "ERROR: cannot find fig3 combined yaml in $CONFIG_DIR" >&2; exit 1
}

resolve_script(){
  local step="$1"; shift
  local ov="${SCRIPT_OVERRIDE[$step]-}"
  if [[ -n "$ov" ]]; then
    [[ "$ov" = /* ]] && { [[ -f "$ov" ]] || { echo "ERROR: script not found: $ov"; exit 1; }; echo "$ov"; return 0; }
    [[ -f "$SCRIPTS_DIR/$ov" ]] || { echo "ERROR: script not found: $SCRIPTS_DIR/$ov"; exit 1; }
    echo "$SCRIPTS_DIR/$ov"; return 0
  fi
  for cand in "$@"; do [[ -f "$SCRIPTS_DIR/$cand" ]] && { echo "$SCRIPTS_DIR/$cand"; return 0; }; done
  echo "ERROR: cannot find script for step '$step' under $SCRIPTS_DIR" >&2; exit 1
}

run_cmd(){
  if [[ "$DRY_RUN" -eq 1 ]]; then printf '+'; printf ' %q' "$@"; printf '\n'; else "$@"; fi
}

# ---- expand short overrides ----
if [[ ${#OUT_DIRS[@]} -gt 0 ]]; then
  od="${OUT_DIRS[-1]}"
  # Convenience: if -o points to a .rds file, treat it as out_obj_rds and
  # derive out_dir from its parent directory.
  if [[ "$od" == *.rds ]]; then
    od_file="$od"
    od="$(dirname "$od_file")"
    OV_SHORT+=("out_obj_rds=$od_file")
  fi
  OV_SHORT+=("out_dir=$od" "output_dir=$od" "results_dir=$od")
  OV_SHORT+=("out_obj_rds=$od/obj_with_pseudotime.rds" "out_deg_csv=$od/deg_graph_test.csv" "out_pr_deg_ids=$od/pr_deg_ids.txt" "out_gene_modules_csv=$od/gene_modules.csv")
  OV_SHORT+=("plot_out_pdf=$od/length_vs_pseudotime.pdf" "plot_out_png=$od/length_vs_pseudotime.png" "plot_out_cor_txt=$od/length_vs_pseudotime_correlation.txt")
  OV_SHORT+=("go_out_csv=$od/GO_BP_merged_modules_simplified.csv")
  OV_SHORT+=("scenic_out_dir=$od/scenic" "scenic_out_csv=$od/scenic/oocyte_1211_withoutMII.csv" "scenic_loom=$od/scenic/oocyte_1211_withoutMII.loom")
  OV_SHORT+=("scenic_adj_tsv=$od/scenic/adj_oocyte_1211_withoutMII.tsv" "scenic_reg_csv=$od/scenic/reg_oocyte_1211_withoutMII.csv" "scenic_result_loom=$od/scenic/oocyte_1211_withoutMII_result.loom")
  OV_SHORT+=("scenic_downstream_out_dir=$od/scenic/downstream" "scenic_rss_csv=$od/scenic/downstream/rss_matrix.csv")
  OV_SHORT+=("scenic_rss_plot_pdf=$od/scenic/downstream/rss_plot.pdf" "scenic_rss_plot_png=$od/scenic/downstream/rss_plot.png")
  OV_SHORT+=("scenic_regulons_rds=$od/scenic/downstream/regulons_list.rds" "scenic_auc_thresholds_rds=$od/scenic/downstream/regulon_auc_thresholds.rds")
  OV_SHORT+=("rds2h5ad_output_h5ad=$od/sctour/seurat_to_h5ad.h5ad" "sctour_out_dir=$od/sctour" "sctour_output_h5ad=$od/sctour/oocyte_sctour_tnode.h5ad" "sctour_ptime_csv=$od/sctour/ptime.csv")
fi
if [[ ${#IN_PATHS[@]} -gt 0 ]]; then
  ip="${IN_PATHS[-1]}"
  OV_SHORT+=("input_rds=$ip" "rds2h5ad_input_rds=$ip")
fi

# ---- extract module + apply overrides ----
extract_module(){
  local combined="$1" module_key="$2" out_cfg="$3"; shift 3
  "$PYTHON_BIN" - <<'PY' "$combined" "$module_key" "$out_cfg" "$@"
import sys
try:
  import yaml
except Exception:
  raise SystemExit('ERROR: need pyyaml (pip install pyyaml)')
combined, mkey, outp = sys.argv[1], sys.argv[2], sys.argv[3]
ovs = sys.argv[4:]
root = yaml.safe_load(open(combined,'r',encoding='utf-8')) or {}
mods = root.get('modules') or {}
if mkey not in mods:
  raise SystemExit(f"ERROR: module '{mkey}' not found under modules in {combined}")
mod = mods[mkey]

def set_path(d, path, val):
  ks=[k for k in path.split('.') if k]
  cur=d
  for k in ks[:-1]:
    if k not in cur or cur[k] is None:
      cur[k] = {}
    cur = cur[k]
  cur[ks[-1]] = val

for ov in ovs:
  # Bash may pass an empty arg in some edge cases; ignore it so we don't error on "got ''".
  if not ov:
    continue
  if '=' not in ov:
    raise SystemExit(f"ERROR: override must be key=value (got {ov})")
  k,v = ov.split('=',1)
  k=k.strip()
  if k.startswith(mkey + '.'): k=k[len(mkey)+1:]
  set_path(mod,k,v)

yaml.safe_dump(mod, open(outp,'w',encoding='utf-8'), sort_keys=False, allow_unicode=True)
PY
}

COMBINED="$(discover_config)"
TMP_CFG="$(mktemp -t fig3_cfg_XXXXXX.yaml)"

# Safely build override args without ever expanding an unset array (set -u),
# and without injecting a spurious empty argument.
EXTRACT_ARGS=()
((${#OV_SHORT[@]})) && EXTRACT_ARGS+=("${OV_SHORT[@]}")
((${#OV_EXPL[@]}))  && EXTRACT_ARGS+=("${OV_EXPL[@]}")
((${#USER_SET_OVERRIDES[@]})) && EXTRACT_ARGS+=("${USER_SET_OVERRIDES[@]}")

extract_module "$COMBINED" "$MODULE_KEY" "$TMP_CFG" "${EXTRACT_ARGS[@]}"

run_one(){
  case "$1" in
    monocle3) s="$(resolve_script monocle3 fig3_monocle3_modules.R)"; run_cmd "$R_BIN" "$s" --config "$TMP_CFG";;
    length)   s="$(resolve_script length   fig3_length_pseudotime.R)"; run_cmd "$R_BIN" "$s" --config "$TMP_CFG";;
    go)       s="$(resolve_script go       fig3_go_enrichment.R)"; run_cmd "$R_BIN" "$s" --config "$TMP_CFG";;
    scenic1)  s="$(resolve_script scenic1  fig3_scenic_rds_to_csv.R)"; run_cmd "$R_BIN" "$s" --config "$TMP_CFG";;
    scenic2)  s="$(resolve_script scenic2  fig3_scenic_pyscenic.py)"; run_cmd "$PYTHON_BIN" "$s" --config "$TMP_CFG";;
    scenic_downstream) s="$(resolve_script scenic_downstream fig3_scenic_downstream.R)"; run_cmd "$R_BIN" "$s" --config "$TMP_CFG";;
    rds2h5ad) s="$(resolve_script rds2h5ad fig3_rds_to_h5ad.R)"; run_cmd "$R_BIN" "$s" --config "$TMP_CFG";;
    sctour)   s="$(resolve_script sctour   fig3_sctour_tnode.py)"; run_cmd "$PYTHON_BIN" "$s" --config "$TMP_CFG";;
    *) echo "ERROR: unknown step '$1'"; usage; exit 2;;
  esac
}

if [[ "$ONLY" == "all" || -z "$ONLY" ]]; then
  run_one monocle3; run_one length; run_one go; run_one scenic1; run_one scenic2; run_one scenic_downstream; run_one rds2h5ad; run_one sctour
else
  IFS=',' read -r -a STEPS <<< "$ONLY"; for st in "${STEPS[@]}"; do run_one "$st"; done
fi

rm -f "$TMP_CFG"
echo "[OK] Fig3 done."
