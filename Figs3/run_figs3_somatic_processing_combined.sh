\
#!/usr/bin/env bash
set -euo pipefail

usage(){
  cat <<'USAGE'
Usage:
  bash Figs3/run_figs3_somatic_processing_combined.sh [--only qc|somatic_processing|all]
                                                      [-o OUTROOT]
                                                      [--set key=value]...
                                                      [--config PATH | --config-dir DIR]
                                                      [--scripts-dir DIR]
                                                      [--dry-run]

Behavior:
  - default is `all`: run QC first, then somatic processing
  - `qc`: run only QC
  - `somatic_processing`: run only downstream processing

Examples:
  # full pipeline (recommended)
  bash Figs3/run_figs3_somatic_processing_combined.sh

  # only QC
  bash Figs3/run_figs3_somatic_processing_combined.sh --only qc

  # only downstream processing
  bash Figs3/run_figs3_somatic_processing_combined.sh --only somatic_processing

  # redirect outputs under a new root
  bash Figs3/run_figs3_somatic_processing_combined.sh -o results/Figs3_alt
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

declare -a OVERRIDES=()
OUTROOT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY="${2:-}"; shift 2;;
    --config) CONFIG="${2:-}"; shift 2;;
    --config-dir) CONFIG_DIR="${2:-}"; shift 2;;
    --scripts-dir) SCRIPTS_DIR="${2:-}"; shift 2;;
    --set) OVERRIDES+=("${2:-}"); shift 2;;
    -o|--out) OUTROOT="${2:-}"; shift 2;;
    --dry-run) DRY_RUN=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

# If user sets a root output dir, sync both modules:
#   qc -> <OUTROOT>/qc
#   processing -> <OUTROOT>/somatic_processing
# and wire processing.input_rds to qc_out_rds automatically
USER_PROC_INPUT_OVERRIDE=""
if [[ ${#OVERRIDES[@]} -gt 0 ]]; then
  for _ov in "${OVERRIDES[@]}"; do
    case "${_ov%%=*}" in
      figs3_somatic_processing.input_rds|input_rds) USER_PROC_INPUT_OVERRIDE="1" ;;
    esac
  done
fi

if [[ -n "$OUTROOT" ]]; then
  OVERRIDES+=(
    "figs3_qc.qc_out_dir=$OUTROOT/qc"
    "figs3_qc.qc_out_rds=$OUTROOT/qc/somatic.qc_processed.rds"
    "figs3_somatic_processing.out_dir=$OUTROOT/somatic_processing"
    "figs3_somatic_processing.out_round1_rds=$OUTROOT/somatic_processing/somatic.round1.rds"
    "figs3_somatic_processing.out_final_rds=$OUTROOT/somatic_processing/somatic.final.rds"
  )

  if [[ "$ONLY" == "all" || -z "$ONLY" ]]; then
    OVERRIDES+=("figs3_somatic_processing.input_rds=$OUTROOT/qc/somatic.qc_processed.rds")
  elif [[ "$ONLY" == "somatic_processing" && -z "$USER_PROC_INPUT_OVERRIDE" && -f "$OUTROOT/qc/somatic.qc_processed.rds" ]]; then
    OVERRIDES+=("figs3_somatic_processing.input_rds=$OUTROOT/qc/somatic.qc_processed.rds")
  fi
fi

find_config(){
  if [[ -n "$CONFIG" ]]; then echo "$CONFIG"; return; fi
  for name in figs3_combined.yaml figs3_combined_v2.yaml; do
    local p="$CONFIG_DIR/$name"
    [[ -f "$p" ]] && { echo "$p"; return; }
  done
  echo "ERROR: cannot find figs3_combined.yaml under $CONFIG_DIR" >&2
  exit 1
}

extract_module(){
  local combined="$1" mkey="$2" out_cfg="$3"; shift 3

  # Filter empty overrides
  local -a clean=()
  local ov
  for ov in "$@"; do
    [[ -n "${ov:-}" ]] && clean+=("$ov")
  done

  "$PYTHON_BIN" - <<'PY' "$combined" "$mkey" "$out_cfg" "${clean[@]}"
import sys
try:
  import yaml
except Exception:
  raise SystemExit('ERROR: missing dependency pyyaml. Install: pip install pyyaml')

combined, mkey, outp = sys.argv[1], sys.argv[2], sys.argv[3]
ovs = sys.argv[4:]
cfg = yaml.safe_load(open(combined,'r',encoding='utf-8')) or {}
mods = cfg.get('modules') or {}
if mkey not in mods:
  raise SystemExit(f"ERROR: module '{mkey}' not found under modules in {combined}")
mod = mods[mkey]

def parse_value(val):
  if val == "":
    return ""
  try:
    return yaml.safe_load(val)
  except Exception:
    return val

def set_path(root, path, val):
  ks = [k for k in path.split('.') if k]
  if not ks:
    raise SystemExit("ERROR: empty override key")

  cur = root
  parent = None
  parent_key = None

  for i, k in enumerate(ks[:-1]):
    nxt = ks[i + 1]
    want_list = nxt.isdigit()

    if k.isdigit():
      idx = int(k)
      if not isinstance(cur, list):
        new_list = []
        if parent is None:
          raise SystemExit(f"ERROR: path '{path}' cannot start with list index")
        parent[parent_key] = new_list
        cur = new_list
      while len(cur) <= idx:
        cur.append(None)
      if cur[idx] is None:
        cur[idx] = [] if want_list else {}
      parent, parent_key, cur = cur, idx, cur[idx]
    else:
      if not isinstance(cur, dict):
        raise SystemExit(f"ERROR: override path '{path}' conflicts with existing non-dict node at '{k}'")
      if k not in cur or cur[k] is None:
        cur[k] = [] if want_list else {}
      parent, parent_key, cur = cur, k, cur[k]

  last = ks[-1]
  value = parse_value(val)
  if last.isdigit():
    idx = int(last)
    if not isinstance(cur, list):
      if parent is None:
        raise SystemExit(f"ERROR: path '{path}' cannot end in list index at root")
      parent[parent_key] = []
      cur = parent[parent_key]
    while len(cur) <= idx:
      cur.append(None)
    cur[idx] = value
  else:
    if not isinstance(cur, dict):
      raise SystemExit(f"ERROR: override path '{path}' conflicts with existing non-dict leaf parent")
    cur[last] = value

for ov in ovs:
  if not ov or '=' not in ov:
    raise SystemExit(f"ERROR: override must be key=value (got: {ov})")
  k,v=ov.split('=',1)
  k=k.strip()
  if k.startswith(mkey+'.'):
    k=k[len(mkey)+1:]
  set_path(mod,k,v)

yaml.safe_dump(mod, open(outp,'w',encoding='utf-8'), sort_keys=False, allow_unicode=True)
PY
}

run_cmd(){
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '+'; printf ' %q' "$@"; printf '\n'
  else
    "$@"
  fi
}

COMBINED="$(find_config)"

run_qc(){
  local tmp; tmp="$(mktemp -t figs3_qc_XXXX.yaml)"
  extract_module "$COMBINED" "figs3_qc" "$tmp" "${OVERRIDES[@]}"
  run_cmd "$R_BIN" "$SCRIPTS_DIR/figs3_somatic_qc.R" --config "$tmp"
  rm -f "$tmp"
}

run_proc(){
  local tmp; tmp="$(mktemp -t figs3_proc_XXXX.yaml)"
  extract_module "$COMBINED" "figs3_somatic_processing" "$tmp" "${OVERRIDES[@]}"
  run_cmd "$R_BIN" "$SCRIPTS_DIR/figs3_somatic_processing.R" --config "$tmp"
  rm -f "$tmp"
}

case "$ONLY" in
  all|"") run_qc; run_proc;;
  qc) run_qc;;
  somatic_processing) run_proc;;
  *) echo "ERROR: unsupported --only '$ONLY'"; usage; exit 2;;
esac

echo "[OK] Figs3 finished."
