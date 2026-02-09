#!/usr/bin/env python3
"""Fig.3: SCENIC (pySCENIC) pipeline (Step 2)

This script:
  1) Converts the SCENIC input CSV (cells x genes) to a .loom file
  2) Runs pySCENIC: grn -> ctx -> aucell

Recommended:
  python Fig3/fig3_scenic_pyscenic.py --config Fig3/configs/fig3_monocle3.yaml

Inputs (YAML; defaults shown):
  scenic_out_csv: "results/Fig3/scenic/oocyte_1211_withoutMII.csv"   # produced by Step 1
  scenic_tf_list: "allTFs_mm.txt"
  scenic_rankings_feather: "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
  scenic_motif_annotations: "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

Outputs (YAML; defaults under scenic_out_dir):
  scenic_loom: "results/Fig3/scenic/oocyte_1211_withoutMII.loom"
  scenic_adj_tsv: "results/Fig3/scenic/adj_oocyte_1211_withoutMII.tsv"
  scenic_reg_csv: "results/Fig3/scenic/reg_oocyte_1211_withoutMII.csv"
  scenic_result_loom: "results/Fig3/scenic/oocyte_1211_withoutMII_result.loom"

Notes:
  - Requires: scanpy, loompy, numpy, pyyaml
  - Requires pySCENIC CLI available as `pyscenic` in PATH (or set scenic_pyscenic_bin in YAML)
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import shlex

import numpy as np

try:
    import yaml
except Exception:
    yaml = None

import scanpy as sc
import loompy as lp


def load_yaml(path: str) -> dict:
    if yaml is None:
        raise RuntimeError("PyYAML is required. Install with: pip install pyyaml")
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def run_cmd(cmd: list[str]) -> None:
    # Print as a shell-escaped string for logs
    cmd_str = " ".join(shlex.quote(c) for c in cmd)
    print(f"[fig3-scenic2] $ {cmd_str}")
    subprocess.run(cmd, check=True)


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig3 pySCENIC step2: CSV->loom and run grn/ctx/aucell")
    ap.add_argument("--config", required=True, help="YAML config (reuse Fig3 config)")
    args = ap.parse_args()

    cfg = load_yaml(args.config)

    # Directories
    scenic_out_dir = Path(cfg.get("scenic_out_dir", "results/Fig3/scenic"))
    scenic_out_dir.mkdir(parents=True, exist_ok=True)

    # Inputs from step1
    csv_path = Path(cfg.get("scenic_out_csv", str(scenic_out_dir / "oocyte_1211_withoutMII.csv")))
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing scenic_out_csv: {csv_path}\n(Please run Step1: fig3_scenic_rds_to_csv.R first.)")

    tf_list = Path(cfg.get("scenic_tf_list", "allTFs_mm.txt"))
    rankings = Path(cfg.get("scenic_rankings_feather", "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
    motif_ann = Path(cfg.get("scenic_motif_annotations", "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"))

    for p, name in [(tf_list, "scenic_tf_list"), (rankings, "scenic_rankings_feather"), (motif_ann, "scenic_motif_annotations")]:
        if not p.exists():
            raise FileNotFoundError(f"Missing {name}: {p}")

    # Outputs
    loom_path = Path(cfg.get("scenic_loom", str(scenic_out_dir / (csv_path.stem + ".loom"))))
    adj_tsv = Path(cfg.get("scenic_adj_tsv", str(scenic_out_dir / f"adj_{csv_path.stem}.tsv")))
    reg_csv = Path(cfg.get("scenic_reg_csv", str(scenic_out_dir / f"reg_{csv_path.stem}.csv")))
    result_loom = Path(cfg.get("scenic_result_loom", str(scenic_out_dir / f"{csv_path.stem}_result.loom")))

    pyscenic_bin = str(cfg.get("scenic_pyscenic_bin", "pyscenic"))

    # Workers / mode
    grn_workers = int(cfg.get("scenic_grn_workers", 14))
    ctx_workers = int(cfg.get("scenic_ctx_workers", 30))
    auc_workers = int(cfg.get("scenic_aucell_workers", 16))
    ctx_mode = str(cfg.get("scenic_ctx_mode", "dask_multiprocessing"))
    ctx_mask_dropouts = bool(cfg.get("scenic_ctx_mask_dropouts", True))
    grn_method = str(cfg.get("scenic_grn_method", "grnboost2"))

    # 1) CSV -> loom
    print(f"[fig3-scenic2] Reading CSV: {csv_path}")
    x = sc.read_csv(str(csv_path))  # expects cells x genes with first col as index
    # Ensure dense array for loompy
    X = x.X
    try:
        import scipy.sparse as sp
        if sp.issparse(X):
            X = X.toarray()
    except Exception:
        pass

    # loom expects genes x cells matrix, plus row/col attrs
    row_attrs = {"Gene": np.array(x.var_names)}
    col_attrs = {"CellID": np.array(x.obs_names)}

    print(f"[fig3-scenic2] Writing loom: {loom_path}")
    loom_path.parent.mkdir(parents=True, exist_ok=True)
    # transpose: cells x genes -> genes x cells
    lp.create(str(loom_path), np.asarray(X).T, row_attrs, col_attrs)

    # 2) pySCENIC GRN
    run_cmd([
        pyscenic_bin, "grn",
        "--num_workers", str(grn_workers),
        "--output", str(adj_tsv),
        "--method", grn_method,
        str(loom_path),
        str(tf_list),
    ])

    # 3) pySCENIC ctx
    ctx_cmd = [
        pyscenic_bin, "ctx",
        str(adj_tsv),
        str(rankings),
        "--annotations_fname", str(motif_ann),
        "--expression_mtx_fname", str(loom_path),
        "--mode", ctx_mode,
        "--output", str(reg_csv),
        "--num_workers", str(ctx_workers),
    ]
    if ctx_mask_dropouts:
        ctx_cmd.append("--mask_dropouts")
    run_cmd(ctx_cmd)

    # 4) pySCENIC aucell
    run_cmd([
        pyscenic_bin, "aucell",
        str(loom_path),
        str(reg_csv),
        "--output", str(result_loom),
        "--num_workers", str(auc_workers),
    ])

    # Minimal provenance
    params_txt = scenic_out_dir / "params_used_fig3_scenic2.txt"
    params_txt.write_text(
        "\n".join([
            f"scenic_out_csv: {csv_path}",
            f"scenic_loom: {loom_path}",
            f"scenic_tf_list: {tf_list}",
            f"scenic_rankings_feather: {rankings}",
            f"scenic_motif_annotations: {motif_ann}",
            f"scenic_adj_tsv: {adj_tsv}",
            f"scenic_reg_csv: {reg_csv}",
            f"scenic_result_loom: {result_loom}",
            f"scenic_pyscenic_bin: {pyscenic_bin}",
            f"scenic_grn_workers: {grn_workers}",
            f"scenic_grn_method: {grn_method}",
            f"scenic_ctx_workers: {ctx_workers}",
            f"scenic_ctx_mode: {ctx_mode}",
            f"scenic_ctx_mask_dropouts: {ctx_mask_dropouts}",
            f"scenic_aucell_workers: {auc_workers}",
        ]) + "\n",
        encoding="utf-8",
    )
    print(f"[fig3-scenic2] Saved params: {params_txt}")
    print("[fig3-scenic2] Done.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[fig3-scenic2] ERROR: {e}", file=sys.stderr)
        raise
