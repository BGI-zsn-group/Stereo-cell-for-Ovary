#!/usr/bin/env python3
"""
Fig.3: scTour (sctour) TNODE pseudotime & latent embedding

Recommended:
  python Fig3/fig3_sctour_tnode.py --config Fig3/configs/fig3_monocle3.yaml

This step is independent from the Monocle3 workflow, but we keep it in the same
Fig3 config for reproducibility.

Inputs (YAML, with defaults):
  sctour_input_h5ad: path to input .h5ad
Outputs (YAML, with defaults under out_dir):
  - sctour_output_h5ad: updated .h5ad containing:
      obs['ptime'], obsm['X_TNODE'], obsm['X_VF']
  - sctour_ptime_csv: per-cell ptime table (sorted)

Notes:
  - Requires: sctour, scanpy, numpy, pandas, pyyaml, scipy
  - If your adata.X is sparse, this script keeps it sparse but casts dtype to float32.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

try:
    import yaml  # PyYAML
except Exception:
    yaml = None

import scanpy as sc
import sctour as sct

try:
    import scipy.sparse as sp
except Exception:
    sp = None


def load_yaml(path: str) -> dict:
    if yaml is None:
        raise RuntimeError(
            "PyYAML is required for --config mode. Install with: pip install pyyaml"
        )
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def ensure_float32_X(adata):
    X = adata.X
    if sp is not None and sp.issparse(X):
        adata.X = X.astype(np.float32)
    else:
        adata.X = np.asarray(X).astype(np.float32)
    return adata


def main():
    ap = argparse.ArgumentParser(description="Fig3 scTour TNODE pseudotime & latents")
    ap.add_argument("--config", required=True, help="YAML config (reuse Fig3 config)")
    args = ap.parse_args()

    cfg = load_yaml(args.config)

    out_dir = Path(cfg.get("out_dir", "results/Fig3"))
    out_dir.mkdir(parents=True, exist_ok=True)

    in_h5ad = Path(cfg.get("sctour_input_h5ad", "oocyte_some_c1_filter_1211.h5ad"))
    if not in_h5ad.exists():
        raise FileNotFoundError(f"Missing sctour_input_h5ad: {in_h5ad}")

    sctour_out_dir = Path(cfg.get("sctour_out_dir", str(out_dir / "sctour")))
    sctour_out_dir.mkdir(parents=True, exist_ok=True)

    out_h5ad = Path(cfg.get("sctour_output_h5ad", str(sctour_out_dir / "oocyte_sctour_tnode.h5ad")))
    out_ptime_csv = Path(cfg.get("sctour_ptime_csv", str(sctour_out_dir / "ptime.csv")))

    # Params (match your snippet defaults)
    use_gpu = bool(cfg.get("sctour_use_gpu", False))
    loss_mode = str(cfg.get("sctour_loss_mode", "mse"))
    alpha_recon_lec = float(cfg.get("sctour_alpha_recon_lec", 0.5))
    alpha_recon_lode = float(cfg.get("sctour_alpha_recon_lode", 0.5))
    hvg_n_top = int(cfg.get("sctour_hvg_n_top_genes", 1000))
    alpha_z = float(cfg.get("sctour_alpha_z", 0.5))
    alpha_predz = float(cfg.get("sctour_alpha_predz", 0.5))

    print(f"[fig3-sctour] Reading: {in_h5ad}")
    adata = sc.read_h5ad(str(in_h5ad))

    # Ensure dtype
    adata = ensure_float32_X(adata)

    # QC + HVG
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=hvg_n_top)

    # Train
    print(f"[fig3-sctour] Training TNODE (use_gpu={use_gpu}, loss_mode={loss_mode}) ...")
    tnode = sct.train.Trainer(
        adata,
        use_gpu=use_gpu,
        loss_mode=loss_mode,
        alpha_recon_lec=alpha_recon_lec,
        alpha_recon_lode=alpha_recon_lode,
    )
    tnode.train()

    # Pseudotime + latents
    adata.obs["ptime"] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=alpha_z, alpha_predz=alpha_predz)
    adata.obsm["X_TNODE"] = mix_zs
    adata.obsm["X_VF"] = tnode.get_vector_field(adata.obs["ptime"].values, adata.obsm["X_TNODE"])

    # Sort by ptime
    order = np.argsort(adata.obs["ptime"].values)
    adata = adata[order, :].copy()

    print(f"[fig3-sctour] Writing: {out_h5ad}")
    adata.write_h5ad(str(out_h5ad))

    # Save ptime table
    ptime_df = pd.DataFrame({"cell": adata.obs_names, "ptime": adata.obs["ptime"].values})
    ptime_df.to_csv(out_ptime_csv, index=False)
    print(f"[fig3-sctour] Saved ptime: {out_ptime_csv} (n={len(ptime_df)})")

    # Save minimal params used
    params_txt = sctour_out_dir / "params_used_fig3_sctour.txt"
    params_txt.write_text(
        "\n".join(
            [
                f"sctour_input_h5ad: {in_h5ad}",
                f"sctour_output_h5ad: {out_h5ad}",
                f"sctour_ptime_csv: {out_ptime_csv}",
                f"sctour_use_gpu: {use_gpu}",
                f"sctour_loss_mode: {loss_mode}",
                f"sctour_alpha_recon_lec: {alpha_recon_lec}",
                f"sctour_alpha_recon_lode: {alpha_recon_lode}",
                f"sctour_hvg_n_top_genes: {hvg_n_top}",
                f"sctour_alpha_z: {alpha_z}",
                f"sctour_alpha_predz: {alpha_predz}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"[fig3-sctour] Saved params: {params_txt}")

    print("[fig3-sctour] Done.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[fig3-sctour] ERROR: {e}", file=sys.stderr)
        raise
