#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fig4: cell2location training (single-cell reference + spatial mapping)

This script keeps the main (reproducible) workflow:
  1) Load single-cell reference AnnData (.h5ad) with obs columns:
       - celltype (labels_key)
       - orig.ident (batch_key; optional)
     Filter genes (cell2location.utils.filtering.filter_genes)
     Train RegressionModel -> export posterior -> infer average expression per celltype (inf_aver)

  2) Load spatial AnnData (.h5ad) with raw counts in .X
     Intersect genes with inf_aver -> train Cell2location -> export posterior abundances
     Create a simple predicted label per spot (argmax of q05 abundances)
     Save trained models + processed AnnData + key tables

Usage:
  python Fig4/fig4_cell2location_train.py --config Fig4/configs/fig4_cell2location_train.yaml

Notes:
- Plotting is optional and disabled by default (for headless reproducibility).
- Paths in YAML are intended to be relative to the working directory (repo root recommended).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

try:
    import yaml
except ImportError as e:
    raise SystemExit("ERROR: Please install pyyaml (pip install pyyaml).") from e


def _now() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _log(msg: str, fh=None):
    line = f"[{_now()}] {msg}"
    print(line)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


def _read_yaml(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return cfg or {}


def _seed_everything(seed: int):
    try:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    except Exception:
        pass
    np.random.seed(seed)
    try:
        import scvi
        scvi.settings.seed = seed
    except Exception:
        pass


def _check_counts_integrality(adata: sc.AnnData, name: str, log_fh):
    from scipy.sparse import issparse
    X = adata.X
    if issparse(X):
        non_int = int(np.sum(X.data != np.floor(X.data)))
        dtype = str(X.dtype)
    else:
        non_int = int(np.sum(X != np.floor(X)))
        dtype = str(X.dtype)
    _log(f"[check] {name}.X dtype={dtype} non_integer_values={non_int}", log_fh)
    return non_int


def _infer_inf_aver(adata_sc: sc.AnnData) -> pd.DataFrame:
    # robust extraction after export_posterior
    if "mod" not in adata_sc.uns or "factor_names" not in adata_sc.uns["mod"]:
        raise RuntimeError("Missing adata_sc.uns['mod']['factor_names'] after export_posterior.")
    factor_names = list(adata_sc.uns["mod"]["factor_names"])

    if "means_per_cluster_mu_fg" in adata_sc.varm.keys():
        cols = [f"means_per_cluster_mu_fg_{i}" for i in factor_names]
        inf_aver = adata_sc.varm["means_per_cluster_mu_fg"][cols].copy()
    else:
        cols = [f"means_per_cluster_mu_fg_{i}" for i in factor_names]
        missing = [c for c in cols if c not in adata_sc.var.columns]
        if missing:
            raise RuntimeError(f"Cannot find means_per_cluster_mu_fg in varm or var columns. Missing: {missing[:5]}")
        inf_aver = adata_sc.var[cols].copy()

    inf_aver.columns = factor_names
    inf_aver.index = adata_sc.var_names.astype(str)
    return inf_aver


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="YAML config path")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    cfg = _read_yaml(cfg_path)

    out_dir = Path(cfg.get("out_dir", "results/cell2location_train"))
    _ensure_dir(out_dir)
    log_fh = open(out_dir / "cell2location_train.log", "w", encoding="utf-8")

    seed = int(cfg.get("seed", 0))
    _seed_everything(seed)
    _log(f"[start] config={cfg_path}", log_fh)
    _log(f"[seed] {seed}", log_fh)

    # -------------------
    # 1) Single-cell reference
    # -------------------
    sc_cfg = cfg.get("single_cell", {}) or {}
    sc_h5ad = Path(sc_cfg.get("input_h5ad", "results/cell2location_sc_ref/somatic_1226.h5ad"))
    labels_key = str(sc_cfg.get("labels_key", "celltype"))
    batch_key = sc_cfg.get("batch_key", "orig.ident")
    drop_mt = bool(sc_cfg.get("drop_mt", True))

    if not sc_h5ad.exists():
        raise SystemExit(f"ERROR: missing single-cell h5ad: {sc_h5ad}")

    adata_sc = sc.read_h5ad(str(sc_h5ad))
    _log(f"[sc] loaded: {sc_h5ad} n_obs={adata_sc.n_obs} n_vars={adata_sc.n_vars}", log_fh)

    if labels_key not in adata_sc.obs.columns:
        raise SystemExit(f"ERROR: labels_key '{labels_key}' not found in adata_sc.obs.")
    if batch_key is not None and batch_key != "" and batch_key not in adata_sc.obs.columns:
        _log(f"[sc] WARNING: batch_key '{batch_key}' not in obs; will set batch_key=None.", log_fh)
        batch_key = None

    _check_counts_integrality(adata_sc, "adata_sc", log_fh)

    if drop_mt:
        mt_mask = adata_sc.var_names.str.upper().str.startswith("MT-")
        n_mt = int(mt_mask.sum())
        _log(f"[sc] drop_mt=True -> removing {n_mt} MT- genes", log_fh)
        if n_mt > 0:
            adata_sc = adata_sc[:, ~mt_mask].copy()

    # Filter genes (cell2location helper)
    from cell2location.utils.filtering import filter_genes
    fg = sc_cfg.get("filter_genes", {}) or {}
    selected = filter_genes(
        adata_sc,
        cell_count_cutoff=int(fg.get("cell_count_cutoff", 5)),
        cell_percentage_cutoff2=float(fg.get("cell_percentage_cutoff2", 0.03)),
        nonz_mean_cutoff=float(fg.get("nonz_mean_cutoff", 1.12)),
    )
    adata_sc = adata_sc[:, selected].copy()
    _log(f"[sc] after filter_genes: n_vars={adata_sc.n_vars}", log_fh)

    # Train RegressionModel
    import cell2location
    from cell2location.models import RegressionModel

    RegressionModel.setup_anndata(
        adata=adata_sc,
        batch_key=batch_key,
        labels_key=labels_key,
        categorical_covariate_keys=None,
    )
    mod = RegressionModel(adata_sc)

    train_cfg = sc_cfg.get("train", {}) or {}
    use_gpu = bool(train_cfg.get("use_gpu", True))
    accelerator = "gpu" if use_gpu else "cpu"
    max_epochs = int(train_cfg.get("max_epochs", 800))
    lr = float(train_cfg.get("lr", 1e-3))
    train_size = float(train_cfg.get("train_size", 1.0))

    _log(f"[sc][train] accelerator={accelerator} max_epochs={max_epochs} lr={lr} train_size={train_size}", log_fh)
    mod.train(max_epochs=max_epochs, lr=lr, train_size=train_size, accelerator=accelerator)

    # Export posterior (for signatures)
    exp_cfg = sc_cfg.get("export_posterior", {}) or {}
    num_samples = int(exp_cfg.get("num_samples", 1000))
    batch_size = int(exp_cfg.get("batch_size", 2500))
    adata_sc = mod.export_posterior(
        adata_sc,
        sample_kwargs={"num_samples": num_samples, "batch_size": batch_size},
    )

    # Save reference signatures
    ref_dir = out_dir / "reference_signatures"
    _ensure_dir(ref_dir)
    mod.save(str(ref_dir), overwrite=True)
    ref_h5ad = ref_dir / "reference_trained.h5ad"
    adata_sc.write_h5ad(str(ref_h5ad))
    _log(f"[sc][save] model_dir={ref_dir} ref_h5ad={ref_h5ad}", log_fh)

    # Quantile export (QC + robust inf_aver)
    adata_sc = mod.export_posterior(
        adata_sc,
        use_quantiles=True,
        add_to_varm=["q05", "q50", "q95", "q0001"],
        sample_kwargs={"batch_size": batch_size},
    )

    inf_aver = _infer_inf_aver(adata_sc)
    inf_aver_csv = out_dir / "inf_aver.csv"
    inf_aver.to_csv(inf_aver_csv)
    _log(f"[sc][save] inf_aver_csv={inf_aver_csv} shape={inf_aver.shape}", log_fh)

    # -------------------
    # 2) Spatial mapping
    # -------------------
    sp_cfg = cfg.get("spatial", {}) or {}
    sp_h5ad = Path(sp_cfg.get("input_h5ad", "results/cell2location_spatial/B04372C211.cell2location_spatial_counts.h5ad"))
    sp_batch_key = sp_cfg.get("batch_key", "orig.ident")
    if sp_batch_key is not None and sp_batch_key != "" and not sp_h5ad.exists():
        pass

    if not sp_h5ad.exists():
        raise SystemExit(f"ERROR: missing spatial h5ad: {sp_h5ad}")

    adata_sp = sc.read_h5ad(str(sp_h5ad))
    _log(f"[sp] loaded: {sp_h5ad} n_obs={adata_sp.n_obs} n_vars={adata_sp.n_vars}", log_fh)
    _check_counts_integrality(adata_sp, "adata_sp", log_fh)

    # Intersect genes
    intersect = np.intersect1d(adata_sp.var_names.astype(str), inf_aver.index.astype(str))
    _log(f"[sp] intersect genes: {len(intersect)}", log_fh)
    adata_sp = adata_sp[:, intersect].copy()
    inf_aver_sp = inf_aver.loc[intersect, :].copy()

    inf_aver_sp_csv = out_dir / "inf_aver_sp.csv"
    inf_aver_sp.to_csv(inf_aver_sp_csv)
    _log(f"[sp][save] inf_aver_sp_csv={inf_aver_sp_csv} shape={inf_aver_sp.shape}", log_fh)

    # Setup and train Cell2location
    from cell2location.models import Cell2location

    if sp_batch_key is not None and sp_batch_key != "" and sp_batch_key not in adata_sp.obs.columns:
        _log(f"[sp] WARNING: batch_key '{sp_batch_key}' not in obs; set batch_key=None.", log_fh)
        sp_batch_key = None

    Cell2location.setup_anndata(adata=adata_sp, batch_key=sp_batch_key)

    map_cfg = sp_cfg.get("map", {}) or {}
    N_cells_per_location = float(map_cfg.get("N_cells_per_location", 1.0))
    detection_alpha = float(map_cfg.get("detection_alpha", 20.0))

    mod_tr = Cell2location(
        adata_sp,
        cell_state_df=inf_aver_sp,
        N_cells_per_location=N_cells_per_location,
        detection_alpha=detection_alpha,
    )

    tr_cfg = sp_cfg.get("train", {}) or {}
    max_epochs_sp = int(tr_cfg.get("max_epochs", 10000))
    train_size_sp = float(tr_cfg.get("train_size", 1.0))
    batch_size_sp = tr_cfg.get("batch_size", None)  # None = use full batch in cell2location
    _log(f"[sp][train] max_epochs={max_epochs_sp} train_size={train_size_sp} batch_size={batch_size_sp}", log_fh)

    mod_tr.train(max_epochs=max_epochs_sp, batch_size=batch_size_sp, train_size=train_size_sp)

    # Export posterior abundances
    exp2 = sp_cfg.get("export_posterior", {}) or {}
    num_samples_sp = int(exp2.get("num_samples", 1000))
    batch_size_export = exp2.get("batch_size", None)
    if batch_size_export is None:
        batch_size_export = int(adata_sp.n_obs)
    adata_sp = mod_tr.export_posterior(
        adata_sp,
        sample_kwargs={"num_samples": num_samples_sp, "batch_size": int(batch_size_export)},
    )

    # Save spatial model + anndata
    sp_model_dir = out_dir / "spatial_model"
    _ensure_dir(sp_model_dir)
    mod_tr.save(str(sp_model_dir), overwrite=True)

    out_sp_h5ad = out_dir / "spatial_trained.h5ad"
    adata_sp.write_h5ad(str(out_sp_h5ad))
    _log(f"[sp][save] model_dir={sp_model_dir} out_h5ad={out_sp_h5ad}", log_fh)

    # predicted label (argmax of q05)
    if "q05_cell_abundance_w_sf" in adata_sp.obsm:
        celltypes = list(adata_sp.uns["mod"]["factor_names"])
        q05 = adata_sp.obsm["q05_cell_abundance_w_sf"]
        if isinstance(q05, pd.DataFrame):
            mat = q05.values
        else:
            mat = np.asarray(q05)
        pred_idx = np.argmax(mat, axis=1)
        pred_labels = [celltypes[i] for i in pred_idx]
        adata_sp.obs["predicted_cell_type"] = pred_labels
        pred_csv = out_dir / "predicted_cell_type.csv"
        pd.DataFrame({"spot": adata_sp.obs_names, "predicted_cell_type": pred_labels}).to_csv(pred_csv, index=False)
        _log(f"[sp][save] predicted_cell_type.csv={pred_csv}", log_fh)
    else:
        _log("[sp] WARNING: missing obsm['q05_cell_abundance_w_sf']; skip predicted_cell_type.", log_fh)

    # Save params used
    cfg_out = dict(cfg)
    cfg_out["resolved"] = {
        "config": str(cfg_path),
        "out_dir": str(out_dir),
        "reference_model_dir": str(ref_dir),
        "reference_h5ad": str(ref_h5ad),
        "inf_aver_csv": str(inf_aver_csv),
        "spatial_model_dir": str(sp_model_dir),
        "spatial_h5ad": str(out_sp_h5ad),
    }
    with open(out_dir / "params_used_fig4_cell2location_train.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg_out, f, sort_keys=False, allow_unicode=True)

    _log("[done] finished", log_fh)
    log_fh.close()


if __name__ == "__main__":
    main()
