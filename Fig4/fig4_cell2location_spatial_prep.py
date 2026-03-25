#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig4 spatial preprocessing for cell2location.
Notebook-aligned, Python 3.8 compatible.

Key behavior:
1) read one adjusted .h5ad
2) aggregate duplicated genes by `gene_name_col`
3) QC filter cells/genes
4) create internal full-count snapshot BEFORE normalize/log/HVG
5) do processing / leiden-based subset selection
6) build final cell2location input from the snapshot (full genes + raw counts)
"""
from __future__ import annotations

import argparse
import datetime as _dt
from pathlib import Path
from typing import List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

try:
    import yaml
except ImportError as e:
    raise SystemExit("ERROR: Please install pyyaml (pip install pyyaml).") from e


def _now() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


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


def _ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def _as_list(x) -> List[str]:
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        return [str(i) for i in x]
    return [str(x)]


def _read_id_list(path: Path) -> List[str]:
    with open(path, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def _aggregate_by_gene_name(adata: ad.AnnData, gene_name_col: str, drop_gene_names: Optional[List[str]], log_fh=None) -> ad.AnnData:
    if gene_name_col not in adata.var.columns:
        raise SystemExit(f"ERROR: gene_name_col '{gene_name_col}' not found in adata.var")

    gene_names = pd.Series(adata.var[gene_name_col]).astype(str).fillna("")
    valid = gene_names != ""
    if drop_gene_names:
        drop_set = set(str(x) for x in drop_gene_names)
        valid &= ~gene_names.isin(drop_set)
        _log(f"[gene] dropping {int((~valid).sum())} genes by name filter: {sorted(list(drop_set))}", log_fh)

    adata = adata[:, valid.to_numpy()].copy()
    gene_names = pd.Series(adata.var[gene_name_col]).astype(str).fillna("")

    X = adata.X.tocsc() if sparse.issparse(adata.X) else sparse.csc_matrix(adata.X)
    unique_genes, inverse_idx = np.unique(gene_names.values, return_inverse=True)
    _log(f"[gene] aggregating duplicated genes: {adata.n_vars} -> {len(unique_genes)}", log_fh)

    X_agg = sparse.lil_matrix((X.shape[0], len(unique_genes)), dtype=X.dtype)
    for i in range(len(unique_genes)):
        cols = np.where(inverse_idx == i)[0]
        if len(cols) == 1:
            X_agg[:, i] = X[:, cols[0]]
        else:
            X_agg[:, i] = X[:, cols].sum(axis=1)
    X_agg = X_agg.tocsr()

    new_var = pd.DataFrame(index=pd.Index(unique_genes.astype(str), name=None))
    adata_out = ad.AnnData(
        X=X_agg,
        obs=adata.obs.copy(),
        var=new_var,
        obsm={k: v.copy() if hasattr(v, "copy") else v for k, v in adata.obsm.items()},
        uns=adata.uns.copy(),
    )
    adata_out.var_names = adata_out.var.index.astype(str)
    adata_out.var = pd.DataFrame(index=adata_out.var_names)
    return adata_out


def _ensure_counts_layer(adata: ad.AnnData):
    X = adata.X.copy() if hasattr(adata.X, "copy") else np.array(adata.X, copy=True)
    adata.layers["counts"] = X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = _read_yaml(Path(args.config))
    in_h5ad = Path(cfg.get("input_h5ad", "data/Fig4/spatial/B04372C211.adjusted.cellbin.h5ad"))
    if not in_h5ad.exists():
        raise SystemExit(f"ERROR: missing input_h5ad: {in_h5ad}")

    out_dir = Path(cfg.get("out_dir", "results/Fig4/cell2location_spatial"))
    _ensure_dir(out_dir)
    prefix = str(cfg.get("prefix", in_h5ad.name.replace(".adjusted.cellbin.h5ad", "").replace(".h5ad", "")))
    log_fh = open(out_dir / f"{prefix}.fig4_cell2location_spatial.log", "w", encoding="utf-8")

    _log(f"[start] config={args.config}", log_fh)
    _log(f"[io] input_h5ad={in_h5ad}", log_fh)
    _log(f"[io] out_dir={out_dir}", log_fh)
    _log(f"[io] prefix={prefix}", log_fh)

    adata = sc.read_h5ad(str(in_h5ad))
    _log(f"[load] adata: n_obs={adata.n_obs}, n_vars={adata.n_vars}", log_fh)

    spatial_key = str(cfg.get("spatial_obsm_key", "spatial"))
    if spatial_key not in adata.obsm.keys():
        _log(f"[warn] obsm['{spatial_key}'] not found; continue without spatial coordinates", log_fh)

    adata = _aggregate_by_gene_name(
        adata,
        gene_name_col=str(cfg.get("gene_name_col", "real_gene_name")),
        drop_gene_names=_as_list(cfg.get("drop_gene_names", ["a"])),
        log_fh=log_fh,
    )
    _log(f"[gene] after aggregation: n_vars={adata.n_vars}", log_fh)

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.contains("^HB[AB]", regex=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    min_genes = int(cfg.get("min_genes", 100))
    min_cells = int(cfg.get("min_cells", 3))
    _log(f"[qc] filter_cells min_genes={min_genes}; filter_genes min_cells={min_cells}", log_fh)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Notebook-aligned snapshot after QC filtering but before normalize/log/HVG
    adata_orig = adata.copy()
    _ensure_counts_layer(adata_orig)
    orig_qc_h5ad = Path(cfg.get("out_orig_qc_h5ad", out_dir / f"{prefix}.orig_qc_counts.h5ad"))
    adata_orig.write_h5ad(str(orig_qc_h5ad))
    _log(f"[save] orig_qc_h5ad={orig_qc_h5ad}", log_fh)

    do_processing = bool(cfg.get("do_processing", True))
    if do_processing:
        sc.pp.normalize_total(adata, target_sum=float(cfg.get("normalize_target_sum", 1e4)))
        sc.pp.log1p(adata)
        adata.raw = adata.copy()

        sc.pp.highly_variable_genes(
            adata,
            flavor=str(cfg.get("hvg_flavor", "seurat_v3")),
            n_top_genes=int(cfg.get("hvg_n_top_genes", 2000)),
        )
        if "highly_variable" in adata.var.columns:
            adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()

        if bool(cfg.get("do_cell_cycle", True)):
            s_genes = [g for g in _as_list(cfg.get("s_genes")) if g in adata.var_names]
            g2m_genes = [g for g in _as_list(cfg.get("g2m_genes")) if g in adata.var_names]
            if len(s_genes) > 0 and len(g2m_genes) > 0:
                sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, copy=False)
                sc.pp.regress_out(adata, ["S_score", "G2M_score"])
                _log(f"[cell-cycle] scored/regressed with {len(s_genes)} S genes and {len(g2m_genes)} G2M genes", log_fh)

        sc.pp.scale(adata, max_value=float(cfg.get("scale_max_value", 10)))
        sc.tl.pca(adata, svd_solver=str(cfg.get("pca_solver", "arpack")))
        sc.pp.neighbors(adata, n_neighbors=int(cfg.get("n_neighbors", 20)), n_pcs=int(cfg.get("n_pcs", 30)))
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=float(cfg.get("leiden_resolution_10", 1.0)), key_added=str(cfg.get("leiden_key_10", "leiden_10")))
        sc.tl.leiden(adata, resolution=float(cfg.get("leiden_resolution_15", 1.5)), key_added=str(cfg.get("leiden_key_15", "leiden_15")))

    processed_h5ad = Path(cfg.get("out_processed_h5ad", out_dir / f"{prefix}.adjusted_processed.h5ad"))
    adata.write_h5ad(str(processed_h5ad))
    _log(f"[save] processed_h5ad={processed_h5ad}", log_fh)

    subset_cells = list(adata.obs_names)

    roi_cfg = cfg.get("roi", {}) or {}
    cell_ids_txt = roi_cfg.get("cell_ids_txt", None)
    if cell_ids_txt not in (None, ""):
        roi_ids = set(_read_id_list(Path(cell_ids_txt)))
        subset_cells = [c for c in subset_cells if c in roi_ids]
        _log(f"[roi] keep by cell_ids_txt -> {len(subset_cells)} cells", log_fh)

    subset_cfg = cfg.get("subset", {}) or {}
    leiden_key = str(subset_cfg.get("leiden_key", cfg.get("leiden_key_15", "leiden_15")))
    exclude_clusters = set(_as_list(subset_cfg.get("exclude_clusters", [])))
    if leiden_key in adata.obs.columns and len(exclude_clusters) > 0:
        subset_cells = [c for c in subset_cells if str(adata.obs.loc[c, leiden_key]) not in exclude_clusters]
        _log(f"[subset] leiden_key={leiden_key} exclude={sorted(exclude_clusters)} -> {len(subset_cells)} cells", log_fh)

    adata_proc_sub = adata[subset_cells, :].copy()
    subset_processed_h5ad = Path(cfg.get("out_subset_processed_h5ad", out_dir / f"{prefix}.adjusted_processed_subset.h5ad"))
    adata_proc_sub.write_h5ad(str(subset_processed_h5ad))
    _log(f"[save] subset_processed_h5ad={subset_processed_h5ad}", log_fh)

    # Final c2l input: take subset cells from the pre-normalization snapshot
    adata_counts = adata_orig[subset_cells, :].copy()
    _ensure_counts_layer(adata_counts)
    if "orig.ident" not in adata_counts.obs.columns:
        adata_counts.obs["orig.ident"] = prefix

    cell2loc_h5ad = Path(cfg.get("out_cell2location_h5ad", out_dir / f"{prefix}.cell2location_spatial_counts.h5ad"))
    adata_counts.write_h5ad(str(cell2loc_h5ad))
    _log(f"[save] cell2location_h5ad={cell2loc_h5ad} n_obs={adata_counts.n_obs} n_vars={adata_counts.n_vars}", log_fh)

    with open(out_dir / "params_used_fig4_cell2location_spatial.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)
    _log("[done] spatial preprocessing finished", log_fh)
    log_fh.close()


if __name__ == "__main__":
    main()
