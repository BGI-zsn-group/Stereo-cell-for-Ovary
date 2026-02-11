#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fig4: Cell2location spatial transcriptomics preprocessing (AnnData)

Goal
- Prepare a spatial AnnData object with:
  - raw counts in .X (and also stored in .layers["counts"])
  - unique gene symbols in .var_names (optional aggregation by real_gene_name)
  - spatial coordinates in .obsm["spatial"]
- (Optional) create a processed object for QC/clustering and select a subset of spots/cells.
- (Optional) apply an ROI selection by a provided list of cell/spot IDs.

Usage
  python fig4_cell2location_spatial_prep.py --config Fig4/configs/fig4_cell2location_spatial_prep.yaml

Notes
- This is a reproducible, script-friendly version of the original notebook-like pipeline.
- Plotting and interactive lasso ROI are intentionally removed; ROI selection is driven by config.
"""

import argparse
import datetime as _dt
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

try:
    import yaml
except ImportError as e:
    raise SystemExit("ERROR: Please install pyyaml (pip install pyyaml).") from e


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _now() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _read_yaml(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return cfg or {}


def _log_setup(log_path: Path):
    _ensure_dir(log_path.parent)
    f = open(log_path, "w", encoding="utf-8")
    return f


def _log(msg: str, fh):
    line = f"[{_now()}] {msg}"
    print(line)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


def _aggregate_by_gene_name(
    adata: ad.AnnData,
    gene_name_col: str,
    drop_gene_names: list[str] | None,
    log_fh,
) -> ad.AnnData:
    """Aggregate duplicated genes by `gene_name_col` using sparse COO summation."""
    if gene_name_col not in adata.var.columns:
        _log(f"[gene] var does not contain '{gene_name_col}', skip aggregation.", log_fh)
        # Ensure var_names unique anyway
        adata.var_names_make_unique()
        return adata

    gene_names = adata.var[gene_name_col].astype(str).values
    if drop_gene_names:
        drop_set = set([str(x) for x in drop_gene_names])
        keep_mask = ~pd.isna(gene_names)
        keep_mask &= np.array([g not in drop_set for g in gene_names], dtype=bool)
        n_drop = int((~keep_mask).sum())
        if n_drop > 0:
            _log(f"[gene] dropping {n_drop} genes by name filter: {sorted(list(drop_set))[:10]}", log_fh)
        adata = adata[:, keep_mask].copy()
        gene_names = adata.var[gene_name_col].astype(str).values

    unique_genes, inverse_idx = np.unique(gene_names, return_inverse=True)
    if len(unique_genes) == adata.n_vars:
        _log("[gene] no duplicated gene names detected; set var_names to gene symbols.", log_fh)
        adata.var_names = pd.Index(unique_genes)
        adata.var = pd.DataFrame(index=adata.var_names)
        adata.var_names_make_unique()
        return adata

    _log(f"[gene] aggregating duplicated genes: {adata.n_vars} -> {len(unique_genes)}", log_fh)

    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    X = X.tocoo()
    new_col = inverse_idx[X.col]
    X_new = sparse.coo_matrix((X.data, (X.row, new_col)), shape=(adata.n_obs, len(unique_genes))).tocsr()

    new_var = pd.DataFrame(index=pd.Index(unique_genes, name="gene_symbol"))
    adata_new = ad.AnnData(
        X=X_new,
        obs=adata.obs.copy(),
        var=new_var,
        obsm=adata.obsm.copy(),
        uns=adata.uns.copy(),
    )
    # carry layers if exist and represent counts (optional)
    for layer_name, layer in getattr(adata, "layers", {}).items():
        try:
            L = layer
            if not sparse.issparse(L):
                L = sparse.csr_matrix(L)
            L = L.tocoo()
            L_new = sparse.coo_matrix((L.data, (L.row, inverse_idx[L.col])), shape=(adata.n_obs, len(unique_genes))).tocsr()
            adata_new.layers[layer_name] = L_new
        except Exception:
            # ignore non-matrix layers
            pass

    adata_new.var_names = adata_new.var.index.astype(str)
    adata_new.var_names_make_unique()
    return adata_new


def _add_qc_flags(adata: ad.AnnData):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB[AB]", regex=True)


def _read_id_list(path: Path) -> list[str]:
    ids = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            ids.append(s)
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="YAML config path")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    cfg = _read_yaml(cfg_path)

    # I/O
    input_h5ad = Path(cfg.get("input_h5ad", "B04372C211.adjusted.cellbin.h5ad"))
    out_dir = Path(cfg.get("out_dir", "results/cell2location_spatial"))
    prefix = str(cfg.get("prefix", input_h5ad.stem))

    _ensure_dir(out_dir)
    log_fh = _log_setup(out_dir / f"{prefix}.spatial_prep.log")

    _log(f"[start] config={cfg_path}", log_fh)
    _log(f"[io] input_h5ad={input_h5ad}", log_fh)
    _log(f"[io] out_dir={out_dir}", log_fh)
    _log(f"[io] prefix={prefix}", log_fh)

    if not input_h5ad.exists():
        raise SystemExit(f"ERROR: missing input_h5ad: {input_h5ad}")

    # Load
    adata = sc.read_h5ad(str(input_h5ad))
    _log(f"[load] adata: n_obs={adata.n_obs}, n_vars={adata.n_vars}", log_fh)

    # Ensure spatial coordinates exist
    spatial_key = cfg.get("spatial_obsm_key", "spatial")
    if spatial_key not in adata.obsm_keys():
        raise SystemExit(f"ERROR: missing adata.obsm['{spatial_key}'] (spatial coordinates).")

    # Aggregate duplicated genes by gene_name_col (optional)
    gene_name_col = cfg.get("gene_name_col", "real_gene_name")
    drop_gene_names = cfg.get("drop_gene_names", ["a"])  # default matches your notebook inspection
    if drop_gene_names is None:
        drop_gene_names = []
    adata = _aggregate_by_gene_name(adata, gene_name_col=gene_name_col, drop_gene_names=drop_gene_names, log_fh=log_fh)
    _log(f"[gene] after aggregation: n_vars={adata.n_vars}", log_fh)

    # QC + filtering on counts
    # Store raw counts
    adata.layers["counts"] = adata.X.copy() if sparse.issparse(adata.X) else np.array(adata.X, copy=True)

    _add_qc_flags(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=False)

    min_genes = int(cfg.get("min_genes", 100))
    min_cells = int(cfg.get("min_cells", 3))
    _log(f"[qc] filter_cells min_genes={min_genes}; filter_genes min_cells={min_cells}", log_fh)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Re-sync counts layer after filtering
    adata.layers["counts"] = adata.layers["counts"][adata.obs_names, adata.var_names]

    orig_qc_h5ad = out_dir / f"{prefix}.orig_qc_counts.h5ad"
    sc.write_h5ad(str(orig_qc_h5ad), adata)
    _log(f"[save] orig_qc_counts={orig_qc_h5ad}", log_fh)

    # Processed (for clustering / selection)
    do_processing = bool(cfg.get("do_processing", True))
    processed_h5ad = out_dir / f"{prefix}.processed.h5ad"
    subset_processed_h5ad = out_dir / f"{prefix}.subset_processed.h5ad"

    subset_cells = None

    if do_processing:
        adata_proc = adata.copy()

        # normalize/log + HVG
        target_sum = float(cfg.get("normalize_target_sum", 1e4))
        sc.pp.normalize_total(adata_proc, target_sum=target_sum)
        sc.pp.log1p(adata_proc)
        adata_proc.raw = adata_proc.copy()

        n_top_genes = int(cfg.get("hvg_n_top_genes", 2000))
        sc.pp.highly_variable_genes(adata_proc, flavor="seurat_v3", n_top_genes=n_top_genes)
        adata_proc = adata_proc[:, adata_proc.var["highly_variable"]].copy()

        # Cell cycle scoring (optional)
        do_cell_cycle = bool(cfg.get("do_cell_cycle", True))
        if do_cell_cycle:
            s_genes = cfg.get("s_genes", [])
            g2m_genes = cfg.get("g2m_genes", [])
            if len(s_genes) > 0 and len(g2m_genes) > 0:
                sc.tl.score_genes_cell_cycle(adata_proc, s_genes=s_genes, g2m_genes=g2m_genes, copy=False)
                sc.pp.regress_out(adata_proc, ["S_score", "G2M_score"])
            else:
                _log("[cellcycle] s_genes/g2m_genes empty; skip cell cycle regression.", log_fh)

        sc.pp.scale(adata_proc, max_value=float(cfg.get("scale_max_value", 10)))
        sc.tl.pca(adata_proc, svd_solver="arpack")

        n_neighbors = int(cfg.get("n_neighbors", 20))
        n_pcs = int(cfg.get("n_pcs", 30))
        sc.pp.neighbors(adata_proc, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(adata_proc)

        # Leiden clustering
        leiden_resolutions = cfg.get("leiden_resolutions", {"leiden_10": 1.0, "leiden_15": 1.5})
        for key, res in leiden_resolutions.items():
            sc.tl.leiden(adata_proc, resolution=float(res), key_added=str(key))
        sc.write_h5ad(str(processed_h5ad), adata_proc)
        _log(f"[save] processed={processed_h5ad}", log_fh)

        # subset by excluding clusters (optional)
        subset_cfg = cfg.get("subset", {}) or {}
        leiden_key = subset_cfg.get("leiden_key", "leiden_15")
        exclude_clusters = [str(x) for x in subset_cfg.get("exclude_clusters", [])]

        if exclude_clusters:
            if leiden_key not in adata_proc.obs.columns:
                raise SystemExit(f"ERROR: subset.leiden_key='{leiden_key}' not found in processed adata.obs.")
            keep_mask = ~adata_proc.obs[leiden_key].astype(str).isin(exclude_clusters)
            adata_proc_sub = adata_proc[keep_mask].copy()
            subset_cells = adata_proc_sub.obs_names.tolist()
            sc.write_h5ad(str(subset_processed_h5ad), adata_proc_sub)
            _log(f"[subset] by clusters: kept {len(subset_cells)} cells/spots; saved {subset_processed_h5ad}", log_fh)
        else:
            subset_cells = adata_proc.obs_names.tolist()
            _log(f"[subset] no exclude_clusters provided; keep all {len(subset_cells)} cells/spots", log_fh)
    else:
        subset_cells = adata.obs_names.tolist()
        _log(f"[proc] do_processing=False; skip processing and keep all {len(subset_cells)}", log_fh)

    # Build a counts AnnData (full genes) for cell2location input
    adata_counts = adata[subset_cells, :].copy()
    adata_counts.X = adata_counts.layers["counts"].copy()
    cell2loc_h5ad = out_dir / f"{prefix}.cell2location_spatial_counts.h5ad"

    # Optional ROI selection by a provided cell/spot ID list
    roi_cfg = cfg.get("roi", {}) or {}
    roi_ids_txt = roi_cfg.get("cell_ids_txt", "")
    if roi_ids_txt:
        roi_path = Path(roi_ids_txt)
        if not roi_path.exists():
            raise SystemExit(f"ERROR: roi.cell_ids_txt not found: {roi_path}")
        roi_ids = set(_read_id_list(roi_path))
        keep = [cid for cid in adata_counts.obs_names if cid in roi_ids]
        _log(f"[roi] keeping {len(keep)} ids from {roi_path}", log_fh)
        adata_counts = adata_counts[keep, :].copy()
        cell2loc_h5ad = out_dir / f"{prefix}.cell2location_spatial_counts_roi.h5ad"

    sc.write_h5ad(str(cell2loc_h5ad), adata_counts)
    _log(f"[save] cell2location spatial counts={cell2loc_h5ad}", log_fh)

    # provenance
    cfg_out = dict(cfg)
    cfg_out["resolved"] = {
        "config": str(cfg_path),
        "input_h5ad": str(input_h5ad),
        "out_dir": str(out_dir),
        "prefix": prefix,
        "orig_qc_h5ad": str(orig_qc_h5ad),
        "processed_h5ad": str(processed_h5ad),
        "cell2location_h5ad": str(cell2loc_h5ad),
    }
    with open(out_dir / "params_used_fig4_cell2location_spatial_prep.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg_out, f, sort_keys=False, allow_unicode=True)

    _write_text(out_dir / "env_versions.txt", "\n".join([
        f"python={sys.version.replace(os.linesep, ' ')}",
        f"scanpy={sc.__version__}",
        f"anndata={ad.__version__}",
        f"numpy={np.__version__}",
        f"pandas={pd.__version__}",
    ]) + "\n")

    _log("[done] finished", log_fh)
    log_fh.close()


if __name__ == "__main__":
    main()
