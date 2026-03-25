#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fig4: Cell2location spatial transcriptomics preprocessing (notebook-aligned version)

Key behavior
- Use `input_h5ad` for QC / normalization / HVG / clustering / subset selection.
- Use `input_h5ad_orig` (if provided or auto-detected) to rebuild the final
  cell2location input with FULL genes + raw counts on the selected subset,
  matching the original notebook logic more closely.
- Compatible with Python 3.8.
"""

import argparse
import datetime as _dt
import os
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
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
    return open(log_path, "w", encoding="utf-8")


def _log(msg: str, fh) -> None:
    line = f"[{_now()}] {msg}"
    print(line)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


def _aggregate_by_gene_name(
    adata: ad.AnnData,
    gene_name_col: str,
    drop_gene_names: Optional[Sequence[str]],
    log_fh,
) -> ad.AnnData:
    if gene_name_col not in adata.var.columns:
        _log(f"[gene] var does not contain '{gene_name_col}', skip aggregation.", log_fh)
        adata.var_names_make_unique()
        return adata

    gene_names = adata.var[gene_name_col].astype(str).values
    if drop_gene_names:
        drop_set = set(str(x) for x in drop_gene_names)
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
    for layer_name, layer in getattr(adata, "layers", {}).items():
        try:
            L = layer
            if not sparse.issparse(L):
                L = sparse.csr_matrix(L)
            L = L.tocoo()
            L_new = sparse.coo_matrix((L.data, (L.row, inverse_idx[L.col])), shape=(adata.n_obs, len(unique_genes))).tocsr()
            adata_new.layers[layer_name] = L_new
        except Exception:
            pass

    adata_new.var_names = adata_new.var.index.astype(str)
    adata_new.var_names_make_unique()
    return adata_new


def _add_qc_flags(adata: ad.AnnData) -> None:
    vn = adata.var_names.astype(str).str.upper()
    adata.var["mt"] = vn.str.startswith("MT-")
    adata.var["ribo"] = vn.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = vn.str.contains(r"^HB[AB]", regex=True)


def _read_id_list(path: Path) -> List[str]:
    ids = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            ids.append(s)
    return ids


def _guess_orig_h5ad(input_h5ad: Path) -> Optional[Path]:
    s = input_h5ad.name
    if s.endswith(".h5ad"):
        candidate = input_h5ad.with_name(s[:-5] + "_orig.h5ad")
        if candidate.exists():
            return candidate
    return None


def _subset_obsm(orig: ad.AnnData, subset_cells: Sequence[str]) -> dict:
    idx = orig.obs_names.get_indexer(subset_cells)
    if np.any(idx < 0):
        missing = [subset_cells[i] for i, v in enumerate(idx) if v < 0][:10]
        raise SystemExit(f"ERROR: {int((idx < 0).sum())} subset cells not found in orig object. Example: {missing}")
    out = {}
    for key in orig.obsm.keys():
        arr = orig.obsm[key]
        if hasattr(arr, "iloc"):
            out[key] = arr.iloc[idx].copy()
        else:
            out[key] = arr[idx].copy()
    return out


def _build_final_counts_from_orig(orig: ad.AnnData, subset_cells: Sequence[str], log_fh) -> ad.AnnData:
    _log(f"[final] rebuilding final counts from orig object with full genes: n_obs={orig.n_obs}, n_vars={orig.n_vars}", log_fh)
    obs = orig.obs.loc[list(subset_cells)].copy()
    var = orig.var.copy()
    obsm = _subset_obsm(orig, list(subset_cells))
    uns = orig.uns.copy()

    if "counts" in orig.layers:
        X = orig.layers["counts"][orig.obs_names.get_indexer(subset_cells), :].copy()
    else:
        X = orig[list(subset_cells), :].X.copy()

    adata_counts = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns)
    adata_counts.layers["counts"] = X.copy() if hasattr(X, "copy") else X
    return adata_counts


def _maybe_apply_roi(adata_counts: ad.AnnData, roi_ids_txt: str, log_fh) -> ad.AnnData:
    if not roi_ids_txt:
        return adata_counts
    roi_path = Path(roi_ids_txt)
    if not roi_path.exists():
        raise SystemExit(f"ERROR: roi.cell_ids_txt not found: {roi_path}")
    roi_ids = set(_read_id_list(roi_path))
    keep = [cid for cid in adata_counts.obs_names if cid in roi_ids]
    _log(f"[roi] keeping {len(keep)} ids from {roi_path}", log_fh)
    return adata_counts[keep, :].copy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="YAML config path")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    cfg = _read_yaml(cfg_path)

    input_h5ad = Path(cfg.get("input_h5ad", "B04372C211.adjusted.cellbin.h5ad"))
    input_h5ad_orig_cfg = str(cfg.get("input_h5ad_orig", "")).strip()
    input_h5ad_orig = Path(input_h5ad_orig_cfg) if input_h5ad_orig_cfg else None
    if input_h5ad_orig is None:
        input_h5ad_orig = _guess_orig_h5ad(input_h5ad)

    out_dir = Path(cfg.get("out_dir", "results/cell2location_spatial"))
    prefix = str(cfg.get("prefix", input_h5ad.stem))

    _ensure_dir(out_dir)
    log_fh = _log_setup(out_dir / f"{prefix}.spatial_prep.log")

    _log(f"[start] config={cfg_path}", log_fh)
    _log(f"[io] input_h5ad={input_h5ad}", log_fh)
    _log(f"[io] input_h5ad_orig={input_h5ad_orig if input_h5ad_orig else 'None'}", log_fh)
    _log(f"[io] out_dir={out_dir}", log_fh)
    _log(f"[io] prefix={prefix}", log_fh)

    if not input_h5ad.exists():
        raise SystemExit(f"ERROR: missing input_h5ad: {input_h5ad}")
    if input_h5ad_orig is not None and not input_h5ad_orig.exists():
        raise SystemExit(f"ERROR: missing input_h5ad_orig: {input_h5ad_orig}")

    spatial_key = cfg.get("spatial_obsm_key", "spatial")

    # Load processing object
    adata = sc.read_h5ad(str(input_h5ad))
    _log(f"[load] adata: n_obs={adata.n_obs}, n_vars={adata.n_vars}", log_fh)
    if spatial_key not in adata.obsm_keys():
        raise SystemExit(f"ERROR: missing adata.obsm['{spatial_key}'] (spatial coordinates).")

    # Aggregate by real_gene_name like notebook
    gene_name_col = cfg.get("gene_name_col", "real_gene_name")
    drop_gene_names = cfg.get("drop_gene_names", ["a"])
    if drop_gene_names is None:
        drop_gene_names = []
    adata = _aggregate_by_gene_name(adata, gene_name_col=gene_name_col, drop_gene_names=drop_gene_names, log_fh=log_fh)
    _log(f"[gene] after aggregation: n_vars={adata.n_vars}", log_fh)

    # Store raw counts before normalization, like notebook's adata_orig.layers['counts']
    adata.layers["counts"] = adata.X.copy() if sparse.issparse(adata.X) else np.array(adata.X, copy=True)

    _add_qc_flags(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=False)

    min_genes = int(cfg.get("min_genes", 100))
    min_cells = int(cfg.get("min_cells", 3))
    _log(f"[qc] filter_cells min_genes={min_genes}; filter_genes min_cells={min_cells}", log_fh)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    orig_qc_h5ad = out_dir / f"{prefix}.orig_qc_counts.h5ad"
    adata.write_h5ad(str(orig_qc_h5ad))
    _log(f"[save] orig_qc_counts={orig_qc_h5ad}", log_fh)

    processed_h5ad = out_dir / f"{prefix}.processed.h5ad"
    subset_processed_h5ad = out_dir / f"{prefix}.subset_processed.h5ad"

    do_processing = bool(cfg.get("do_processing", True))
    if do_processing:
        adata_proc = adata.copy()

        target_sum = float(cfg.get("normalize_target_sum", 1e4))
        sc.pp.normalize_total(adata_proc, target_sum=target_sum)
        sc.pp.log1p(adata_proc)
        adata_proc.raw = adata_proc.copy()

        n_top_genes = int(cfg.get("hvg_n_top_genes", 2000))
        sc.pp.highly_variable_genes(adata_proc, flavor="seurat_v3", n_top_genes=n_top_genes)
        adata_proc = adata_proc[:, adata_proc.var["highly_variable"]].copy()

        do_cell_cycle = bool(cfg.get("do_cell_cycle", True))
        if do_cell_cycle:
            s_genes = list(cfg.get("s_genes", []))
            g2m_genes = list(cfg.get("g2m_genes", []))
            s_use = [g for g in s_genes if g in adata_proc.var_names]
            g2m_use = [g for g in g2m_genes if g in adata_proc.var_names]
            if s_use and g2m_use:
                sc.tl.score_genes_cell_cycle(adata_proc, s_genes=s_use, g2m_genes=g2m_use)
                _log(f"[cc] scored cell cycle genes: S={len(s_use)}, G2M={len(g2m_use)}", log_fh)
                sc.pp.regress_out(adata_proc, ["S_score", "G2M_score"])
                _log("[cc] regressed out S_score and G2M_score", log_fh)
            else:
                _log("[cc] skipped cell cycle scoring: no overlapping genes found", log_fh)

        scale_max_value = float(cfg.get("scale_max_value", 10))
        sc.pp.scale(adata_proc, max_value=scale_max_value)

        n_pcs = int(cfg.get("n_pcs", 30))
        n_neighbors = int(cfg.get("n_neighbors", 20))
        sc.tl.pca(adata_proc, svd_solver="arpack")
        sc.pp.neighbors(adata_proc, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(adata_proc)

        leiden_resolutions = cfg.get("leiden_resolutions", {"leiden_10": 1.0, "leiden_15": 1.5}) or {}
        for key_added, resolution in leiden_resolutions.items():
            sc.tl.leiden(adata_proc, resolution=float(resolution), key_added=str(key_added))
            _log(f"[cluster] leiden key={key_added}, resolution={resolution}", log_fh)

        adata_proc.write_h5ad(str(processed_h5ad))
        _log(f"[save] processed={processed_h5ad}", log_fh)

        subset_cfg = cfg.get("subset", {}) or {}
        leiden_key = str(subset_cfg.get("leiden_key", "leiden_15"))
        exclude_clusters = [str(x) for x in subset_cfg.get("exclude_clusters", [])]
        if exclude_clusters:
            if leiden_key not in adata_proc.obs.columns:
                raise SystemExit(f"ERROR: subset.leiden_key='{leiden_key}' not found in processed adata.obs.")
            keep_mask = ~adata_proc.obs[leiden_key].astype(str).isin(exclude_clusters)
            adata_proc_sub = adata_proc[keep_mask].copy()
            subset_cells = adata_proc_sub.obs_names.tolist()
            adata_proc_sub.write_h5ad(str(subset_processed_h5ad))
            _log(f"[subset] by clusters on {leiden_key}: exclude={exclude_clusters}; kept {len(subset_cells)} cells/spots; saved {subset_processed_h5ad}", log_fh)
        else:
            subset_cells = adata_proc.obs_names.tolist()
            _log(f"[subset] no exclude_clusters provided; keep all {len(subset_cells)} cells/spots", log_fh)
    else:
        subset_cells = adata.obs_names.tolist()
        _log(f"[proc] do_processing=False; skip processing and keep all {len(subset_cells)}", log_fh)

    # Build final c2l object
    use_orig_full_counts = bool(cfg.get("use_orig_full_counts", True))
    roi_cfg = cfg.get("roi", {}) or {}
    roi_ids_txt = str(roi_cfg.get("cell_ids_txt", ""))

    if use_orig_full_counts and input_h5ad_orig is not None:
        adata_orig = sc.read_h5ad(str(input_h5ad_orig))
        _log(f"[load] adata_orig: n_obs={adata_orig.n_obs}, n_vars={adata_orig.n_vars}", log_fh)
        if spatial_key not in adata_orig.obsm_keys():
            raise SystemExit(f"ERROR: missing adata_orig.obsm['{spatial_key}'] (spatial coordinates).")
        adata_counts = _build_final_counts_from_orig(adata_orig, subset_cells, log_fh)
    else:
        if use_orig_full_counts and input_h5ad_orig is None:
            _log("[final] input_h5ad_orig not provided/found; fall back to current object for final counts", log_fh)
        adata_counts = adata[subset_cells, :].copy()
        if "counts" in adata_counts.layers:
            adata_counts.X = adata_counts.layers["counts"].copy()
        else:
            adata_counts.layers["counts"] = adata_counts.X.copy() if sparse.issparse(adata_counts.X) else np.array(adata_counts.X, copy=True)

    adata_counts = _maybe_apply_roi(adata_counts, roi_ids_txt=roi_ids_txt, log_fh=log_fh)

    cell2loc_h5ad = out_dir / f"{prefix}.cell2location_spatial_counts.h5ad"
    if roi_ids_txt:
        cell2loc_h5ad = out_dir / f"{prefix}.cell2location_spatial_counts_roi.h5ad"

    adata_counts.write_h5ad(str(cell2loc_h5ad))
    _log(f"[save] cell2location spatial counts={cell2loc_h5ad}", log_fh)
    _log(f"[final] cell2location object: n_obs={adata_counts.n_obs}, n_vars={adata_counts.n_vars}", log_fh)

    cfg_out = dict(cfg)
    cfg_out["resolved"] = {
        "config": str(cfg_path),
        "input_h5ad": str(input_h5ad),
        "input_h5ad_orig": str(input_h5ad_orig) if input_h5ad_orig else "",
        "out_dir": str(out_dir),
        "prefix": prefix,
        "orig_qc_h5ad": str(orig_qc_h5ad),
        "processed_h5ad": str(processed_h5ad),
        "subset_processed_h5ad": str(subset_processed_h5ad),
        "cell2location_h5ad": str(cell2loc_h5ad),
    }
    with open(out_dir / "params_used_fig4_cell2location_spatial_prep.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg_out, f, sort_keys=False, allow_unicode=True)

    _write_text(
        out_dir / "env_versions.txt",
        "\n".join([
            f"python={sys.version.replace(os.linesep, ' ')}",
            f"scanpy={sc.__version__}",
            f"anndata={ad.__version__}",
            f"numpy={np.__version__}",
            f"pandas={pd.__version__}",
        ]) + "\n",
    )

    _log("[done] finished", log_fh)
    log_fh.close()


if __name__ == "__main__":
    main()
