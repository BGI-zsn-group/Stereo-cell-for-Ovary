#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig4 spatial preprocessing for cell2location.
Aligned to the user's notebook logic (ROI ignored, image loading ignored), Python 3.8 compatible.

Notebook-aligned logic:
1) read one adjusted .h5ad
2) aggregate duplicated genes by `gene_name_col`
3) QC metrics + filter_cells/filter_genes
4) snapshot `adata_orig` after QC filtering but before normalization
5) round1 on all cells:
   normalize_total -> log1p -> raw copy -> HVG -> score/regress cell cycle -> scale -> PCA -> spatial neighbors -> neighbors -> UMAP -> leiden_10/leiden_15
6) subset by round1 leiden_15 exclusion (and optional ROI cell list)
7) round2 on subset:
   adata_subset = adata_subset.raw.to_adata()   # back to round1 raw matrix
   adata_subset.raw = adata_subset.copy()
   HVG -> scale -> PCA -> spatial neighbors -> neighbors -> UMAP -> leiden_10/leiden_15
   (no second normalize/log1p, no second cell-cycle scoring/regression)
8) final c2l input is built from adata_orig[adata_subset.obs_names, :] with full genes + copied obs/var/obsm/uns
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

try:
    import squidpy as sq
except Exception:
    sq = None


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
        raise SystemExit("ERROR: gene_name_col '{}' not found in adata.var".format(gene_name_col))

    gene_names = pd.Series(adata.var[gene_name_col]).astype(str).fillna("")
    valid = gene_names != ""
    drop_gene_names = [str(x) for x in (drop_gene_names or []) if str(x) != ""]
    if drop_gene_names:
        drop_set = set(drop_gene_names)
        valid &= ~gene_names.isin(drop_set)
        _log("[gene] dropping {} genes by name filter: {}".format(int((~valid).sum()), sorted(drop_set)), log_fh)

    adata = adata[:, valid.to_numpy()].copy()
    gene_names = pd.Series(adata.var[gene_name_col]).astype(str).fillna("")

    X = adata.X.tocsc() if sparse.issparse(adata.X) else sparse.csc_matrix(adata.X)
    unique_genes, inverse_idx = np.unique(gene_names.values, return_inverse=True)
    _log("[gene] aggregating duplicated genes: {} -> {}".format(adata.n_vars, len(unique_genes)), log_fh)

    X_agg = sparse.lil_matrix((X.shape[0], len(unique_genes)), dtype=X.dtype)
    for i in range(len(unique_genes)):
        cols = np.where(inverse_idx == i)[0]
        if len(cols) == 1:
            X_agg[:, i] = X[:, cols[0]]
        else:
            X_agg[:, i] = X[:, cols].sum(axis=1)
    X_agg = X_agg.tocsr()

    new_var = pd.DataFrame(index=pd.Index(unique_genes.astype(str), name=None))
    out = ad.AnnData(
        X=X_agg,
        obs=adata.obs.copy(),
        var=new_var,
        obsm={k: (v.copy() if hasattr(v, "copy") else v) for k, v in adata.obsm.items()},
        uns=adata.uns.copy(),
    )
    out.var_names = out.var.index.astype(str)
    return out


def _ensure_counts_layer(adata: ad.AnnData):
    X = adata.X.copy() if hasattr(adata.X, "copy") else np.array(adata.X, copy=True)
    adata.layers["counts"] = X


def _write_h5ad(adata: ad.AnnData, path: Path):
    _ensure_dir(path.parent)
    adata.write_h5ad(str(path))


def _run_neighbors_umap_leiden(adata: ad.AnnData, cfg: dict, log_fh=None, prefix: str = ""):
    sc.pp.scale(adata, max_value=float(cfg.get("scale_max_value", 10)))
    sc.tl.pca(adata, svd_solver=str(cfg.get("pca_solver", "arpack")))

    if sq is not None and "spatial" in adata.obsm.keys() and bool(cfg.get("do_spatial_neighbors", True)):
        try:
            sq.gr.spatial_neighbors(
                adata,
                coord_type=str(cfg.get("spatial_coord_type", "generic")),
                delaunay=bool(cfg.get("spatial_delaunay", True)),
            )
            _log("{}[graph] built squidpy spatial neighbors".format(prefix), log_fh)
        except Exception as e:
            _log("{}[warn] squidpy spatial_neighbors failed: {}".format(prefix, e), log_fh)

    sc.pp.neighbors(adata, n_neighbors=int(cfg.get("n_neighbors", 20)), n_pcs=int(cfg.get("n_pcs", 30)))
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=float(cfg.get("leiden_resolution_10", 1.0)), key_added=str(cfg.get("leiden_key_10", "leiden_10")))
    sc.tl.leiden(adata, resolution=float(cfg.get("leiden_resolution_15", 1.5)), key_added=str(cfg.get("leiden_key_15", "leiden_15")))


def _score_and_regress_cell_cycle(adata: ad.AnnData, cfg: dict, log_fh=None, prefix: str = ""):
    s_genes = [g for g in _as_list(cfg.get("s_genes")) if g in adata.var_names]
    g2m_genes = [g for g in _as_list(cfg.get("g2m_genes")) if g in adata.var_names]
    if len(s_genes) == 0 or len(g2m_genes) == 0:
        _log("{}[cell-cycle] skip: no matched S/G2M genes".format(prefix), log_fh)
        return
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, copy=False)
    sc.pp.regress_out(adata, ["S_score", "G2M_score"])
    _log("{}[cell-cycle] scored/regressed with {} S genes and {} G2M genes".format(prefix, len(s_genes), len(g2m_genes)), log_fh)


def _round1_process(adata: ad.AnnData, cfg: dict, log_fh=None) -> ad.AnnData:
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
        _score_and_regress_cell_cycle(adata, cfg, log_fh, prefix="[round1]")
    _run_neighbors_umap_leiden(adata, cfg, log_fh, prefix="[round1]")
    return adata


def _round2_process_notebook(adata_subset: ad.AnnData, cfg: dict, log_fh=None) -> ad.AnnData:
    if adata_subset.raw is None:
        raise SystemExit("ERROR: adata_subset.raw is None; notebook-style round2 requires raw.to_adata()")
    adata_subset = adata_subset.raw.to_adata()
    _log("[subset-round2] returned subset to raw counts before subset reprocessing", log_fh)
    adata_subset.raw = adata_subset.copy()
    sc.pp.highly_variable_genes(
        adata_subset,
        flavor=str(cfg.get("subset_hvg_flavor", cfg.get("hvg_flavor", "seurat_v3"))),
        n_top_genes=int(cfg.get("subset_hvg_n_top_genes", cfg.get("hvg_n_top_genes", 2000))),
    )
    if "highly_variable" in adata_subset.var.columns:
        adata_subset = adata_subset[:, adata_subset.var["highly_variable"].to_numpy()].copy()
    # Notebook round2 does NOT normalize/log1p again and does NOT do cell-cycle regress again.
    _run_neighbors_umap_leiden(adata_subset, cfg, log_fh, prefix="[round2]")
    return adata_subset


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = _read_yaml(Path(args.config))
    in_h5ad = Path(cfg.get("input_h5ad", "data/Fig4/spatial/B04372C211.adjusted.cellbin.h5ad"))
    if not in_h5ad.exists():
        raise SystemExit("ERROR: missing input_h5ad: {}".format(in_h5ad))

    out_dir = Path(cfg.get("out_dir", "results/Fig4/cell2location_spatial"))
    _ensure_dir(out_dir)
    prefix = str(cfg.get("prefix", in_h5ad.name.replace(".adjusted.cellbin.h5ad", "").replace(".h5ad", "")))
    log_fh = open(out_dir / (prefix + ".fig4_cell2location_spatial.log"), "w", encoding="utf-8")

    _log("[start] config={}".format(args.config), log_fh)
    _log("[io] input_h5ad={}".format(in_h5ad), log_fh)
    _log("[io] out_dir={}".format(out_dir), log_fh)
    _log("[io] prefix={}".format(prefix), log_fh)

    adata = sc.read_h5ad(str(in_h5ad))
    _log("[load] adata: n_obs={}, n_vars={}".format(adata.n_obs, adata.n_vars), log_fh)

    adata = _aggregate_by_gene_name(
        adata,
        gene_name_col=str(cfg.get("gene_name_col", "real_gene_name")),
        drop_gene_names=_as_list(cfg.get("drop_gene_names", [])),
        log_fh=log_fh,
    )
    _log("[gene] after aggregation: n_vars={}".format(adata.n_vars), log_fh)

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.contains("^HB[AB]", regex=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    min_genes = int(cfg.get("min_genes", 100))
    min_cells = int(cfg.get("min_cells", 3))
    _log("[qc] filter_cells min_genes={}; filter_genes min_cells={}".format(min_genes, min_cells), log_fh)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    adata_orig = adata.copy()
    _ensure_counts_layer(adata_orig)
    orig_qc_h5ad_cfg = cfg.get("out_orig_qc_h5ad", None)
    orig_qc_h5ad = Path(orig_qc_h5ad_cfg) if orig_qc_h5ad_cfg not in (None, "") else (out_dir / (prefix + ".orig_qc_counts.h5ad"))
    _write_h5ad(adata_orig, orig_qc_h5ad)
    _log("[save] orig_qc_h5ad={} n_obs={} n_vars={}".format(orig_qc_h5ad, adata_orig.n_obs, adata_orig.n_vars), log_fh)

    adata = _round1_process(adata, cfg, log_fh)
    processed_h5ad_cfg = cfg.get("out_processed_h5ad", None)
    processed_h5ad = Path(processed_h5ad_cfg) if processed_h5ad_cfg not in (None, "") else (out_dir / (prefix + ".adjusted_processed.h5ad"))
    _write_h5ad(adata, processed_h5ad)
    _log("[save] processed_h5ad={} n_obs={} n_vars={}".format(processed_h5ad, adata.n_obs, adata.n_vars), log_fh)

    subset_cells = list(adata.obs_names)
    roi_cfg = cfg.get("roi", {}) or {}
    cell_ids_txt = roi_cfg.get("cell_ids_txt", None)
    if cell_ids_txt not in (None, ""):
        roi_ids = set(_read_id_list(Path(cell_ids_txt)))
        subset_cells = [c for c in subset_cells if c in roi_ids]
        _log("[roi] keep by cell_ids_txt -> {} cells".format(len(subset_cells)), log_fh)

    subset_cfg = cfg.get("subset", {}) or {}
    leiden_key = str(subset_cfg.get("leiden_key", cfg.get("leiden_key_15", "leiden_15")))
    exclude_clusters = set(_as_list(subset_cfg.get("exclude_clusters", [])))
    if leiden_key in adata.obs.columns and exclude_clusters:
        subset_cells = [c for c in subset_cells if str(adata.obs.loc[c, leiden_key]) not in exclude_clusters]
        _log("[subset-round1] leiden_key={} exclude={} -> {} cells".format(leiden_key, sorted(exclude_clusters), len(subset_cells)), log_fh)
    adata_subset = adata[subset_cells, :].copy()

    subset_raw_h5ad_cfg = cfg.get("out_subset_raw_h5ad", None)
    subset_raw_h5ad = Path(subset_raw_h5ad_cfg) if subset_raw_h5ad_cfg not in (None, "") else (out_dir / (prefix + ".adjusted_subset_raw.h5ad"))

    if bool(cfg.get("subset_reprocess", True)):
        adata_subset = _round2_process_notebook(adata_subset, cfg, log_fh)
        # Save raw-count snapshot after raw.to_adata and before HVG? notebook saves the raw subset before HVG.
        # To preserve notebook behavior, reconstruct and save that snapshot separately.
        # We can't recover the pre-HVG round2 snapshot from adata_subset now, so create from round1 subset raw.
        tmp_subset_raw = adata[subset_cells, :].copy().raw.to_adata()
        tmp_subset_raw.raw = tmp_subset_raw.copy()
        _write_h5ad(tmp_subset_raw, subset_raw_h5ad)
        _log("[save] subset_raw_h5ad={} n_obs={} n_vars={}".format(subset_raw_h5ad, tmp_subset_raw.n_obs, tmp_subset_raw.n_vars), log_fh)
    else:
        _write_h5ad(adata_subset, subset_raw_h5ad)
        _log("[save] subset_raw_h5ad={} n_obs={} n_vars={}".format(subset_raw_h5ad, adata_subset.n_obs, adata_subset.n_vars), log_fh)

    subset_processed_h5ad_cfg = cfg.get("out_subset_processed_h5ad", None)
    subset_processed_h5ad = Path(subset_processed_h5ad_cfg) if subset_processed_h5ad_cfg not in (None, "") else (out_dir / (prefix + ".adjusted_processed_subset.h5ad"))
    _write_h5ad(adata_subset, subset_processed_h5ad)
    _log("[save] subset_processed_h5ad={} n_obs={} n_vars={}".format(subset_processed_h5ad, adata_subset.n_obs, adata_subset.n_vars), log_fh)

    subset_cells_final = list(adata_subset.obs_names)
    obs_subset = adata_orig.obs.loc[subset_cells_final].copy()
    var_full = adata_orig.var.copy()
    X_full = adata_orig[subset_cells_final, :].X.copy()

    obsm_subset = {}
    for key in adata_orig.obsm.keys():
        arr = adata_orig.obsm[key]
        idx = adata_orig.obs_names.get_indexer_for(subset_cells_final)
        if hasattr(arr, "__getitem__"):
            piece = arr[idx]
            obsm_subset[key] = piece.copy() if hasattr(piece, "copy") else piece

    adata_counts = ad.AnnData(
        X=X_full,
        obs=obs_subset,
        var=var_full,
        obsm=obsm_subset,
        uns=adata_orig.uns.copy(),
    )
    _ensure_counts_layer(adata_counts)
    if "orig.ident" not in adata_counts.obs.columns:
        adata_counts.obs["orig.ident"] = prefix

    cell2loc_h5ad_cfg = cfg.get("out_cell2location_h5ad", None)
    cell2loc_h5ad = Path(cell2loc_h5ad_cfg) if cell2loc_h5ad_cfg not in (None, "") else (out_dir / (prefix + ".cell2location_spatial_counts.h5ad"))
    _write_h5ad(adata_counts, cell2loc_h5ad)
    _log("[save] cell2location_h5ad={} n_obs={} n_vars={}".format(cell2loc_h5ad, adata_counts.n_obs, adata_counts.n_vars), log_fh)

    with open(out_dir / "params_used_fig4_cell2location_spatial.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)
    _log("[done] spatial preprocessing finished", log_fh)
    log_fh.close()


if __name__ == "__main__":
    main()
