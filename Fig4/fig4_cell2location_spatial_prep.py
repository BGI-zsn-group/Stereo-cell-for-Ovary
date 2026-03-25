#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig4 spatial preprocessing for cell2location.
Notebook-aligned, Python 3.8 compatible.

Key alignment points:
1) QC snapshot adata_orig is taken AFTER filter_cells/filter_genes and BEFORE normalize/log1p.
2) Round1 uses notebook S/G2M gene lists by default and performs score_genes_cell_cycle + regress_out.
3) Round2 follows the notebook exactly: raw.to_adata() -> HVG -> scale -> PCA -> neighbors -> UMAP -> leiden
   (NO normalize_total/log1p, NO cell-cycle regression).
4) Final full-count object is rebuilt from adata_orig[subset_cells, :] and copies obs/var/obsm/uns.
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


NOTEBOOK_S_GENES = [
    "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6",
    "Cdca7", "Dtl", "Prim1", "Uhrf1", "Cenpu", "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1",
    "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2",
    "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1",
    "Chaf1b", "Brip1", "E2f8"
]
NOTEBOOK_G2M_GENES = [
    "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2",
    "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Smc4", "Ccnb2",
    "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b",
    "Hjurp", "Cdca3", "Jpt1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2",
    "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln",
    "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Cbx5", "Cenpa"
]


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

    adata_out = ad.AnnData(
        X=X_agg,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Index(unique_genes.astype(str), name=None)),
        obsm={k: (v.copy() if hasattr(v, "copy") else v) for k, v in adata.obsm.items()},
        uns=adata.uns.copy(),
    )
    adata_out.var_names = adata_out.var.index.astype(str)
    return adata_out


def _ensure_counts_layer(adata: ad.AnnData):
    X = adata.X.copy() if hasattr(adata.X, "copy") else np.array(adata.X, copy=True)
    adata.layers["counts"] = X


def _run_graph(adata: ad.AnnData, cfg: dict, log_fh=None, prefix: str = ""):
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


def _round1_process(adata: ad.AnnData, cfg: dict, log_fh=None) -> ad.AnnData:
    sc.pp.normalize_total(adata, target_sum=float(cfg.get("normalize_target_sum", 1e4)))
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata, flavor=str(cfg.get("hvg_flavor", "seurat_v3")), n_top_genes=int(cfg.get("hvg_n_top_genes", 2000)))
    if "highly_variable" in adata.var.columns:
        adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()

    s_genes = [g for g in _as_list(cfg.get("s_genes", NOTEBOOK_S_GENES)) if g in adata.var_names]
    g2m_genes = [g for g in _as_list(cfg.get("g2m_genes", NOTEBOOK_G2M_GENES)) if g in adata.var_names]
    if len(s_genes) == 0 or len(g2m_genes) == 0:
        raise SystemExit("ERROR: round1 cell-cycle genes matched 0 genes after HVG selection; this is not notebook-equivalent.")
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, copy=False)
    sc.pp.regress_out(adata, ["S_score", "G2M_score"])
    _log("[round1][cell-cycle] scored/regressed with {} S genes and {} G2M genes".format(len(s_genes), len(g2m_genes)), log_fh)
    _run_graph(adata, cfg, log_fh, prefix="[round1]")
    return adata


def _round2_process_notebook(adata_subset: ad.AnnData, cfg: dict, log_fh=None) -> ad.AnnData:
    # Notebook round2 does NOT normalize/log1p and does NOT run cell-cycle regression again.
    sc.pp.highly_variable_genes(adata_subset, flavor=str(cfg.get("hvg_flavor", "seurat_v3")), n_top_genes=int(cfg.get("hvg_n_top_genes", 2000)))
    if "highly_variable" in adata_subset.var.columns:
        adata_subset = adata_subset[:, adata_subset.var["highly_variable"].to_numpy()].copy()
    _run_graph(adata_subset, cfg, log_fh, prefix="[round2]")
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
    sc.pp.filter_cells(adata, min_genes=int(cfg.get("min_genes", 100)))
    sc.pp.filter_genes(adata, min_cells=int(cfg.get("min_cells", 3)))
    _log("[qc] filter_cells min_genes={}; filter_genes min_cells={}".format(int(cfg.get("min_genes", 100)), int(cfg.get("min_cells", 3))), log_fh)

    adata_orig = adata.copy()
    _ensure_counts_layer(adata_orig)
    orig_qc_h5ad = out_dir / (prefix + ".orig_qc_counts.h5ad")
    adata_orig.write_h5ad(str(orig_qc_h5ad))
    _log("[save] orig_qc_h5ad={} n_obs={} n_vars={}".format(orig_qc_h5ad, adata_orig.n_obs, adata_orig.n_vars), log_fh)

    adata = _round1_process(adata, cfg, log_fh)
    processed_h5ad = out_dir / (prefix + ".adjusted_processed.h5ad")
    adata.write_h5ad(str(processed_h5ad))
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

    if adata_subset.raw is None:
        raise SystemExit("ERROR: adata_subset.raw is None; cannot run notebook-style subset raw->to_adata()")
    adata_subset = adata_subset.raw.to_adata()
    _log("[subset-round2] returned subset to raw counts before subset reprocessing", log_fh)
    adata_subset.raw = adata_subset.copy()
    subset_raw_h5ad = out_dir / (prefix + ".adjusted_subset_raw.h5ad")
    adata_subset.write_h5ad(str(subset_raw_h5ad))
    _log("[save] subset_raw_h5ad={} n_obs={} n_vars={}".format(subset_raw_h5ad, adata_subset.n_obs, adata_subset.n_vars), log_fh)

    adata_subset = _round2_process_notebook(adata_subset, cfg, log_fh)
    subset_processed_h5ad = out_dir / (prefix + ".adjusted_processed_subset.h5ad")
    adata_subset.write_h5ad(str(subset_processed_h5ad))
    _log("[save] subset_processed_h5ad={} n_obs={} n_vars={}".format(subset_processed_h5ad, adata_subset.n_obs, adata_subset.n_vars), log_fh)

    subset_cells_final = list(adata_subset.obs_names)
    obs_subset = adata_orig.obs.loc[subset_cells_final].copy()
    var_full = adata_orig.var.copy()
    X_full = adata_orig[subset_cells_final, :].layers["counts"].copy()
    obsm_subset = {key: adata_orig.obsm[key][adata_orig.obs_names.get_indexer_for(subset_cells_final)] for key in adata_orig.obsm.keys()}
    adata_counts = ad.AnnData(X=X_full, obs=obs_subset, var=var_full, obsm=obsm_subset, uns=adata_orig.uns.copy())
    _ensure_counts_layer(adata_counts)
    if "orig.ident" not in adata_counts.obs.columns:
        adata_counts.obs["orig.ident"] = prefix
    cell2loc_h5ad = out_dir / (prefix + ".cell2location_spatial_counts.h5ad")
    adata_counts.write_h5ad(str(cell2loc_h5ad))
    _log("[save] cell2location_h5ad={} n_obs={} n_vars={}".format(cell2loc_h5ad, adata_counts.n_obs, adata_counts.n_vars), log_fh)
    _log("[done] spatial preprocessing finished", log_fh)
    log_fh.close()


if __name__ == "__main__":
    main()
