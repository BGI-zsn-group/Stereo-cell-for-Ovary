# Fig3 — Monocle3 / Length / GO / SCENIC / scTour

本目录用于复现论文 **Fig.3** 的轨迹与调控分析流程。推荐入口是 `run_fig3_combined.sh`。

## 目录结构

```text
Fig3/
  configs/
    fig3_combined.yaml
  fig3_monocle3_modules.R
  fig3_length_pseudotime.R
  fig3_go_enrichment.R
  fig3_scenic_rds_to_csv.R
  fig3_scenic_pyscenic.py
  fig3_scenic_downstream.R
  fig3_rds_to_h5ad.R
  fig3_sctour_tnode.py
  run_fig3_combined.sh
```

## 主流程

- `monocle3`：Monocle3 轨迹、pseudotime、graph_test、gene modules
- `length`：Length_px vs pseudotime 曲线与相关性
- `go`：gene modules 的 GO enrichment
- `scenic1,scenic2,scenic_downstream`：SCENIC 表达矩阵导出、pySCENIC、RSS downstream
- `rds2h5ad,sctour`：Seurat 转 h5ad，并运行 scTour

## Monocle3 当前逻辑

`fig3_monocle3_modules.R` 现已同步为你最新的分析逻辑：

1. `obj_oo.rds -> as.cell_data_set()`
2. `preprocess_cds(num_dim=50)` 与 `reduce_dimension()`
3. 强制单 partition，并复用 Seurat 的 `seurat_clusters` 与 UMAP 坐标
4. 以 `root_cluster` 作为起点执行 `learn_graph()` 与 `order_cells()`
5. `graph_test()` 得到随主图变化的 DEG
6. 先做一次 **预模块划分**（默认 `prefilter_module_resolution=0.001`）
7. 基于 `stage` 的 `FindAllMarkers()`：
   - early stages：默认 `EGO, GO1`
   - late stages：默认 `GO3, FGO`
8. 对目标模块（默认 module `3`）执行阶段感知过滤：
   - 找出该模块中既不是 early/late marker，或同时属于 early+late marker 的基因
   - 将这些基因从 `pr_deg_ids` 中移除
9. 用过滤后的 `pr_deg_ids` 再跑一次 **最终模块划分**（默认 `module_resolution=0.0009`）
10. 输出更新后的 `obj_with_pseudotime.rds`、`deg_graph_test.csv`、`pr_deg_ids.txt`、`gene_modules.csv`

额外还会保存：

- `gene_modules_prefilter.csv`：预模块划分结果
- `params_used_fig3_monocle3.yaml`
- `sessionInfo_fig3.txt`

## 常用运行方式

### 只跑 Monocle3

```bash
bash Fig3/run_fig3_combined.sh   --only monocle3   -i /path/to/obj_oo.rds   -o results/Fig3
```

### Monocle3 + Length + GO

```bash
bash Fig3/run_fig3_combined.sh   --only monocle3,length,go   -i /path/to/obj_oo.rds   -o results/Fig3
```

### 只跑 SCENIC

```bash
bash Fig3/run_fig3_combined.sh   --only scenic1,scenic2,scenic_downstream   -i results/Fig3/obj_with_pseudotime.rds   -o results/Fig3
```

### 只跑 scTour

```bash
bash Fig3/run_fig3_combined.sh   --only rds2h5ad,sctour   -i results/Fig3/obj_with_pseudotime.rds   -o results/Fig3
```

## 新增 Monocle3 相关配置

这些键已经加入 `fig3_combined.yaml`：

- `prefilter_enabled`
- `prefilter_module_resolution`
- `prefilter_target_module`
- `stage_col`
- `early_stage_levels`
- `late_stage_levels`
- `out_prefilter_modules_csv`

默认值对应你现在这版分析逻辑：

```yaml
prefilter_enabled: true
prefilter_module_resolution: 0.001
prefilter_target_module: '3'
stage_col: stage
early_stage_levels: [EGO, GO1]
late_stage_levels: [GO3, FGO]
out_prefilter_modules_csv: results/Fig3/gene_modules_prefilter.csv
```
