# Fig4 — GC processing + Monocle3 + Cell2location (sc reference + spatial mapping)

本目录用于复现论文 **Fig4** 的主要分析流程，按模块组织为 6 个步骤：

1. `gc`：从全量注释 Seurat 对象中提取 Granulosa cell (GC) 并进行两轮 Harmony 整合/清理
2. `ssgsea`：对 GC 各 celltype 的平均表达做 Hallmark ssGSEA（GSVA）并输出热图/矩阵
3. `monocle3`：在 GC 子集上进行 Monocle3 轨迹推断（pseudotime）
4. `c2l_sc`：构建 Cell2location 单细胞参考（Seurat → AnnData `.h5ad`）
5. `c2l_spatial`：空间转录组 `.h5ad` 预处理（生成可用于 cell2location 的 counts AnnData）
6. `c2l_train`：cell2location 训练：Regression (sc reference) + spatial mapping

> 推荐通过 `run_fig4_combined.sh` 从 combined YAML 中提取对应 module 配置并运行每一步。

## 目录结构

```
Fig4/
  README.md
  configs/
    fig4_combined.yaml
  run_fig4_combined.sh
  fig4_gc_processing.R
  fig4_monocle3.R
  fig4_ssgsea_hallmark.R
  fig4_cell2location_sc_prep.R
  fig4_cell2location_spatial_prep.py
  fig4_cell2location_train.py
```

## 依赖环境

### R（用于 `gc / monocle3 / c2l_sc`）
- R ≥ 4.2
- R packages：`Seurat`, `harmony`, `yaml`, `dplyr`, `ggplot2`, `patchwork`
- Monocle3（按你的安装方式：CRAN/Conda/GitHub）

### Python（用于 `c2l_spatial / c2l_train`）
- Python ≥ 3.9
- 建议用 conda 环境
- 依赖：`anndata`, `scanpy`, `numpy`, `pandas`, `scikit-learn`, `cell2location`, `scvi-tools`, `torch`
  - 具体版本请以你论文复现环境为准

## 运行方式

从仓库根目录执行：

```bash
bash Fig4/run_fig4_combined.sh
```

### 只跑某一步

```bash
bash Fig4/run_fig4_combined.sh --only gc
bash Fig4/run_fig4_combined.sh --only monocle3
bash Fig4/run_fig4_combined.sh --only c2l_sc
bash Fig4/run_fig4_combined.sh --only c2l_spatial
bash Fig4/run_fig4_combined.sh --only c2l_train
```

### 重定向输出目录（推荐用于复现实验分组）

`-o/--out` 传入一个 **base dir**，脚本会自动把每一步输出写到子目录，避免相互覆盖：

```bash
bash Fig4/run_fig4_combined.sh -o results/Fig4_alt
```

将会生成：
- `results/Fig4_alt/gc_processing/`
- `results/Fig4_alt/monocle3/`
- `results/Fig4_alt/cell2location_sc_ref/`
- `results/Fig4_alt/cell2location_spatial/`
- `results/Fig4_alt/cell2location_train/`

### 覆盖配置（--set）

- 支持嵌套 key：`pipeline_round2.cluster_resolution=1.2`
- 支持 module-scoped：`fig4_gc_processing.pipeline_round2.theta=2.0`（只影响该 module）

示例：

```bash
bash Fig4/run_fig4_combined.sh --set fig4_gc_processing.pipeline_round2.cluster_resolution=1.2
```

## 配置文件说明（`configs/fig4_combined.yaml`）

发布到 GitHub 后，你通常需要改的只有输入路径：

- `modules.fig4_gc_processing.input_rds`
- `modules.fig4_cell2location_spatial_prep.input_h5ad`

其余输出路径默认写到 `results/Fig4/...`，并且各步骤之间的依赖关系已在 YAML 中串起来（如 monocle3 使用 gc 的输出）。

## 输出概览

- GC processing：`results/Fig4/gc_processing/obj_gr_newanno.rds` 等
- Monocle3：`results/Fig4/monocle3/gc.obj_with_pseudotime.rds` 等
- Cell2location sc ref：`results/Fig4/cell2location_sc_ref/*.h5ad`
- Spatial prep：`results/Fig4/cell2location_spatial/*.h5ad`
- Train：`results/Fig4/cell2location_train/` 下的模型与后验导出文件（具体见脚本日志）
