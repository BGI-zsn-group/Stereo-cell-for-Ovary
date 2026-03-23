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

```text
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

### R（用于 `gc / ssgsea / monocle3 / c2l_sc`）
- R ≥ 4.2
- R packages：`Seurat`, `harmony`, `yaml`, `dplyr`, `ggplot2`, `patchwork`
- `GSVA`, `msigdbr`, `pheatmap`
- Monocle3（按你的安装方式：CRAN/Conda/GitHub）

### Python（用于 `c2l_spatial / c2l_train`）
- Python ≥ 3.9
- 建议用 conda 环境
- 依赖：`anndata`, `scanpy`, `numpy`, `pandas`, `scikit-learn`, `cell2location`, `scvi-tools`, `torch`

## 最常用跑法

从仓库根目录执行：

```bash
bash Fig4/run_fig4_combined.sh
```

### 只跑某一步

```bash
bash Fig4/run_fig4_combined.sh --only gc
bash Fig4/run_fig4_combined.sh --only ssgsea
bash Fig4/run_fig4_combined.sh --only monocle3
bash Fig4/run_fig4_combined.sh --only c2l_sc
bash Fig4/run_fig4_combined.sh --only c2l_spatial
bash Fig4/run_fig4_combined.sh --only c2l_train
```

### 只跑中间几步，并给一个自定义输入

对于 `ssgsea / monocle3 / c2l_sc` 这条中间链，它们都可以共享同一个 GC refined RDS 作为输入：

```bash
bash Fig4/run_fig4_combined.sh --only ssgsea,monocle3,c2l_sc   -i results/Fig4/gc_processing/obj_gr_newanno.rds
```

这里的 `-i` 是快捷主输入，会自动喂给：
- `fig4_ssgsea_hallmark.input_rds`
- `fig4_monocle3.input_rds`
- `fig4_cell2location_sc_prep.obj_gr_rds`

### 单步自定义输入

```bash
bash Fig4/run_fig4_combined.sh --only gc -i data/Fig4/obj_all_withanno.rds
bash Fig4/run_fig4_combined.sh --only c2l_spatial -i data/Fig4/spatial/B04372C211.adjusted.cellbin.h5ad
```

### 复杂情况：命名输入覆盖

当一个步骤不止一个输入，或者你想显式指定不同模块的输入时，用 `-i alias=PATH`：

```bash
bash Fig4/run_fig4_combined.sh --only c2l_train   -i train_sc_h5ad=results/Fig4/cell2location_sc_ref/somatic_1226.h5ad   -i train_spatial_h5ad=results/Fig4/cell2location_spatial/B04372C211.cell2location_spatial_counts.h5ad
```

支持的 alias：

- `gc_rds`
- `gc_result_rds`（同时作用于 `ssgsea / monocle3 / c2l_sc`）
- `ssgsea_rds`
- `monocle3_rds`
- `c2l_sc_gc_rds`
- `c2l_sc_total_rds`
- `spatial_h5ad`
- `train_sc_h5ad`
- `train_spatial_h5ad`

### 重定向输出目录（推荐）

`-o/--out` 传入一个 **base dir**，脚本会自动把每一步输出写到子目录，并且把下游步骤输入自动改到这些新输出上：

```bash
bash Fig4/run_fig4_combined.sh -o results/Fig4_alt
```

将会生成：

- `results/Fig4_alt/gc_processing/`
- `results/Fig4_alt/ssgsea_hallmark/`
- `results/Fig4_alt/monocle3/`
- `results/Fig4_alt/cell2location_sc_ref/`
- `results/Fig4_alt/cell2location_spatial/`
- `results/Fig4_alt/cell2location_train/`

### `-o` 和 `-i` 一起用

优先级是：

1. YAML 默认输入  
2. `-o` 自动串联下游输入  
3. `-i` 显式输入覆盖（最高优先级）

也就是说，下面这种是允许的：

```bash
bash Fig4/run_fig4_combined.sh --only ssgsea,monocle3,c2l_sc   -o results/Fig4_alt   -i /path/to/another_obj_gr_newanno.rds
```

这时输出仍写到 `results/Fig4_alt/...`，但中间步骤会用你手动指定的 RDS。

### 覆盖配置（`--set`）

- 支持嵌套 key：`pipeline_round2.cluster_resolution=1.2`
- 支持 module-scoped：`fig4_gc_processing.pipeline_round2.theta=2.0`（只影响该 module）

示例：

```bash
bash Fig4/run_fig4_combined.sh   --set fig4_gc_processing.pipeline_round2.cluster_resolution=1.2
```

## 配置文件说明（`configs/fig4_combined.yaml`）

发布到 GitHub 后，你通常需要改的只有这些默认输入路径：

- `modules.fig4_gc_processing.input_rds`
- `modules.fig4_cell2location_spatial_prep.input_h5ad`
- `modules.fig4_cell2location_sc_prep.obj_total_rds`

如果你不想改 YAML，也可以直接在命令行用 `-i` 覆盖。

## 输出概览

- GC processing：`results/Fig4/gc_processing/obj_gr_newanno.rds` 等
- ssGSEA：`results/Fig4/ssgsea_hallmark/`
- Monocle3：`results/Fig4/monocle3/gc.obj_with_pseudotime.rds` 等
- Cell2location sc ref：`results/Fig4/cell2location_sc_ref/*.h5ad`
- Spatial prep：`results/Fig4/cell2location_spatial/*.h5ad`
- Train：`results/Fig4/cell2location_train/` 下的模型与后验导出文件
