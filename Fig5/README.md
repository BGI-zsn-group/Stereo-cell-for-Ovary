# Fig5 — CellChat 细胞通讯分析（卵母细胞阶段 × 颗粒细胞状态）

本目录用于复现论文 **Fig5** 的 CellChat 分析流程：将 **颗粒细胞（GC）Seurat 对象** 与 **卵母细胞（oocyte）Seurat 对象** 合并后，
按不同卵母细胞阶段（EGO/GO1/GO2/GO3/FGO）与对应的 GC 亚型组合构建多个子对象，分别运行完整 CellChat pipeline，
并导出每个阶段的 CellChat 对象与配体-受体通讯表（带 interaction annotation）。

该流程同时支持两种“数据库模式”：
- `mode: all`：使用 `CellChatDB.mouse` 全库（包含 Secreted / ECM / Cell-Cell Contact 等类别）
- `mode: cellcell_contact`：只使用 `subsetDB(CellChatDB.mouse, search="Cell-Cell Contact")`（仅 Cell-Cell Contact）

---

## 目录结构

```
Fig5/
  configs/
    fig5_combined.yaml
  fig5_cellchat_build.R
  run_fig5_combined.sh
```

- `run_fig5_combined.sh`：**推荐入口**。从 combined YAML 中提取 `modules.fig5_cellchat_build`，再调用 R 脚本运行。
- `fig5_cellchat_build.R`：主流程脚本（合并对象 → 分阶段子集 → CellChat pipeline → liftCellChat → 导出）。
- `configs/fig5_combined.yaml`：配置文件（发布到 GitHub 后请把输入路径改为你的本地路径）。

---

## 依赖环境

- Bash（Linux/macOS 或 WSL）
- Python（仅用于解析 combined YAML）：需要 `pyyaml`
  - `pip install pyyaml`
- R（建议 ≥ 4.2）
  - 主要 R 包：`Seurat`, `CellChat`, `yaml`, `magrittr`, `ggplot2`, `patchwork`, `reshape2`, `grid`

> CellChat 会使用内置的 `PPI.mouse` 网络（脚本中 `smoothData(adj = PPI.mouse)`）。

---

## 运行方式

从仓库根目录执行：

```bash
bash Fig5/run_fig5_combined.sh -i Fig5/configs/fig5_combined.yaml
```

### 常用参数

- `-i/--in`：combined YAML 路径（默认 `Fig5/configs/fig5_combined.yaml`）
- `-o/--out`：覆盖输出目录；或传入 `xxx.rds` 形式用于推导 `{out_dir, prefix}`
- `--set key=value`：覆盖 module 配置（可重复；支持点号路径）

示例：

```bash
# 只跑 Cell-Cell Contact（对应 CCC-only 分析）
bash Fig5/run_fig5_combined.sh --set mode=cellcell_contact

# 覆盖输入对象路径
bash Fig5/run_fig5_combined.sh   --set obj_gr_rds=data/Fig5/obj_gr_newanno.rds   --set obj_oo_rds=data/Fig5/oocyte_edit.rds

# 覆盖输出目录
bash Fig5/run_fig5_combined.sh -o results/Fig5/cellchat_debug
```

---

## 配置说明（`modules.fig5_cellchat_build`）

你最需要改的是输入路径：

- `obj_gr_rds`：GC Seurat 对象（需要包含 `seurat_clusters`、`sample` 等 meta 信息）
- `obj_oo_rds`：oocyte Seurat 对象（需要包含 `stage` 等 meta 信息）

### 关键字段

- `mode`：
  - `all`：全库
  - `cellcell_contact`：仅 Cell-Cell Contact 子库
- `cluster_to_celltype`：**可选**。用于将 GC 的 `seurat_clusters` 映射为 `celltype`
  - 脚本会按该映射重写/生成 `gc_celltype_col`（默认 `celltype`）
- `gc_exclude_celltypes`：默认剔除 `GC_Atretic`
- `oocyte_exclude_stages`：默认剔除 `MII`
- `group_new_levels`：用于 `liftCellChat()` 的统一分组顺序（保证不同阶段对象的 group 对齐）
- `compute_type / seed_use / min_cells`：CellChat 计算参数（默认 `triMean`, `seed=888`, `min_cells=10`）

### 阶段子集定义（高级）

脚本支持自定义阶段子集 `stage_subsets`（YAML 中可选填写）；若不填写，会根据 `mode` 使用内置默认规则：
- `mode=all`：`GO3/FGO` 默认包含 `GC_Antral_Mural`
- `mode=cellcell_contact`：`GO3/FGO` 默认只保留 `GC_Cumulus_2`（与你原始 CCC-only 代码一致）

---

## 输出文件

全部输出位于 `out_dir`（默认 `results/Fig5/cellchat/`）：

- `cellchat_<STAGE>_triMean.rds` 或 `cellchat_<STAGE>_triMean_CCC.rds`  
  例如：`cellchat_EGO_triMean.rds`, `cellchat_GO1_triMean_CCC.rds`
- `comm_<STAGE>_triMean.csv` 或 `comm_<STAGE>_triMean_CCC.csv`（可选，默认导出）
  - 该表来自 `subsetCommunication()`，并合并了 `interaction_name` 对应的 `annotation`（通讯类别）
- `<prefix>.cellchat_all_<tag>.rds`：一个汇总 RDS（包含各阶段的 CellChat 对象与通讯表）
- `<prefix>.cellchat_build.log`：运行日志

---

## 常见问题

1. **输入对象缺少 meta 列**
   - GC 需要：`seurat_clusters`（或你在 `gc_cluster_col` 指定的列）、`sample`（或 `gc_sample_col`）
   - oocyte 需要：`stage`（或你在 `oo_stage_col` 指定的列）

2. **跑得很慢 / 内存占用大**
   - CellChat 对数据规模较敏感，建议先按目标细胞类型/阶段完成过滤。
   - 可以提高 `min_cells` 或减少阶段子集中的 group 数量以加速。

3. **结果中某些 group 没有通讯**
   - 常见原因：该 group 的细胞数不足（`min_cells`）、或数据库子集（CCC-only）不包含该类互作。

---
