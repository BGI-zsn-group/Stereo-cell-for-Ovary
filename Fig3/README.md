# Fig3 — Pseudotime / Gene modules / GO / SCENIC / scTour（整合脚本）

本目录用于复现论文 **Fig.3** 的一整套分析流程，核心围绕：
- **Monocle3**：从整合后的 Seurat 对象推断轨迹与伪时序（pseudotime），并识别轨迹相关基因与模块；
- **Length（Diameter）~ Pseudotime**：平滑拟合与相关性计算；
- **GO 富集**：对合并后的模块做 GO 富集并简化冗余条目；
- **SCENIC / pySCENIC**：转录因子调控网络推断（表达矩阵导出 → pySCENIC → RSS downstream）；
- **scTour**：TNODE 伪时序与潜变量（以 `.h5ad` 为输入）。

> 说明：Fig3 默认输入是 Fig2 的整合 Seurat 对象（例如 `obj_oo.rds`）。你也可以用 `-i` 显式指定输入 RDS。

---

## 目录结构

```
Fig3/
  configs/
    fig3_combined.yaml
  run_fig3_combined.sh
  fig3_monocle3_modules.R
  fig3_length_pseudotime.R
  fig3_go_enrichment.R
  fig3_scenic_rds_to_csv.R
  fig3_scenic_pyscenic.py
  fig3_scenic_downstream.R
  fig3_rds_to_h5ad.R
  fig3_sctour_tnode.py
```

- `run_fig3_combined.sh`：**推荐入口**。从 combined YAML 中提取 `modules.<module-key>` 这一段配置（默认 `module-key=fig3`），再按步骤调用对应 R / Python 脚本。
- `configs/fig3_combined.yaml`：combined 配置文件，包含多个模块（fig3、fig3_monocle3_with_go、…）。一般用默认 `fig3` 即可。

---

## 输入与输出（Overview）

### 主要输入
- `input_rds`：整合后的 Seurat 对象（建议来自 Fig2 输出）。
- （SCENIC）pySCENIC 资源文件：
  - `scenic_tf_list`
  - `scenic_rankings_feather`
  - `scenic_motif_annotations`
- （scTour）`sctour_input_h5ad`：scTour 输入 AnnData 文件（可由 `fig3_rds_to_h5ad.R` 生成）

### 主要输出（默认写入 `results/Fig3/`）
- `obj_with_pseudotime.rds`：包含 pseudotime + pseudotime bins 的 Seurat 对象
- `deg_graph_test.csv`：Monocle3 graph_test 结果
- `gene_modules.csv`：基因模块归属
- `length_vs_pseudotime.pdf/.png` + `length_vs_pseudotime_correlation.txt`
- `GO_BP_merged_modules_simplified.csv`
- `results/Fig3/scenic/`：SCENIC 输出（loom/adj/reg/result loom 以及 downstream RSS）
- `results/Fig3/sctour/`：scTour 输出（`oocyte_sctour_tnode.h5ad`, `ptime.csv` 等）

---

## 快速开始（推荐）

在仓库根目录运行：

```bash
# 1) 全流程（包含 Monocle3/Length/GO/SCENIC/scTour）
bash Fig3/run_fig3_combined.sh   -i /path/to/obj_oo.rds   -o results/Fig3
```

你也可以只跑部分步骤（逗号分隔）：

```bash
# 只跑 Monocle3 + Length + GO
bash Fig3/run_fig3_combined.sh   --only monocle3,length,go   -i /path/to/obj_oo.rds   -o results/Fig3
```

通过 `--set key=value` 覆盖任意配置项（可重复）：

```bash
bash Fig3/run_fig3_combined.sh   --only monocle3   -i /path/to/obj_oo.rds   -o results/Fig3   --set root_cluster=6   --set pseudotime_bins=13
```

---

## 步骤说明（Steps）

脚本支持的步骤列表：

- `monocle3` → `fig3_monocle3_modules.R`
- `length` → `fig3_length_pseudotime.R`
- `go` → `fig3_go_enrichment.R`
- `scenic1` → `fig3_scenic_rds_to_csv.R`（Seurat → CSV）
- `scenic2` → `fig3_scenic_pyscenic.py`（pySCENIC：grn/ctx/aucell）
- `scenic_downstream` → `fig3_scenic_downstream.R`（RSS 计算与绘图）
- `rds2h5ad` → `fig3_rds_to_h5ad.R`（Seurat → h5ad）
- `sctour` → `fig3_sctour_tnode.py`

> 提示：SCENIC 的 `scenic2` 依赖 `scenic1` 输出的表达矩阵 CSV；如果你已经有现成 CSV，只需在 YAML 中设置 `scenic_out_csv` 指向它，然后直接跑 `scenic2` 即可。

---

## 配置说明（`configs/fig3_combined.yaml`）

`run_fig3_combined.sh` 默认会抽取 `modules.fig3` 这一段作为 module-level config。常用键：

### Monocle3 / Modules
| 键 | 含义 |
|---|---|
| `input_rds` | 输入 Seurat RDS（Fig2 输出） |
| `out_dir` | Fig3 输出目录 |
| `out_obj_rds` | 输出 pseudotime Seurat RDS |
| `root_cluster` | 轨迹 root cluster（决定 pseudotime 起点） |
| `pseudotime_bins` | 将 pseudotime 离散分箱数量 |
| `out_deg_csv` / `out_gene_modules_csv` | 中间文件输出路径 |

### Length vs Pseudotime
| 键 | 含义 |
|---|---|
| `plot_out_pdf` / `plot_out_png` | 输出图 |
| `plot_out_cor_txt` | 相关性输出 |
| `plot_loess_span` | LOESS 平滑参数 |

### GO enrichment
| 键 | 含义 |
|---|---|
| `go_merge_map` | 模块合并映射（如 M1 由哪些 module 组成） |
| `go_out_csv` | GO 输出表 |

### SCENIC
| 键 | 含义 |
|---|---|
| `scenic_out_csv` | SCENIC 输入表达矩阵（cells x genes） |
| `scenic_tf_list` | TF 列表 |
| `scenic_rankings_feather` | rankings feather |
| `scenic_motif_annotations` | motif 注释表 |
| `scenic_result_loom` | pySCENIC 输出 loom（供 downstream 使用） |
| `scenic_rss_var` | downstream 分组变量（默认 `pseudotime_bin_q`） |

### scTour
| 键 | 含义 |
|---|---|
| `rds2h5ad_output_h5ad` | rds→h5ad 输出 |
| `sctour_input_h5ad` | scTour 输入 h5ad |
| `sctour_output_h5ad` / `sctour_ptime_csv` | scTour 输出 |

---

## 环境依赖（Software requirements）

### Bash wrapper
- Bash
- Python 3 + `pyyaml`（wrapper 解析 YAML 用）  
  `pip install pyyaml`

### R（建议 R >= 4.2）
- Monocle3 / Seurat 相关：`Seurat`, `SeuratWrappers`, `monocle3`, `dplyr`, `ggplot2`, `patchwork`, `yaml`
- GO 富集：`clusterProfiler`, `org.Mm.eg.db`（Bioconductor）
- SCENIC downstream：`SCopeLoomR`, `AUCell`, `SCENIC`

### Python（SCENIC / scTour）
- SCENIC step2：`scanpy`, `loompy`, `numpy`, `pyyaml`，并且需要 `pyscenic` CLI 在 PATH 中
- scTour：`sctour`, `scanpy`, `numpy`, `pandas`, `scipy`, `pyyaml`

---

## 常见问题（Troubleshooting）

- **`override must be key=value` / `unbound variable`**：请使用本仓库提供的 `run_fig3_combined.sh`（已做 `set -u` 兼容与参数校验）。
- **SCENIC 资源文件缺失**：在 YAML 中把 `scenic_tf_list / scenic_rankings_feather / scenic_motif_annotations` 指到你本地下载的资源文件位置。
- **sceasy 转换失败**：通常是 reticulate/python 环境问题；确保你的 R 能找到一个装有 `anndata` 的 Python 环境。

---

## 可复现性说明

- R 脚本会在输出目录写出 `params_used_*.yaml` 与 `sessionInfo_*.txt`（如果该脚本实现了这部分）。
- 设置 `seed`（若对应脚本支持）可以尽量减少随机性差异；不同机器/并行策略可能仍会带来轻微差异。
