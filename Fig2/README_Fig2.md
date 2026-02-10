# Fig2 — SCTransform + Harmony 多样本整合（两轮流程）

本目录用于复现论文 **Fig.2** 的整合分析：将多个样本的 Seurat 对象（`.rds`）进行 **SCTransform + Harmony** 整合，并输出一个最终的整合 Seurat 对象（含 Harmony reduction、UMAP、聚类结果等）。

## 目录结构

```
Fig2/
  configs/
    fig2_combined.yaml
  fig2_harmony_integration.R
  run_fig2_combined.sh
```

- `run_fig2_combined.sh`：**推荐入口**。从 combined YAML 中提取 `modules.fig2_harmony` 这一段配置（module），然后调用 R 脚本运行。
- `fig2_harmony_integration.R`：核心 R pipeline（也支持直接运行，但要求传入“module-level”的 YAML，见下文）。
- `configs/fig2_combined.yaml`：combined 配置文件（包含 `figure:` 与 `modules:`；本目录使用的模块为 `fig2_harmony`）。

## 输入（Inputs）

- 一个目录，里面包含 **每个样本一个** 的 Seurat 对象 `.rds` 文件。
- 每个 `.rds` 应为标准 Seurat object，至少包含 `RNA` assay。
- 脚本会用文件名（去掉扩展名）作为样本名，并写入元数据列 `sample`。

## 输出（Outputs）

- `out`：整合后的 Seurat object（RDS）。
- 额外写出两个可复现文件（与 `out` 同目录）：
  - `params_used_fig2_harmony.yaml`：本次运行实际使用的参数
  - `sessionInfo_fig2.txt`：R 运行环境与包版本（`sessionInfo()`）

输出 Seurat object 通常包含：
- Harmony reduction（名字为 `harmony`）
- UMAP（名字为 `umap`）
- 聚类标签 `seurat_clusters`

## 快速开始（推荐）

在仓库根目录运行：

```bash
bash Fig2/run_fig2_combined.sh \
  -i /path/to/per_sample_rds_dir \
  -o /path/to/output_dir
```

说明：
- `-i/--in` 会写入配置项 `rds_dir`。
- `-o/--out` 是便捷参数：
  - 如果以 `.rds` 结尾：视为**输出文件路径**，直接写入 `out=/.../...rds`
  - 否则：视为**输出目录**，默认写入 `out=<out_dir>/obj_oo.rds`

你也可以用 `--set key=value` 覆盖任意配置项（可重复）：

```bash
bash Fig2/run_fig2_combined.sh \
  -i /path/to/per_sample_rds_dir \
  -o /path/to/output_dir \
  --set seed=1 \
  --set res_pass2=1.4 \
  --set theta_pass2=1
```

## 配置说明（`configs/fig2_combined.yaml`）

本目录使用 `modules.fig2_harmony` 下的参数。常用键如下：

| 键 | 含义 |
|---|---|
| `rds_dir` | 输入目录：每个样本一个 `.rds` Seurat object |
| `out` | 输出整合对象（RDS）路径 |
| `pattern` | 用于筛选输入文件的正则（默认 `\\.rds$`） |
| `seed` | 随机种子 |
| `sct_vfeatures_n` | per-sample `SCTransform(variable.features.n=...)` |
| `nfeatures_integrate` | `SelectIntegrationFeatures()` 选取 feature 数 |
| `exclude_regex` | 从 features 中排除的基因正则（默认 `^(Rp|mt)`） |
| `exclude_genes` | 额外手动排除的基因列表 |
| `dims_pass1` | Pass1 PCA/Harmony 使用的维度（如 `1:20`） |
| `res_pass1` | Pass1 聚类分辨率 |
| `remove_clusters` | Pass1 后要删除的 cluster id（用于去污染/异常簇） |
| `harmony_vars_pass1` | Pass1 Harmony 的 batch 变量（默认 `sample`） |
| `dims_pass2` | Pass2 PCA/Harmony 使用的维度（如 `1:30`） |
| `res_pass2` | Pass2 聚类分辨率 |
| `theta_pass2` | Pass2 Harmony 的 `theta` |
| `harmony_vars_pass2` | Pass2 Harmony 的 batch 变量（默认 `sample`） |

## 直接运行 R 脚本（高级用法）

`fig2_harmony_integration.R` 需要的是 **module-level YAML**（顶层包含 `rds_dir/out/...`），而不是 combined YAML（`figure/modules/...`）。

如果你不想用 bash wrapper，可以先把 module 配置抽出来：

```bash
python - <<'PY'
import yaml
root = yaml.safe_load(open('Fig2/configs/fig2_combined.yaml', 'r', encoding='utf-8'))
mod = root['modules']['fig2_harmony']
yaml.safe_dump(mod, open('Fig2/configs/fig2_harmony.yaml', 'w', encoding='utf-8'),
               sort_keys=False, allow_unicode=True)
print('wrote: Fig2/configs/fig2_harmony.yaml')
PY

Rscript Fig2/fig2_harmony_integration.R --config Fig2/configs/fig2_harmony.yaml
```

## 环境依赖（Software requirements）

- Bash
- Python 3（wrapper 读取 YAML 用）
  - `pyyaml`：`pip install pyyaml`
- R（建议 R >= 4.2）
  - R 包：`Seurat`, `harmony`, `dplyr`, `patchwork`, `ggplot2`, `ggrepel`, `yaml`

## 常见问题（Troubleshooting）

- **`No RDS files found`**：检查 `rds_dir` 与 `pattern`。
- **`missing pyyaml`**：在运行 wrapper 的 Python 环境里安装 `pyyaml`。
- **R 提示缺少 `yaml`**：`install.packages('yaml')`。
- Seurat/Harmony 在 merge 或 SCTransform 报错：通常是输入对象不规范或 assay/基因命名不一致；请确保每个输入 RDS 都是可用的 Seurat object，并且基因标识一致（同一参考基因集合）。

## 可复现性说明

- 脚本会在输出目录写出 `params_used_fig2_harmony.yaml` 与 `sessionInfo_fig2.txt`。
- 设置 `seed` 有助于尽量复现结果；但不同机器/BLAS/并行策略可能仍会带来轻微随机性差异。
