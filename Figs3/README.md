# Figs3 — Somatic 细胞整合与注释（Harmony + 两轮清理）

本目录用于复现论文 **Supplementary Fig.S3（Figs3）**：将多个时间点/批次的 QC 后 Seurat 对象合并，进行 **SCTransform/PCA + Harmony** 整合、聚类、两轮剔除低质量簇，并根据 `celltype_map` 完成细胞类型注释，输出用于作图的 RDS / PDF / 表格。

## 目录结构

```
Figs3/
  configs/
    figs3_combined.yaml
  figs3_somatic_processing.R
  run_figs3_somatic_processing_combined.sh
```

- `run_figs3_somatic_processing_combined.sh`：**推荐入口**。从 combined YAML 中提取 `modules.figs3_somatic_processing` 这一段（module），再调用 R 脚本运行。
- `figs3_somatic_processing.R`：主流程脚本（整合 / 聚类 / 两轮清理 / 注释 / 导出）。
- `configs/figs3_combined.yaml`：配置文件（注意：发布到 GitHub 后需要你把输入 `.rds` 路径改成自己的本地路径）。

## 依赖环境

- Bash（Linux/macOS 或 WSL）
- Python（仅用于解析 combined YAML）：需要 `pyyaml`
  - `pip install pyyaml`
- R（建议 ≥ 4.2）
  - R 包：`Seurat`, `harmony`, `yaml`, `dplyr`, `ggplot2`, `patchwork`

## 运行方式

从仓库根目录执行：

```bash
bash Figs3/run_figs3_somatic_processing_combined.sh \
  -i Figs3/configs/figs3_combined.yaml
```

### 常用参数

- `-i/--in`：combined YAML 路径（默认 `Figs3/configs/figs3_combined.yaml`）
- `-o/--out`：覆盖输出路径
  - 传 **目录**：所有输出写到该目录
  - 传 **以 `.rds` 结尾的文件路径**：作为最终 `out_final_rds`
- `--set key=value`：覆盖配置（可重复；支持点号路径写法）
  - 例：`--set pipeline_round2.theta=2.0`

示例：

```bash
# 改输出目录
bash Figs3/run_figs3_somatic_processing_combined.sh -o results/Figs3/somatic_processing_alt

# 覆盖二轮 Harmony theta
bash Figs3/run_figs3_somatic_processing_combined.sh --set pipeline_round2.theta=2.0
```

## 配置说明（`modules.figs3_somatic_processing`）

你最需要改的是输入文件路径：

- `groups[*].items[*].rds`：每个样本/重复的 QC 后 Seurat `.rds` 路径
- `exclude_regex`：在挑选高变基因时排除的基因前缀（默认 `^(Rp|mt)`）
- `pipeline_round1 / pipeline_round2`：两轮整合参数（PCA 维度、Harmony 参数、UMAP/邻居、聚类分辨率等）
- `exclude_clusters_round1 / exclude_clusters_round2`：两轮需要剔除的 cluster id 列表
- `celltype_map`：最终 cluster -> cell type 的映射（字符串形式的 cluster id）

## 输出文件

全部输出位于 `out_dir` 下（默认 `results/Figs3/somatic_processing/`）：

- `somatic.round1.rds`：第一轮整合/清理后的 Seurat 对象
- `somatic.final.rds`：最终 Seurat 对象（用于 Fig.S3 作图/统计）
- `somatic.round1.umap.pdf` / `somatic.round2.umap.pdf` / `somatic.final.umap.pdf`
- `somatic.cell_counts_by_sample_celltype.csv`：按 sample × celltype 的细胞计数表
- `somatic.processing.log`：运行日志

> 文件名前缀由配置 `prefix` 控制（默认 `somatic`）。

## 常见问题

1. **找不到输入 RDS**
   - 请检查 `configs/figs3_combined.yaml` 中 `groups[*].items[*].rds` 是否已改成你本地实际路径。

2. **缺少 R 包**
   - 在 R 中安装：`install.packages(c("yaml","dplyr","ggplot2","patchwork"))`
   - Seurat/harmony：按你环境（CRAN/Conda）安装对应版本。

3. **PyYAML 缺失**
   - `pip install pyyaml`
