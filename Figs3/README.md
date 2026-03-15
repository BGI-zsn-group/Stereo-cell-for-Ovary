# Figs3 — Somatic QC → Harmony integration → annotation

本目录用于复现 **Supplementary Fig.S3（Figs3）**。  
现在的默认流程已经调整为：

1. **sc_qc / somatic QC**
2. **somatic processing（Harmony + 两轮清理 + 注释）**

也就是说，推荐入口 `run_figs3_somatic_processing_combined.sh` 会先跑 QC，再把 QC 输出自动作为下游整合脚本的输入。

---

## 目录结构

```text
Figs3/
  Copy-sc_qc.r                    # 你原始的 sc_qc 脚本（参考）
  figs3_somatic_qc.R              # repo 中实际调用的 QC 脚本
  figs3_somatic_processing.R      # 下游整合/聚类/注释
  run_figs3_somatic_processing_combined.sh
  configs/
    figs3_combined.yaml
```

---

## 推荐运行方式

从仓库根目录执行：

```bash
bash Figs3/run_figs3_somatic_processing_combined.sh
```

默认等价于：

```bash
bash Figs3/run_figs3_somatic_processing_combined.sh --only all
```

执行顺序为：

```text
figs3_qc  ->  figs3_somatic_processing
```

---

## 常用命令

### 1) 只跑 QC
```bash
bash Figs3/run_figs3_somatic_processing_combined.sh --only qc
```

### 2) 只跑后处理
适用于你已经有 QC 输出时：

```bash
bash Figs3/run_figs3_somatic_processing_combined.sh --only somatic_processing
```

### 3) 改输出根目录
```bash
bash Figs3/run_figs3_somatic_processing_combined.sh -o results/Figs3_alt
```

这会自动同步为：

- QC 输出：`results/Figs3_alt/qc/`
- 后处理输出：`results/Figs3_alt/somatic_processing/`
- 同时把 `figs3_somatic_processing.input_rds` 指到新的 QC 输出 RDS

### 4) 精确改配置项
```bash
bash Figs3/run_figs3_somatic_processing_combined.sh   --set figs3_qc.qc_input_rds=/path/to/raw_merged_somatic.rds   --set figs3_somatic_processing.out_dir=results/Figs3_custom/somatic_processing
```

---

## YAML 逻辑（关键更新）

`configs/figs3_combined.yaml` 现在已经改成默认串联：

### QC 模块
```yaml
modules:
  figs3_qc:
    qc_input_rds: data/Figs3/qc/somatic_seurat_object.rds
    qc_out_rds: results/Figs3/qc/somatic.qc_processed.rds
```

### Processing 模块
```yaml
modules:
  figs3_somatic_processing:
    input_rds: results/Figs3/qc/somatic.qc_processed.rds
```

也就是说，**processing 默认直接吃 QC 的输出**。

---

## 关于 `sample_map`

由于 QC 后输入通常是一个合并后的 Seurat 对象，`figs3_somatic_processing.R` 现在支持：

- 直接读取 `input_rds`
- 如果对象里还没有合适的 `sample` 列，则用 `sample_map` 把 `orig.ident` 映射到时间点标签

示例：

```yaml
sample_map:
  PD14-1: P14
  PD14-2: P14
  PD21-1: P21
  PD21-3: P21
  PD49-1: P49
  PD49-2: P49
```

---

## 主要输出

### QC
默认在 `results/Figs3/qc/`：
- `somatic.qc_processed.rds`
- `somatic.Raw.QC.VlnPlot.pdf`
- `somatic.MTQCed.VlnPlot.pdf`
- `somatic.UMIQCed.VlnPlot.pdf`
- `somatic.Final.QC.VlnPlot.pdf`
- `somatic.qc.log`

### Somatic processing
默认在 `results/Figs3/somatic_processing/`：
- `somatic.round1.rds`
- `somatic.final.rds`
- `somatic.round1.umap.pdf`
- `somatic.round2.umap.pdf`
- `somatic.final.umap.pdf`
- `somatic.cell_counts_by_sample_celltype.csv`
- `somatic.processing.log`

---

## 这次同步更新了哪些文件

本次已经按“先 QC，再处理”的逻辑同步更新：

1. `run_figs3_somatic_processing_combined.sh`
2. `figs3_combined.yaml`
3. `figs3_somatic_processing.R`
4. `README.md`
