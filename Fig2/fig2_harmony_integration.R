#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Fig2 | Harmony integration (two-pass SCTransform + Harmony)
#
# What this script does
#   - Reads per-sample Seurat objects (*.rds) from a directory.
#   - Adds a `sample` metadata column based on file name.
#   - Two-pass workflow:
#       Pass 1) per-sample SCTransform -> merge -> SCTransform(residual.features)
#               -> PCA -> Harmony -> neighbors/clusters/UMAP -> remove clusters
#       Pass 2) split by sample -> per-sample SCTransform -> SCTransform(residual.features)
#               -> PCA -> Harmony(theta) -> neighbors/clusters/UMAP
#   - Saves the integrated Seurat object as an RDS.
#
# Recommended way to run (repo wrapper)
#   bash Fig2/run_fig2_combined.sh -i <rds_dir> -o <out_dir>
#
# Run this script directly (module-level YAML)
#   Rscript Fig2/fig2_harmony_integration.R --config <fig2_harmony.yaml>
#   # NOTE: This script expects a *module* YAML (keys like rds_dir/out).
#         If your repo uses a combined YAML (figure/modules/...), use the wrapper
#         which extracts the module config automatically.
#
# Outputs
#   - out (RDS): integrated Seurat object (contains `harmony` reduction, UMAP, clusters)
#   - params_used_fig2_harmony.yaml: parameters actually used (next to `out`)
#   - sessionInfo_fig2.txt: R session/package versions (next to `out`)
#
# Dependencies
#   Seurat, harmony, dplyr, patchwork, ggplot2, ggrepel, yaml
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(harmony)
  library(ggplot2)
  library(ggrepel)
})

# ---------- CLI ----------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    k <- args[[i]]
    if (!startsWith(k, "--")) stop("Unknown arg: ", k)
    key <- sub("^--", "", k)
    if (key %in% c("help")) {
      out[[key]] <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) stop("Missing value for ", k)
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

help_msg <- function() {
  cat(
"fig2_harmony_integration.R

Primary (recommended):
  --config <yaml>          YAML config path

Alternative:
  --rds-dir <dir>          Directory containing per-sample .rds
  --out <path>             Output merged Seurat RDS path

Example:
  Rscript fig2_harmony_integration.R --config configs/fig2_harmony.yaml
", sep = "")
}

# ---------- helpers ----------
need_yaml <- function() {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for --config mode. Please run: install.packages('yaml')", call. = FALSE)
  }
}

as_vars <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.character(x) && length(x) == 1) {
    v <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
    return(v[nzchar(v)])
  }
  # already a vector/list from YAML
  v <- unlist(x)
  v <- trimws(as.character(v))
  v[nzchar(v)]
}

as_int_vec <- function(x) as.integer(as_vars(x))

parse_dims <- function(x, default_seq) {
  if (is.null(x)) return(default_seq)
  if (is.numeric(x)) return(as.integer(x))
  x <- as.character(x)
  if (grepl(":", x, fixed = TRUE)) {
    p <- as.integer(trimws(unlist(strsplit(x, ":", fixed = TRUE))))
    return(seq(p[1], p[2]))
  } else {
    return(as.integer(as_vars(x)))
  }
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

# ---------- config loading ----------
args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }

cfg <- list()

if (!is.null(args$config)) {
  need_yaml()
  cfg <- yaml::read_yaml(args$config)
}

# CLI overrides (if provided)
if (!is.null(args$`rds-dir`)) cfg$rds_dir <- args$`rds-dir`
if (!is.null(args$out)) cfg$out <- args$out

# Required
rds_dir <- cfg$rds_dir
out_path <- cfg$out
if (is.null(rds_dir) || is.null(out_path)) {
  help_msg()
  stop("Missing required rds_dir/out in config (or via --rds-dir/--out).", call. = FALSE)
}

# Defaults
pattern <- if (!is.null(cfg$pattern)) cfg$pattern else "\\.rds$"
seed <- if (!is.null(cfg$seed)) as.integer(cfg$seed) else 1
nfeatures_integrate <- if (!is.null(cfg$nfeatures_integrate)) as.integer(cfg$nfeatures_integrate) else 1000
sct_vfeatures_n <- if (!is.null(cfg$sct_vfeatures_n)) as.integer(cfg$sct_vfeatures_n) else 2000
exclude_regex <- if (!is.null(cfg$exclude_regex)) cfg$exclude_regex else "^(Rp|mt)"
exclude_genes <- as_vars(cfg$exclude_genes)

dims_pass1 <- parse_dims(cfg$dims_pass1, 1:20)
res_pass1 <- if (!is.null(cfg$res_pass1)) as.numeric(cfg$res_pass1) else 1
remove_clusters <- if (!is.null(cfg$remove_clusters)) as_int_vec(cfg$remove_clusters) else c(8, 9)

dims_pass2 <- parse_dims(cfg$dims_pass2, 1:30)
res_pass2 <- if (!is.null(cfg$res_pass2)) as.numeric(cfg$res_pass2) else 1.4
theta_pass2 <- if (!is.null(cfg$theta_pass2)) as.numeric(cfg$theta_pass2) else 1

harmony_vars_pass1 <- if (!is.null(cfg$harmony_vars_pass1)) as_vars(cfg$harmony_vars_pass1) else c("sample")
harmony_vars_pass2 <- if (!is.null(cfg$harmony_vars_pass2)) as_vars(cfg$harmony_vars_pass2) else c("sample")

# Prepare output dir
out_dir <- dirname(out_path)
if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- run ----------
message("[fig2] Reading RDS from: ", rds_dir)
rds_files <- list.files(path = rds_dir, pattern = pattern, full.names = TRUE)
if (length(rds_files) == 0) stop("[fig2] No RDS files found in: ", rds_dir)

obj_list <- lapply(rds_files, readRDS)
sample_ids <- tools::file_path_sans_ext(basename(rds_files))
names(obj_list) <- sample_ids

# Ensure each object has 'sample'
obj_list <- mapply(function(obj, sid) {
  obj$sample <- sid
  obj
}, obj_list, sample_ids, SIMPLIFY = FALSE)

message("[fig2] Loaded samples: ", paste(sample_ids, collapse = ", "))
set.seed(seed)

# Pass 1: per-sample SCT -> merge -> SCT(residual.features) -> Harmony -> cluster/UMAP
message("[fig2] Pass1: per-sample SCTransform ...")
obj_list <- lapply(obj_list, function(x) SCTransform(x, variable.features.n = sct_vfeatures_n, verbose = FALSE))

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures_integrate)

# Use a reference gene list (obj not defined yet)
ref_genes <- rownames(obj_list[[1]])
fl_genes <- grep(exclude_regex, ref_genes, value = TRUE)

features <- setdiff(features, fl_genes)
if (length(exclude_genes) > 0) features <- setdiff(features, exclude_genes)
message("[fig2] Pass1: features kept = ", length(features))

message("[fig2] Merging ...")
obj <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list), project = "Fig2")
DefaultAssay(obj) <- "RNA"

message("[fig2] Pass1: SCTransform on merged (residual.features) ...")
obj <- SCTransform(obj, residual.features = features, verbose = FALSE)

obj <- RunPCA(obj, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = harmony_vars_pass1, reduction = "pca", reduction.save = "harmony")

obj <- FindNeighbors(obj, reduction = "harmony", dims = dims_pass1)
obj <- FindClusters(obj, reduction = "harmony", resolution = res_pass1)
obj <- RunUMAP(obj, reduction = "harmony", dims = dims_pass1)

# Remove contamination clusters
if (length(remove_clusters) > 0) {
  message("[fig2] Removing clusters after pass1: ", paste(remove_clusters, collapse = ", "))
  obj <- subset(obj, subset = !(seurat_clusters %in% remove_clusters))
}

# Pass 2: re-run after filtering
message("[fig2] Pass2: split by sample and re-run SCTransform ...")
DefaultAssay(obj) <- "RNA"
obj_list2 <- SplitObject(obj, split.by = "sample")

obj_list2 <- lapply(obj_list2, function(x) SCTransform(x, variable.features.n = sct_vfeatures_n, verbose = FALSE))
features2 <- SelectIntegrationFeatures(object.list = obj_list2, nfeatures = nfeatures_integrate)

fl_genes2 <- grep(exclude_regex, rownames(obj), value = TRUE)
features2 <- setdiff(features2, fl_genes2)
if (length(exclude_genes) > 0) features2 <- setdiff(features2, exclude_genes)
message("[fig2] Pass2: features kept = ", length(features2))

message("[fig2] Pass2: SCTransform on filtered merged (residual.features) ...")
obj <- SCTransform(obj, residual.features = features2, verbose = FALSE)

obj <- RunPCA(obj, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = harmony_vars_pass2, reduction = "pca", reduction.save = "harmony", theta = theta_pass2)

obj <- FindNeighbors(obj, reduction = "harmony", dims = dims_pass2)
obj <- FindClusters(obj, reduction = "harmony", resolution = res_pass2)
obj <- RunUMAP(obj, reduction = "harmony", dims = dims_pass2)

# Save results
saveRDS(obj, out_path, compress = TRUE)
message("[fig2] Saved object: ", out_path)

# Write reproducibility artifacts
params_used <- list(
  rds_dir = rds_dir,
  out = out_path,
  pattern = pattern,
  seed = seed,
  nfeatures_integrate = nfeatures_integrate,
  sct_vfeatures_n = sct_vfeatures_n,
  exclude_regex = exclude_regex,
  exclude_genes = exclude_genes,
  dims_pass1 = paste0(min(dims_pass1), ":", max(dims_pass1)),
  res_pass1 = res_pass1,
  remove_clusters = remove_clusters,
  dims_pass2 = paste0(min(dims_pass2), ":", max(dims_pass2)),
  res_pass2 = res_pass2,
  theta_pass2 = theta_pass2,
  harmony_vars_pass1 = harmony_vars_pass1,
  harmony_vars_pass2 = harmony_vars_pass2
)

if (!requireNamespace("yaml", quietly = TRUE)) {
  message("[fig2] NOTE: package 'yaml' not installed; skipping params_used YAML write.")
} else {
  yaml_path <- file.path(out_dir, "params_used_fig2_harmony.yaml")
  yaml::write_yaml(params_used, yaml_path)
  message("[fig2] Saved params: ", yaml_path)
}

si_path <- file.path(out_dir, "sessionInfo_fig2.txt")
write_text(si_path, capture.output(sessionInfo()))
message("[fig2] Saved sessionInfo: ", si_path)

message("[fig2] Done.")
