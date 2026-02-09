#!/usr/bin/env Rscript
# Fig.3: Monocle3 pseudotime + graph_test DEG + gene modules (from Seurat object)
#
# Recommended:
#   Rscript fig3_monocle3_modules.R --config configs/fig3_monocle3.yaml
#
# Outputs (written under out_dir):
#   - obj_with_pseudotime.rds
#   - deg_graph_test.csv
#   - pr_deg_ids.txt
#   - gene_modules.csv
#   - params_used_fig3_monocle3.yaml
#   - sessionInfo_fig3.txt
#
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(ggrepel)
  library(tidyverse)
  library(monocle3)
  library(SeuratWrappers)
})

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
"fig3_monocle3_modules.R

Required:
  --config <yaml>          YAML config path

Example:
  Rscript fig3_monocle3_modules.R --config configs/fig3_monocle3.yaml
", sep = "")
}

need_yaml <- function() {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Please run: install.packages('yaml')", call. = FALSE)
  }
}

as_vars <- function(x) {
  if (is.null(x)) return(character(0))
  v <- unlist(x)
  v <- trimws(as.character(v))
  v[nzchar(v)]
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

need_yaml()
cfg <- yaml::read_yaml(args$config)

in_rds <- cfg$input_rds
if (is.null(in_rds)) stop("Missing input_rds in YAML", call. = FALSE)

out_dir <- if (!is.null(cfg$out_dir)) cfg$out_dir else "results/Fig3"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_obj_rds <- if (!is.null(cfg$out_obj_rds)) cfg$out_obj_rds else file.path(out_dir, "obj_with_pseudotime.rds")

num_dim <- if (!is.null(cfg$num_dim)) as.integer(cfg$num_dim) else 50
cores <- if (!is.null(cfg$cores)) as.integer(cfg$cores) else 8

root_cluster <- if (!is.null(cfg$root_cluster)) as.character(cfg$root_cluster) else "6"

q_value_thr <- if (!is.null(cfg$q_value_thr)) as.numeric(cfg$q_value_thr) else 0.05
morans_I_thr <- if (!is.null(cfg$morans_I_thr)) as.numeric(cfg$morans_I_thr) else 0.2

pseudotime_bins <- if (!is.null(cfg$pseudotime_bins)) as.integer(cfg$pseudotime_bins) else 13

module_resolution <- if (!is.null(cfg$module_resolution)) as.numeric(cfg$module_resolution) else 0.0009
module_seed <- if (!is.null(cfg$module_seed)) as.integer(cfg$module_seed) else 888

deg_csv <- if (!is.null(cfg$out_deg_csv)) cfg$out_deg_csv else file.path(out_dir, "deg_graph_test.csv")
pr_ids_txt <- if (!is.null(cfg$out_pr_deg_ids)) cfg$out_pr_deg_ids else file.path(out_dir, "pr_deg_ids.txt")
modules_csv <- if (!is.null(cfg$out_gene_modules_csv)) cfg$out_gene_modules_csv else file.path(out_dir, "gene_modules.csv")

message("[fig3] Loading Seurat object: ", in_rds)
obj <- readRDS(in_rds)
if (is.null(obj$seurat_clusters)) stop("Seurat object missing meta.data column: seurat_clusters", call. = FALSE)
if (is.null(obj@reductions$umap)) stop("Seurat object missing UMAP reduction: obj@reductions$umap", call. = FALSE)

message("[fig3] Converting to Monocle3 cell_data_set ...")
cds <- as.cell_data_set(obj)

cds <- preprocess_cds(cds, num_dim = num_dim)
cds <- reduce_dimension(cds)

# Set gene_short_name (safer than fData in newer monocle3)
try({
  rowData(cds)$gene_short_name <- rownames(rowData(cds))
}, silent = TRUE)

# Force single partition
recreate.partition <- factor(rep(1, length(colnames(cds))))
names(recreate.partition) <- colnames(cds)
cds@clusters$UMAP$partitions <- recreate.partition

# Use Seurat clusters as monocle clusters
cds@clusters$UMAP$clusters <- obj$seurat_clusters

# Overwrite Monocle UMAP with Seurat UMAP
umap_emb <- obj@reductions$umap@cell.embeddings
common_cells <- intersect(colnames(cds), rownames(umap_emb))
cds <- cds[, common_cells]
umap_emb <- umap_emb[common_cells, , drop = FALSE]
cds@int_colData@listData$reducedDims$UMAP <- umap_emb

message("[fig3] Learning principal graph ...")
cds <- learn_graph(cds)

root_cells <- colnames(cds)[as.character(colData(cds)$seurat_clusters) == root_cluster]
if (length(root_cells) == 0) {
  stop("No root cells found for root_cluster=", root_cluster,
       ". Check seurat_clusters values.", call. = FALSE)
}

message("[fig3] Ordering cells (root_cluster=", root_cluster, ", n_root=", length(root_cells), ") ...")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

cds$monocle3_pseudotime <- pseudotime(cds)
obj$pseudotime <- colData(cds)$monocle3_pseudotime

obj$pseudotime_bin_q <- dplyr::ntile(obj$pseudotime, pseudotime_bins)
obj$pseudotime_bin_q <- factor(obj$pseudotime_bin_q, levels = 1:pseudotime_bins, ordered = TRUE)

message("[fig3] graph_test (cores=", cores, ") ...")
deg_oocytes <- graph_test(cds, neighbor_graph = "principal_graph", cores = cores)
write.csv(deg_oocytes, deg_csv, row.names = TRUE)
message("[fig3] Saved DEG table: ", deg_csv)

pr_deg_ids <- rownames(subset(deg_oocytes, q_value < q_value_thr & morans_I > morans_I_thr))
writeLines(pr_deg_ids, con = pr_ids_txt)
message("[fig3] Selected pr_deg_ids: n=", length(pr_deg_ids), " -> ", pr_ids_txt)

if (length(pr_deg_ids) == 0) {
  stop("No genes passed thresholds. Consider relaxing q_value_thr/morans_I_thr.", call. = FALSE)
}

message("[fig3] find_gene_modules (resolution=", module_resolution, ", seed=", module_seed, ") ...")
gene_module_df <- find_gene_modules(cds[pr_deg_ids, ], resolution = module_resolution, random_seed = module_seed)
write.csv(gene_module_df, modules_csv, row.names = FALSE)
message("[fig3] Saved gene modules: ", modules_csv)

saveRDS(obj, out_obj_rds, compress = TRUE)
message("[fig3] Saved Seurat with pseudotime: ", out_obj_rds)

params_used <- list(
  input_rds = in_rds,
  out_dir = out_dir,
  out_obj_rds = out_obj_rds,
  num_dim = num_dim,
  cores = cores,
  root_cluster = root_cluster,
  q_value_thr = q_value_thr,
  morans_I_thr = morans_I_thr,
  pseudotime_bins = pseudotime_bins,
  module_resolution = module_resolution,
  module_seed = module_seed,
  out_deg_csv = deg_csv,
  out_pr_deg_ids = pr_ids_txt,
  out_gene_modules_csv = modules_csv
)

yaml_path <- file.path(out_dir, "params_used_fig3_monocle3.yaml")
yaml::write_yaml(params_used, yaml_path)
message("[fig3] Saved params: ", yaml_path)

si_path <- file.path(out_dir, "sessionInfo_fig3.txt")
write_text(si_path, capture.output(sessionInfo()))
message("[fig3] Saved sessionInfo: ", si_path)

message("[fig3] Done.")
