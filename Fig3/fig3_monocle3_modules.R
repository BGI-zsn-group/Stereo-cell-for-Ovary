#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Fig3 | Monocle3 pseudotime + graph_test DEGs + gene modules
#
# What this script does
#   - Reads an integrated Seurat object (RDS; e.g. Fig2 output).
#   - Converts to Monocle3 cell_data_set, learns trajectory graph, and orders cells
#     in pseudotime using a specified root cluster.
#   - Runs Monocle3 graph_test to identify genes varying along the trajectory.
#   - Finds gene modules from significant DEGs.
#   - Optionally applies a stage-aware prefilter step before the final module call:
#       1) run an initial find_gene_modules()
#       2) find positive stage markers from early/late stages via FindAllMarkers()
#       3) for a target module (default module 3), remove genes that are neither
#          early nor late markers OR are simultaneously both early and late
#       4) rerun find_gene_modules() on the filtered DEG list
#   - Writes the updated Seurat object (with `pseudotime` and binned pseudotime)
#     plus CSV/TXT artifacts used by downstream Fig3 panels.
#
# Recommended way to run (repo wrapper)
#   bash Fig3/run_fig3_combined.sh --only monocle3 -i <obj_oo.rds> -o <out_dir>
#
# Run this script directly (module-level YAML)
#   Rscript Fig3/fig3_monocle3_modules.R --config <fig3_module.yaml>
#
# Key config fields (module-level YAML)
#   input_rds, out_dir, out_obj_rds, root_cluster, num_dim, cores,
#   q_value_thr, morans_I_thr, pseudotime_bins,
#   module_resolution, prefilter_module_resolution, module_seed,
#   stage_col, early_stage_levels, late_stage_levels,
#   prefilter_enabled, prefilter_target_module,
#   out_deg_csv, out_pr_deg_ids, out_gene_modules_csv
#
# Outputs
#   - out_obj_rds: Seurat object with pseudotime (+ bins)
#   - out_deg_csv: DEG table from Monocle3 graph_test
#   - out_pr_deg_ids: final gene IDs used for the final module call
#   - out_gene_modules_csv: final gene-module assignments
#   - optional: out_prefilter_modules_csv (initial module call before filtering)
#   - params_used_fig3_monocle3.yaml, sessionInfo_fig3.txt
#
# Dependencies
#   Seurat, SeuratWrappers, monocle3, dplyr/tidyverse, ggplot2, patchwork, yaml
# -----------------------------------------------------------------------------
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
  if (length(x) == 1 && is.character(x) && grepl(",", x, fixed = TRUE)) {
    v <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  } else {
    v <- unlist(x)
    v <- trimws(as.character(v))
  }
  v[nzchar(v)]
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

pick_lfc_col <- function(df) {
  cand <- c("avg_log2FC", "avg_logFC")
  hit <- cand[cand %in% colnames(df)]
  if (length(hit) == 0) {
    stop("FindAllMarkers() output missing expected logFC column (avg_log2FC/avg_logFC).", call. = FALSE)
  }
  hit[[1]]
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
prefilter_module_resolution <- if (!is.null(cfg$prefilter_module_resolution)) as.numeric(cfg$prefilter_module_resolution) else 0.001
module_seed <- if (!is.null(cfg$module_seed)) as.integer(cfg$module_seed) else 888

prefilter_enabled <- if (!is.null(cfg$prefilter_enabled)) isTRUE(as.logical(cfg$prefilter_enabled)) else TRUE
prefilter_target_module <- if (!is.null(cfg$prefilter_target_module)) as.character(cfg$prefilter_target_module) else "3"
stage_col <- if (!is.null(cfg$stage_col)) as.character(cfg$stage_col) else "stage"
early_stage_levels <- if (!is.null(cfg$early_stage_levels)) as_vars(cfg$early_stage_levels) else c("EGO", "GO1")
late_stage_levels <- if (!is.null(cfg$late_stage_levels)) as_vars(cfg$late_stage_levels) else c("GO3", "FGO")


deg_csv <- if (!is.null(cfg$out_deg_csv)) cfg$out_deg_csv else file.path(out_dir, "deg_graph_test.csv")
pr_ids_txt <- if (!is.null(cfg$out_pr_deg_ids)) cfg$out_pr_deg_ids else file.path(out_dir, "pr_deg_ids.txt")
modules_csv <- if (!is.null(cfg$out_gene_modules_csv)) cfg$out_gene_modules_csv else file.path(out_dir, "gene_modules.csv")
prefilter_modules_csv <- if (!is.null(cfg$out_prefilter_modules_csv)) cfg$out_prefilter_modules_csv else file.path(out_dir, "gene_modules_prefilter.csv")

message("[fig3] Loading Seurat object: ", in_rds)
obj <- readRDS(in_rds)
if (is.null(obj$seurat_clusters)) stop("Seurat object missing meta.data column: seurat_clusters", call. = FALSE)
if (is.null(obj@reductions$umap)) stop("Seurat object missing UMAP reduction: obj@reductions$umap", call. = FALSE)

message("[fig3] Converting to Monocle3 cell_data_set ...")
cds <- as.cell_data_set(obj)
cds <- preprocess_cds(cds, num_dim = num_dim)
cds <- reduce_dimension(cds)

# Set gene_short_name robustly
if ("gene_short_name" %in% colnames(rowData(cds))) {
  rowData(cds)$gene_short_name <- rownames(rowData(cds))
} else {
  try({
    rowData(cds)$gene_short_name <- rownames(rowData(cds))
  }, silent = TRUE)
  try({
    fData(cds)$gene_short_name <- rownames(fData(cds))
  }, silent = TRUE)
}

# Force single partition and reuse Seurat clusters/UMAP
recreate.partition <- factor(rep(1, length(colnames(cds))))
names(recreate.partition) <- colnames(cds)
cds@clusters$UMAP$partitions <- recreate.partition
cds@clusters$UMAP$clusters <- obj$seurat_clusters

umap_emb <- obj@reductions$umap@cell.embeddings
common_cells <- intersect(colnames(cds), rownames(umap_emb))
cds <- cds[, common_cells]
umap_emb <- umap_emb[common_cells, , drop = FALSE]
cds@int_colData@listData$reducedDims$UMAP <- umap_emb

message("[fig3] Learning principal graph ...")
cds <- learn_graph(cds)

root_cells <- colnames(cds)[as.character(colData(cds)$seurat_clusters) %in% root_cluster]
if (length(root_cells) == 0) {
  stop("No root cells found for root_cluster=", paste(root_cluster, collapse = ","),
       ". Check seurat_clusters values.", call. = FALSE)
}
message("[fig3] Ordering cells (root_cluster=", paste(root_cluster, collapse = ","), ", n_root=", length(root_cells), ") ...")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

cds$monocle3_pseudotime <- pseudotime(cds)
obj$pseudotime <- colData(cds)$monocle3_pseudotime
obj$pseudotime_bin_q <- dplyr::ntile(obj$pseudotime, pseudotime_bins)
obj$pseudotime_bin_q <- factor(obj$pseudotime_bin_q, levels = 1:pseudotime_bins, ordered = TRUE)

message("[fig3] graph_test (cores=", cores, ") ...")
deg_oocytes <- graph_test(cds, neighbor_graph = "principal_graph", cores = cores)
write.csv(deg_oocytes, deg_csv, row.names = TRUE)
message("[fig3] Saved DEG table: ", deg_csv)

pr_deg_ids_initial <- rownames(subset(deg_oocytes, q_value < q_value_thr & morans_I > morans_I_thr))
if (length(pr_deg_ids_initial) == 0) {
  stop("No genes passed thresholds. Consider relaxing q_value_thr/morans_I_thr.", call. = FALSE)
}
message("[fig3] Initial pr_deg_ids: n=", length(pr_deg_ids_initial))

pr_deg_ids_final <- pr_deg_ids_initial
prefilter_summary <- list(enabled = prefilter_enabled)

if (prefilter_enabled) {
  if (!stage_col %in% colnames(obj@meta.data)) {
    stop("prefilter_enabled=TRUE but stage_col not found in meta.data: ", stage_col, call. = FALSE)
  }

  message("[fig3] Initial find_gene_modules for prefilter (resolution=", prefilter_module_resolution,
          ", seed=", module_seed, ") ...")
  gene_module_df_prefilter <- find_gene_modules(
    cds[pr_deg_ids_initial, ],
    resolution = prefilter_module_resolution,
    random_seed = module_seed
  )
  write.csv(gene_module_df_prefilter, prefilter_modules_csv, row.names = FALSE)
  message("[fig3] Saved prefilter modules: ", prefilter_modules_csv)

  old_ident <- Idents(obj)
  on.exit({
    try(Idents(obj) <- old_ident, silent = TRUE)
  }, add = TRUE)
  Idents(obj) <- stage_col
  marker_df <- FindAllMarkers(obj)
  lfc_col <- pick_lfc_col(marker_df)

  early_deg <- marker_df %>%
    filter(.data[[lfc_col]] > 0, cluster %in% early_stage_levels) %>%
    pull(gene) %>%
    unique()

  late_deg <- marker_df %>%
    filter(.data[[lfc_col]] > 0, cluster %in% late_stage_levels) %>%
    pull(gene) %>%
    unique()

  target_module_genes <- gene_module_df_prefilter %>%
    filter(as.character(module) == prefilter_target_module) %>%
    pull(id) %>%
    unique()

  target_early <- intersect(target_module_genes, early_deg)
  target_late <- intersect(target_module_genes, late_deg)
  target_both <- intersect(target_early, target_late)

  target_keep <- union(
    setdiff(target_module_genes, union(target_early, target_late)),
    target_both
  )

  pr_deg_ids_final <- setdiff(pr_deg_ids_initial, target_keep)
  if (length(pr_deg_ids_final) == 0) {
    stop("All genes were filtered out by stage-aware module prefilter. Disable prefilter or adjust stage/module settings.", call. = FALSE)
  }

  prefilter_summary <- list(
    enabled = TRUE,
    stage_col = stage_col,
    early_stage_levels = early_stage_levels,
    late_stage_levels = late_stage_levels,
    prefilter_target_module = prefilter_target_module,
    prefilter_module_resolution = prefilter_module_resolution,
    initial_pr_deg_n = length(pr_deg_ids_initial),
    target_module_gene_n = length(target_module_genes),
    early_marker_n = length(early_deg),
    late_marker_n = length(late_deg),
    removed_gene_n = length(target_keep),
    final_pr_deg_n = length(pr_deg_ids_final)
  )

  message("[fig3] Stage-aware prefilter applied: target module=", prefilter_target_module,
          ", removed=", length(target_keep),
          ", final pr_deg_ids=", length(pr_deg_ids_final))
}

writeLines(pr_deg_ids_final, con = pr_ids_txt)
message("[fig3] Saved final pr_deg_ids: n=", length(pr_deg_ids_final), " -> ", pr_ids_txt)

message("[fig3] Final find_gene_modules (resolution=", module_resolution, ", seed=", module_seed, ") ...")
gene_module_df <- find_gene_modules(cds[pr_deg_ids_final, ], resolution = module_resolution, random_seed = module_seed)
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
  prefilter_module_resolution = prefilter_module_resolution,
  module_seed = module_seed,
  prefilter_enabled = prefilter_enabled,
  prefilter_target_module = prefilter_target_module,
  stage_col = stage_col,
  early_stage_levels = early_stage_levels,
  late_stage_levels = late_stage_levels,
  out_deg_csv = deg_csv,
  out_pr_deg_ids = pr_ids_txt,
  out_gene_modules_csv = modules_csv,
  out_prefilter_modules_csv = prefilter_modules_csv,
  prefilter_summary = prefilter_summary
)

yaml_path <- file.path(out_dir, "params_used_fig3_monocle3.yaml")
yaml::write_yaml(params_used, yaml_path)
message("[fig3] Saved params: ", yaml_path)

si_path <- file.path(out_dir, "sessionInfo_fig3.txt")
write_text(si_path, capture.output(sessionInfo()))
message("[fig3] Saved sessionInfo: ", si_path)

message("[fig3] Done.")
