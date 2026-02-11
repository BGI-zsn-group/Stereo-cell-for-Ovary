#!/usr/bin/env Rscript
# Fig4: Monocle3 trajectory inference for GC subset (standalone, YAML-driven)
#
# Inputs:
#   - A Seurat object with UMAP computed and `seurat_clusters` in meta.data.
#
# Outputs (under out_dir):
#   - cds_gc.rds
#   - obj_gc_with_pseudotime.rds
#   - pseudotime.csv
#   - monocle3_pseudotime_umap.pdf
#   - params_used_fig4_monocle3.yaml
#   - sessionInfo_fig4_monocle3.txt
#   - fig4_monocle3.log
#
# Usage:
#   Rscript fig4_monocle3.R --config Fig4/configs/fig4_monocle3.yaml
#
suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

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
"fig4_monocle3.R

Required:
  --config <yaml>     YAML config path

Example:
  Rscript fig4_monocle3.R --config Fig4/configs/fig4_monocle3.yaml
", sep = "")
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

cfg <- yaml::read_yaml(args$config)

in_rds <- as.character(cfg$input_rds %||% "results/Fig4/gc_processing/obj_gr_newanno.rds")
if (!file.exists(in_rds)) stop("Missing input_rds: ", in_rds, call. = FALSE)

out_dir <- as.character(cfg$out_dir %||% "results/Fig4/monocle3")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- as.character(cfg$prefix %||% "gc")
log_path <- file.path(out_dir, paste0(prefix, ".fig4_monocle3.log"))

log_file <- file(log_path, open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

message("[fig4-monocle3] start: ", Sys.time())
message("[fig4-monocle3] input_rds: ", in_rds)
message("[fig4-monocle3] out_dir: ", out_dir)

assay <- as.character(cfg$assay %||% "RNA")
num_dim <- as.integer(cfg$num_dim %||% 50)
root_cluster <- as.character(cfg$root_cluster %||% "8")
use_seurat_umap <- as.logical(cfg$use_seurat_umap %||% TRUE)
reduction_method <- as.character(cfg$reduction_method %||% "UMAP")

out_cds_rds <- as.character(cfg$out_cds_rds %||% file.path(out_dir, paste0(prefix, ".cds.rds")))
out_obj_rds <- as.character(cfg$out_obj_rds %||% file.path(out_dir, paste0(prefix, ".obj_with_pseudotime.rds")))
out_pt_csv <- as.character(cfg$out_pseudotime_csv %||% file.path(out_dir, paste0(prefix, ".pseudotime.csv")))
out_plot_pdf <- as.character(cfg$out_plot_pdf %||% file.path(out_dir, paste0(prefix, ".monocle3_pseudotime_umap.pdf")))

obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Expected a Seurat object in input_rds.", call. = FALSE)
if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("Seurat object missing `seurat_clusters` in meta.data.", call. = FALSE)

# Convert to monocle3 cell_data_set
cds <- as.cell_data_set(obj, assay = assay)
cds <- estimate_size_factors(cds)

# Preprocess + reduce dimension
cds <- preprocess_cds(cds, num_dim = num_dim)
cds <- reduce_dimension(cds)

# Ensure gene_short_name exists
try({
  fData(cds)$gene_short_name <- rownames(fData(cds))
}, silent = TRUE)

# Force a single partition (to match your code)
recreate.partition <- rep(1, ncol(cds))
names(recreate.partition) <- colnames(cds)
cds@clusters$UMAP$partitions <- as.factor(recreate.partition)

# Use Seurat clusters as Monocle clusters
cds@clusters$UMAP$clusters <- obj$seurat_clusters

# Optionally inject Seurat UMAP into Monocle reducedDims
if (use_seurat_umap) {
  if (!"umap" %in% names(obj@reductions)) stop("use_seurat_umap=TRUE but Seurat object has no UMAP reduction.", call. = FALSE)
  cds@int_colData@listData$reducedDims$UMAP <- obj@reductions$umap@cell.embeddings
}

# Learn graph + order cells
cds <- learn_graph(cds)

root_cells <- colnames(obj)[obj$seurat_clusters == root_cluster]
if (length(root_cells) == 0) {
  stop("No root cells found for root_cluster=", root_cluster, ". Check your cluster IDs.", call. = FALSE)
}

cds <- order_cells(cds, reduction_method = reduction_method, root_cells = root_cells)

# Save pseudotime back into Seurat object
pt <- pseudotime(cds)
obj$monocle3_pseudotime <- as.numeric(pt[colnames(obj)])
# (optional convenience alias)
obj$pseudotime <- obj$monocle3_pseudotime

# Write outputs
saveRDS(cds, out_cds_rds)
saveRDS(obj, out_obj_rds)

pt_df <- data.frame(Cell = names(pt), pseudotime = as.numeric(pt), stringsAsFactors = FALSE)
write.csv(pt_df, out_pt_csv, row.names = FALSE)

# Plot: monocle3 pseudotime on UMAP
p <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
ggsave(out_plot_pdf, p, width = 6.5, height = 5.5)

# Provenance
params_used <- cfg
params_used$resolved <- list(
  input_rds = in_rds,
  out_dir = out_dir,
  out_cds_rds = out_cds_rds,
  out_obj_rds = out_obj_rds,
  out_pseudotime_csv = out_pt_csv,
  out_plot_pdf = out_plot_pdf,
  log_path = log_path
)
yaml::write_yaml(params_used, file.path(out_dir, "params_used_fig4_monocle3.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig4_monocle3.txt"), capture.output(sessionInfo()))

message("[fig4-monocle3] done: ", Sys.time())

sink()
sink(type = "message")
close(log_file)
