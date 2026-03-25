#!/usr/bin/env Rscript
# Fig4: Monocle3 trajectory inference for GC subset (updated)
# Supports rooting either by seurat cluster ID or by a metadata column/value pair.

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

assay <- as.character(cfg$assay %||% "RNA")
num_dim <- as.integer(cfg$num_dim %||% 50)
use_seurat_umap <- as.logical(cfg$use_seurat_umap %||% TRUE)
out_cds_rds <- as.character(cfg$out_cds_rds %||% file.path(out_dir, paste0(prefix, ".cds.rds")))
out_obj_rds <- as.character(cfg$out_obj_rds %||% file.path(out_dir, paste0(prefix, ".obj_with_pseudotime.rds")))
out_pt_csv <- as.character(cfg$out_pseudotime_csv %||% file.path(out_dir, paste0(prefix, ".pseudotime.csv")))
out_plot_pdf <- as.character(cfg$out_plot_pdf %||% file.path(out_dir, paste0(prefix, ".monocle3_pseudotime_umap.pdf")))

obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Expected a Seurat object in input_rds.", call. = FALSE)
if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("Seurat object missing `seurat_clusters` in meta.data.", call. = FALSE)

cds <- as.cell_data_set(obj, assay = assay)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = num_dim)
cds <- reduce_dimension(cds)
try({ fData(cds)$gene_short_name <- rownames(fData(cds)) }, silent = TRUE)

# inject cluster and optionally seurat UMAP
colData(cds)$seurat_clusters <- obj$seurat_clusters
if (use_seurat_umap && "umap" %in% names(obj@reductions)) {
  reducedDims(cds)$UMAP <- obj@reductions$umap@cell.embeddings[colnames(cds), , drop = FALSE]
}

cds@clusters$UMAP$clusters <- factor(colData(cds)$seurat_clusters)
cds@clusters$UMAP$partitions <- factor(rep(1, ncol(cds)))
cds <- learn_graph(cds, use_partition = FALSE)

root_cells <- NULL
root_by_col <- as.character(cfg$root_by_col %||% "")
root_value <- as.character(cfg$root_value %||% "")
if (nzchar(root_by_col) && nzchar(root_value)) {
  if (!root_by_col %in% colnames(obj@meta.data)) stop("root_by_col not found in obj@meta.data: ", root_by_col, call. = FALSE)
  root_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[root_by_col]]) == root_value]
  message("[fig4-monocle3] root by column: ", root_by_col, " == ", root_value, "; n_root_cells=", length(root_cells))
} else {
  root_cluster <- as.character(cfg$root_cluster %||% "8")
  root_cells <- rownames(obj@meta.data)[as.character(obj@meta.data$seurat_clusters) == root_cluster]
  message("[fig4-monocle3] root by cluster: ", root_cluster, "; n_root_cells=", length(root_cells))
}
if (length(root_cells) == 0) stop("No root cells found for monocle3 ordering.", call. = FALSE)

cds <- order_cells(cds, root_cells = intersect(root_cells, colnames(cds)))
pseudotime_vals <- monocle3::pseudotime(cds)
obj$pseudotime <- pseudotime_vals[colnames(obj)]

saveRDS(cds, out_cds_rds)
saveRDS(obj, out_obj_rds)
write.csv(data.frame(cell = names(obj$pseudotime), pseudotime = obj$pseudotime), out_pt_csv, row.names = FALSE)

pdf(out_plot_pdf, width = 7, height = 6)
print(plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE))
dev.off()

yaml::write_yaml(cfg, file.path(out_dir, "params_used_fig4_monocle3.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig4_monocle3.txt"), capture.output(sessionInfo()))
message("[fig4-monocle3] done: ", Sys.time())

sink()
sink(type = "message")
close(log_file)
