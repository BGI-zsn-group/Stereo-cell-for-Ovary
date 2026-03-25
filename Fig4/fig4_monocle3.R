#!/usr/bin/env Rscript
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
  out <- list(); i <- 1
  while (i <= length(args)) {
    k <- args[[i]]
    if (!startsWith(k, "--")) stop("Unknown arg: ", k, call. = FALSE)
    key <- sub("^--", "", k)
    if (key == "help") { out[[key]] <- TRUE; i <- i + 1; next }
    if (i == length(args)) stop("Missing value for ", k, call. = FALSE)
    out[[key]] <- args[[i+1]]; i <- i + 2
  }
  out
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) quit(status = 0)
if (is.null(args$config)) stop("Missing --config", call. = FALSE)

cfg <- yaml::read_yaml(args$config)
in_rds <- as.character(cfg$input_rds %||% "results/Fig4/gc_processing/obj_gr_newanno.rds")
out_dir <- as.character(cfg$out_dir %||% "results/Fig4/monocle3")
prefix <- as.character(cfg$prefix %||% "gc")
if (!file.exists(in_rds)) stop("Missing input_rds: ", in_rds, call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, paste0(prefix, ".fig4_monocle3.log"))
log_file <- file(log_path, open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log_file), silent = TRUE)
}, add = TRUE)

assay <- as.character(cfg$assay %||% "RNA")
num_dim <- as.integer(cfg$num_dim %||% 50)
use_seurat_umap <- as.logical(cfg$use_seurat_umap %||% TRUE)
reduction_method <- as.character(cfg$reduction_method %||% "UMAP")
root_cluster <- cfg$root_cluster %||% NULL
root_by_col <- as.character(cfg$root_by_col %||% "")
root_value <- as.character(cfg$root_value %||% "")

out_cds_rds <- as.character(cfg$out_cds_rds %||% file.path(out_dir, paste0(prefix, ".cds.rds")))
out_obj_rds <- as.character(cfg$out_obj_rds %||% file.path(out_dir, paste0(prefix, ".obj_with_pseudotime.rds")))
out_pt_csv <- as.character(cfg$out_pseudotime_csv %||% file.path(out_dir, paste0(prefix, ".pseudotime.csv")))
out_plot_pdf <- as.character(cfg$out_plot_pdf %||% file.path(out_dir, paste0(prefix, ".monocle3_pseudotime_umap.pdf")))

message("[fig4-monocle3] start: ", Sys.time())
obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Expected Seurat object.", call. = FALSE)
if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("Missing seurat_clusters in meta.data.", call. = FALSE)

cds <- as.cell_data_set(obj, assay = assay)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = num_dim)
cds <- reduce_dimension(cds)
try({ fData(cds)$gene_short_name <- rownames(fData(cds)) }, silent = TRUE)
recreate.partition <- rep(1, ncol(cds))
names(recreate.partition) <- colnames(cds)
cds@clusters$UMAP$partitions <- as.factor(recreate.partition)
cds@clusters$UMAP$clusters <- obj$seurat_clusters
if (use_seurat_umap) {
  if (!"umap" %in% names(obj@reductions)) stop("Seurat object has no UMAP reduction.", call. = FALSE)
  cds@int_colData@listData$reducedDims$UMAP <- obj@reductions$umap@cell.embeddings
}
cds <- learn_graph(cds)

root_cells <- character(0)
if (!is.null(root_by_col) && nzchar(root_by_col) && root_by_col %in% colnames(obj@meta.data)) {
  root_cells <- colnames(obj)[as.character(obj@meta.data[[root_by_col]]) == root_value]
}
if (length(root_cells) == 0 && !is.null(root_cluster)) {
  root_cells <- colnames(obj)[as.character(obj$seurat_clusters) == as.character(root_cluster)]
}
if (length(root_cells) == 0) stop("No root cells found. Check root_by_col/root_value or root_cluster.", call. = FALSE)

cds <- order_cells(cds, reduction_method = reduction_method, root_cells = root_cells)
pt <- pseudotime(cds)
obj$monocle3_pseudotime <- as.numeric(pt[colnames(obj)])
obj$pseudotime <- obj$monocle3_pseudotime

saveRDS(cds, out_cds_rds)
saveRDS(obj, out_obj_rds)
write.csv(data.frame(Cell = names(pt), pseudotime = as.numeric(pt), stringsAsFactors = FALSE), out_pt_csv, row.names = FALSE)

p <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE,
                label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave(out_plot_pdf, p, width = 6.5, height = 5.5)

yaml::write_yaml(cfg, file.path(out_dir, "params_used_fig4_monocle3.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig4_monocle3.txt"), capture.output(sessionInfo()))
message("[fig4-monocle3] done: ", Sys.time())
