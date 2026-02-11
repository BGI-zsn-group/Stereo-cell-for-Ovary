#!/usr/bin/env Rscript
# Fig4: Cell2location single-cell reference preparation (Seurat -> AnnData)
#
# Purpose
# - Replace the GC population in a chosen sample (default: P49) with a refined GC subset
#   (e.g., obj_gr_newanoo_15rs.rds), then merge with the remaining somatic cells from
#   the full annotated object (obj_all_withanno.rds), and export to .h5ad for Cell2location.
#
# Usage
#   Rscript fig4_cell2location_sc_prep.R --config Fig4/configs/fig4_cell2location_sc_prep.yaml
#
suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(sceasy)
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
"fig4_cell2location_sc_prep.R

Required:
  --config <yaml>     YAML config path

Example:
  Rscript fig4_cell2location_sc_prep.R --config Fig4/configs/fig4_cell2location_sc_prep.yaml
", sep = "")
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

apply_celltype_map <- function(obj, cluster_to_celltype, levels_vec = NULL, overwrite_col = "celltype") {
  # cluster_to_celltype: named list mapping celltype -> vector of clusters (strings)
  if (is.null(cluster_to_celltype) || length(cluster_to_celltype) == 0) return(obj)
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("Missing `seurat_clusters` in meta.data.", call. = FALSE)

  obj[[overwrite_col]] <- NA_character_
  for (ct in names(cluster_to_celltype)) {
    cls <- as.character(unlist(cluster_to_celltype[[ct]]))
    obj[[overwrite_col]][obj$seurat_clusters %in% cls] <- ct
  }
  if (!is.null(levels_vec) && length(levels_vec) > 0) {
    obj[[overwrite_col]] <- factor(obj[[overwrite_col]], levels = as.character(levels_vec))
  } else {
    obj[[overwrite_col]] <- factor(obj[[overwrite_col]])
  }
  obj
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

cfg <- yaml::read_yaml(args$config)

# I/O
obj_gr_rds <- as.character(cfg$obj_gr_rds %||% "obj_gr_newanoo_15rs.rds")
obj_total_rds <- as.character(cfg$obj_total_rds %||% "obj_all_withanno.rds")

if (!file.exists(obj_gr_rds)) stop("Missing obj_gr_rds: ", obj_gr_rds, call. = FALSE)
if (!file.exists(obj_total_rds)) stop("Missing obj_total_rds: ", obj_total_rds, call. = FALSE)

out_dir <- as.character(cfg$out_dir %||% "results/Fig4/cell2location_sc_ref")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- as.character(cfg$prefix %||% "somatic")
out_h5ad <- as.character(cfg$out_h5ad %||% file.path(out_dir, paste0(prefix, ".h5ad")))
out_rds <- as.character(cfg$out_merged_rds %||% file.path(out_dir, paste0(prefix, ".merged_for_cell2location.rds")))
main_layer <- as.character(cfg$main_layer %||% "counts")

target_sample <- as.character(cfg$target_sample %||% "P49")
sample_col <- as.character(cfg$sample_col %||% "sample")
celltype_col_total <- as.character(cfg$celltype_col_total %||% "celltype")

log_path <- file.path(out_dir, paste0(prefix, ".cell2location_sc_prep.log"))
log_file <- file(log_path, open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

message("[fig4-c2l-sc] start: ", Sys.time())
message("[fig4-c2l-sc] target_sample: ", target_sample)
message("[fig4-c2l-sc] obj_gr_rds: ", obj_gr_rds)
message("[fig4-c2l-sc] obj_total_rds: ", obj_total_rds)
message("[fig4-c2l-sc] out_h5ad: ", out_h5ad)
message("[fig4-c2l-sc] main_layer: ", main_layer)

# Load GC refined object (obj_gr_*)
obj_gr <- readRDS(obj_gr_rds)
if (!inherits(obj_gr, "Seurat")) stop("obj_gr_rds must be a Seurat object.", call. = FALSE)
if (!sample_col %in% colnames(obj_gr@meta.data)) stop("obj_gr missing sample column: ", sample_col, call. = FALSE)

# Optionally overwrite celltype in obj_gr using cluster mapping (recommended for reproducibility)
if (!is.null(cfg$gc_celltype_levels) && !is.null(cfg$gc_cluster_to_celltype)) {
  message("[fig4-c2l-sc] applying GC cluster->celltype mapping on obj_gr ...")
  obj_gr <- apply_celltype_map(
    obj_gr,
    cluster_to_celltype = cfg$gc_cluster_to_celltype,
    levels_vec = cfg$gc_celltype_levels,
    overwrite_col = "celltype"
  )
}

# Keep only target sample's GC
obj_gr_target <- subset(obj_gr, subset = get(sample_col) == target_sample)
message("[fig4-c2l-sc] obj_gr_target cells: ", ncol(obj_gr_target))

# Load full annotated object and remove GC-like labels for the same target sample
obj_total <- readRDS(obj_total_rds)
if (!inherits(obj_total, "Seurat")) stop("obj_total_rds must be a Seurat object.", call. = FALSE)
if (!sample_col %in% colnames(obj_total@meta.data)) stop("obj_total missing sample column: ", sample_col, call. = FALSE)
if (!celltype_col_total %in% colnames(obj_total@meta.data)) stop("obj_total missing celltype column: ", celltype_col_total, call. = FALSE)

exclude_celltypes <- as.character(unlist(cfg$exclude_celltypes_in_total %||% c("GC_Atretic","GC_Preantral","GC_Antral_Mural","GC_Column")))
message("[fig4-c2l-sc] exclude celltypes in obj_total (within target sample): ", paste(exclude_celltypes, collapse = ", "))

obj_me <- subset(
  obj_total,
  subset = (get(sample_col) == target_sample) & !(get(celltype_col_total) %in% exclude_celltypes)
)
message("[fig4-c2l-sc] obj_me cells: ", ncol(obj_me))

# Merge: refined GC (target sample) + remaining somatic cells (target sample)
obj_total_new <- merge(obj_gr_target, y = obj_me)
message("[fig4-c2l-sc] merged cells: ", ncol(obj_total_new))

# Save merged Seurat (optional but recommended)
saveRDS(obj_total_new, out_rds)
message("[fig4-c2l-sc] saved merged Seurat: ", out_rds)

# Export to AnnData (.h5ad)
sceasy::convertFormat(
  obj_total_new,
  from = "seurat",
  to = "anndata",
  outFile = out_h5ad,
  main_layer = main_layer
)
message("[fig4-c2l-sc] saved h5ad: ", out_h5ad)

# Provenance
params_used <- cfg
params_used$resolved <- list(
  obj_gr_rds = obj_gr_rds,
  obj_total_rds = obj_total_rds,
  out_dir = out_dir,
  out_h5ad = out_h5ad,
  out_merged_rds = out_rds,
  log_path = log_path
)
yaml::write_yaml(params_used, file.path(out_dir, "params_used_fig4_cell2location_sc_prep.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig4_cell2location_sc_prep.txt"), capture.output(sessionInfo()))

message("[fig4-c2l-sc] done: ", Sys.time())

sink()
sink(type = "message")
close(log_file)
