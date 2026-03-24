#!/usr/bin/env Rscript
# Fig4: Cell2location single-cell reference preparation (patched v3)
#
# Purpose
# - Replace the GC population in a chosen sample (default: P49) with a refined GC subset
#   and merge it with the remaining non-GC somatic cells from the full annotated object.
# - Export the merged Seurat object to .h5ad for Cell2location.
#
# Key patches in v3
# - robust metadata subsetting (avoid subset/get pitfalls)
# - robust cluster->celltype mapping in obj@meta.data
# - explicit reticulate/python binding before sceasy::convertFormat
# - never try to auto-download uv / create ephemeral python envs

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
    if (!startsWith(k, "--")) stop("Unknown arg: ", k, call. = FALSE)
    key <- sub("^--", "", k)
    if (key %in% c("help")) {
      out[[key]] <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) stop("Missing value for ", k, call. = FALSE)
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

Optional YAML keys:
  python_bin: /path/to/python

Env fallback:
  FIG4_C2L_SC_PYTHON or RETICULATE_PYTHON

Example:
  Rscript fig4_cell2location_sc_prep.R --config Fig4/configs/fig4_cell2location_sc_prep.yaml
", sep = "")
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

normalize_cluster_map <- function(cluster_to_celltype) {
  if (is.null(cluster_to_celltype) || length(cluster_to_celltype) == 0) return(list())
  out <- list()
  for (nm in names(cluster_to_celltype)) {
    out[[nm]] <- as.character(unlist(cluster_to_celltype[[nm]], use.names = FALSE))
  }
  out
}

ensure_fresh_meta_col <- function(obj, colname) {
  md <- obj@meta.data
  if (sum(colnames(md) == colname) > 0) {
    md <- md[, colnames(md) != colname, drop = FALSE]
  }
  md[[colname]] <- NA_character_
  obj@meta.data <- md
  obj
}

apply_celltype_map <- function(obj, cluster_to_celltype, levels_vec = NULL, overwrite_col = "celltype") {
  if (is.null(cluster_to_celltype) || length(cluster_to_celltype) == 0) return(obj)
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    stop("Missing `seurat_clusters` in meta.data.", call. = FALSE)
  }

  cluster_to_celltype <- normalize_cluster_map(cluster_to_celltype)
  obj <- ensure_fresh_meta_col(obj, overwrite_col)
  md <- obj@meta.data
  clusters_chr <- as.character(md$seurat_clusters)

  for (ct in names(cluster_to_celltype)) {
    cls <- as.character(cluster_to_celltype[[ct]])
    idx <- which(clusters_chr %in% cls)
    if (length(idx) > 0) md[idx, overwrite_col] <- ct
  }

  if (!is.null(levels_vec) && length(levels_vec) > 0) {
    md[[overwrite_col]] <- factor(md[[overwrite_col]], levels = as.character(levels_vec))
  } else {
    md[[overwrite_col]] <- factor(md[[overwrite_col]])
  }

  obj@meta.data <- md
  obj
}

subset_by_meta_equal <- function(obj, colname, value) {
  md <- obj@meta.data
  if (!colname %in% colnames(md)) stop("Missing metadata column: ", colname, call. = FALSE)
  keep_cells <- rownames(md)[as.character(md[[colname]]) == as.character(value)]
  if (length(keep_cells) == 0) {
    stop("No cells matched ", colname, " == ", value, call. = FALSE)
  }
  subset(obj, cells = keep_cells)
}

subset_total_non_gc <- function(obj, sample_col, target_sample, celltype_col, exclude_celltypes) {
  md <- obj@meta.data
  miss <- setdiff(c(sample_col, celltype_col), colnames(md))
  if (length(miss) > 0) stop("Missing metadata columns in obj_total: ", paste(miss, collapse = ", "), call. = FALSE)

  keep <- as.character(md[[sample_col]]) == as.character(target_sample) &
    !(as.character(md[[celltype_col]]) %in% as.character(exclude_celltypes))

  keep_cells <- rownames(md)[keep]
  if (length(keep_cells) == 0) {
    stop("No non-GC cells kept from obj_total for target sample ", target_sample, call. = FALSE)
  }
  subset(obj, cells = keep_cells)
}

resolve_python_bin <- function(cfg) {
  candidates <- c(
    as.character(cfg$python_bin %||% ""),
    Sys.getenv("FIG4_C2L_SC_PYTHON", unset = ""),
    Sys.getenv("RETICULATE_PYTHON", unset = ""),
    Sys.which("python3"),
    Sys.which("python")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }
  ""
}

bind_python_for_sceasy <- function(cfg) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required by sceasy but is not installed.", call. = FALSE)
  }

  py <- resolve_python_bin(cfg)
  if (!nzchar(py)) {
    stop(
      paste0(
        "Could not find a usable Python binary for sceasy. ",
        "Set one of: cfg$python_bin, FIG4_C2L_SC_PYTHON, or RETICULATE_PYTHON."
      ),
      call. = FALSE
    )
  }

  # prevent reticulate from trying to auto-provision a fresh env via uv
  Sys.setenv(RETICULATE_PYTHON = py)
  Sys.setenv(RETICULATE_AUTOINSTALL = "FALSE")

  reticulate::use_python(py, required = TRUE)
  cfg_py <- reticulate::py_config()
  message("[fig4-c2l-sc] reticulate python: ", cfg_py$python)
  message("[fig4-c2l-sc] python version: ", cfg_py$version)

  needed <- c("anndata", "numpy", "pandas")
  missing <- needed[!vapply(needed, reticulate::py_module_available, logical(1))]
  if (length(missing) > 0) {
    stop(
      paste0(
        "Selected Python is missing required modules for h5ad export: ",
        paste(missing, collapse = ", "),
        ". Use a Python env that already has these packages."
      ),
      call. = FALSE
    )
  }

  invisible(py)
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
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_file), silent = TRUE)
}, add = TRUE)

message("[fig4-c2l-sc] start: ", Sys.time())
message("[fig4-c2l-sc] target_sample: ", target_sample)
message("[fig4-c2l-sc] obj_gr_rds: ", obj_gr_rds)
message("[fig4-c2l-sc] obj_total_rds: ", obj_total_rds)
message("[fig4-c2l-sc] out_h5ad: ", out_h5ad)
message("[fig4-c2l-sc] main_layer: ", main_layer)

# Load GC refined object
obj_gr <- readRDS(obj_gr_rds)
if (!inherits(obj_gr, "Seurat")) stop("obj_gr_rds must be a Seurat object.", call. = FALSE)
if (!sample_col %in% colnames(obj_gr@meta.data)) stop("obj_gr missing sample column: ", sample_col, call. = FALSE)

# Apply cluster -> celltype mapping on obj_gr when provided
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
obj_gr_target <- subset_by_meta_equal(obj_gr, sample_col, target_sample)
message("[fig4-c2l-sc] obj_gr_target cells: ", ncol(obj_gr_target))

# Load total annotated object
obj_total <- readRDS(obj_total_rds)
if (!inherits(obj_total, "Seurat")) stop("obj_total_rds must be a Seurat object.", call. = FALSE)

exclude_celltypes <- as.character(unlist(cfg$exclude_celltypes_in_total %||% c(
  "GC_Atretic", "GC_Preantral", "GC_Antral_Mural", "GC_Column"
)))
message("[fig4-c2l-sc] exclude celltypes in obj_total (within target sample): ", paste(exclude_celltypes, collapse = ", "))

obj_me <- subset_total_non_gc(
  obj_total,
  sample_col = sample_col,
  target_sample = target_sample,
  celltype_col = celltype_col_total,
  exclude_celltypes = exclude_celltypes
)
message("[fig4-c2l-sc] obj_me cells: ", ncol(obj_me))

# Merge refined GC + remaining non-GC somatic cells
obj_total_new <- merge(obj_gr_target, y = obj_me)
message("[fig4-c2l-sc] merged cells: ", ncol(obj_total_new))

saveRDS(obj_total_new, out_rds)
message("[fig4-c2l-sc] saved merged Seurat: ", out_rds)

# Explicitly bind to an existing python before any anndata conversion.
py_used <- bind_python_for_sceasy(cfg)
message("[fig4-c2l-sc] exporting via sceasy using python: ", py_used)

sceasy::convertFormat(
  obj_total_new,
  from = "seurat",
  to = "anndata",
  outFile = out_h5ad,
  main_layer = main_layer
)
message("[fig4-c2l-sc] saved h5ad: ", out_h5ad)

params_used <- cfg
params_used$resolved <- list(
  obj_gr_rds = obj_gr_rds,
  obj_total_rds = obj_total_rds,
  python_bin = py_used,
  out_dir = out_dir,
  out_h5ad = out_h5ad,
  out_merged_rds = out_rds,
  log_path = log_path
)
yaml::write_yaml(params_used, file.path(out_dir, "params_used_fig4_cell2location_sc_prep.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig4_cell2location_sc_prep.txt"), capture.output(sessionInfo()))

message("[fig4-c2l-sc] done: ", Sys.time())
