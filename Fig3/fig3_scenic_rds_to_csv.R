#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Fig3 | SCENIC Step 1: Seurat (.rds) -> expression matrix CSV (cells x genes)
#
# What this script does
#   - Reads a Seurat object (RDS) and optionally filters cells by stage.
#   - Exports an expression matrix for pySCENIC:
#       * rows   = cells
#       * cols   = genes
#     The output CSV is used by Fig3/fig3_scenic_pyscenic.py (Step 2).
#
# Recommended way to run (repo wrapper)
#   bash Fig3/run_fig3_combined.sh --only scenic1 -i <seurat.rds> -o <out_dir>
#
# Run this script directly (module-level YAML)
#   Rscript Fig3/fig3_scenic_rds_to_csv.R --config <fig3_module.yaml>
#
# Key config fields (module-level YAML)
#   input_rds / out_obj_rds:
#       Seurat object used as SCENIC input (default picks out_obj_rds -> input_rds)
#   scenic_stage_column:
#       meta.data column used for stage filtering (default: "stage")
#   scenic_exclude_stage_value:
#       excluded stage value (default: "MII")
#   scenic_assay / scenic_slot:
#       assay/slot to export (default: RNA / counts)
#   scenic_out_csv:
#       output expression CSV (default: results/Fig3/scenic/oocyte_1211_withoutMII.csv)
#   scenic_transpose:
#       if TRUE (default), exports cells x genes; if FALSE, exports genes x cells.
#
# Outputs
#   - scenic_out_csv
#   - params_used_fig3_scenic1.yaml
#   - sessionInfo_fig3_scenic1.txt
#
# Dependencies
#   Seurat, Matrix, dplyr, yaml
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
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
"fig3_scenic_rds_to_csv.R

Required:
  --config <yaml>     YAML config path (reuse Fig3 config module)

Example:
  Rscript Fig3/fig3_scenic_rds_to_csv.R --config Fig3/configs/fig3_module.yaml
", sep = "")
}

need_yaml <- function() {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Please run: install.packages('yaml')", call. = FALSE)
  }
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

# ---- Inputs ----
in_rds <- cfg$scenic_input_rds
if (is.null(in_rds) || !nzchar(as.character(in_rds))) {
  if (!is.null(cfg$out_obj_rds) && nzchar(as.character(cfg$out_obj_rds))) {
    in_rds <- cfg$out_obj_rds
  } else if (!is.null(cfg$input_rds) && nzchar(as.character(cfg$input_rds))) {
    in_rds <- cfg$input_rds
  } else {
    in_rds <- "obj_oo.rds"
  }
}
if (!file.exists(in_rds)) stop("Missing Seurat input (scenic_input_rds/out_obj_rds/input_rds): ", in_rds, call. = FALSE)

stage_col <- if (!is.null(cfg$scenic_stage_column)) as.character(cfg$scenic_stage_column) else "stage"
exclude_stage <- if (!is.null(cfg$scenic_exclude_stage_value)) as.character(cfg$scenic_exclude_stage_value) else "MII"

assay <- if (!is.null(cfg$scenic_assay)) as.character(cfg$scenic_assay) else "RNA"
slot <- if (!is.null(cfg$scenic_slot)) as.character(cfg$scenic_slot) else "counts"
transpose <- if (!is.null(cfg$scenic_transpose)) as.logical(cfg$scenic_transpose) else TRUE

# ---- Outputs ----
scenic_out_dir <- if (!is.null(cfg$scenic_out_dir)) as.character(cfg$scenic_out_dir) else "results/Fig3/scenic"
if (!dir.exists(scenic_out_dir)) dir.create(scenic_out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- if (!is.null(cfg$scenic_out_csv)) as.character(cfg$scenic_out_csv) else file.path(scenic_out_dir, "oocyte_1211_withoutMII.csv")
out_csv_dir <- dirname(out_csv)
if (!dir.exists(out_csv_dir)) dir.create(out_csv_dir, recursive = TRUE, showWarnings = FALSE)

message("[fig3-scenic1] Load Seurat: ", in_rds)
obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Expected a Seurat object in input RDS.", call. = FALSE)

# Filter stage
if (stage_col %in% colnames(obj@meta.data)) {
  message("[fig3-scenic1] Filter stage: ", stage_col, " != ", exclude_stage)
  obj <- subset(obj, subset = get(stage_col) != exclude_stage)
} else {
  message("[fig3-scenic1] stage column not found: ", stage_col, " (skip filtering)")
}

if (!(assay %in% names(obj@assays))) stop("Missing assay in Seurat object: ", assay, call. = FALSE)

mat <- GetAssayData(obj, assay = assay, slot = slot)
if (is.null(mat)) stop("GetAssayData returned NULL for assay=", assay, " slot=", slot, call. = FALSE)

# Export
# Seurat matrices are genes x cells; pySCENIC prefers cells x genes.
out_mtx <- mat
if (transpose) {
  out_mtx <- Matrix::t(out_mtx)
}

message("[fig3-scenic1] Write CSV: ", out_csv, " (transpose=", transpose, ")")
# For reproducibility, always include row names (cell IDs or gene IDs)
write.csv(as.matrix(out_mtx), file = out_csv, quote = FALSE)

# Record params + session
params_used <- list(
  scenic_input_rds = in_rds,
  scenic_stage_column = stage_col,
  scenic_exclude_stage_value = exclude_stage,
  scenic_assay = assay,
  scenic_slot = slot,
  scenic_transpose = transpose,
  scenic_out_dir = scenic_out_dir,
  scenic_out_csv = out_csv,
  n_cells = ncol(obj),
  n_features = nrow(mat)
)
yaml::write_yaml(params_used, file.path(scenic_out_dir, "params_used_fig3_scenic1.yaml"))
write_text(file.path(scenic_out_dir, "sessionInfo_fig3_scenic1.txt"), capture.output(sessionInfo()))

message("[fig3-scenic1] Done.")
