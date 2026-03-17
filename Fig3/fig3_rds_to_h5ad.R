#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Fig3 | Convert Seurat .rds -> AnnData .h5ad (for scTour)
#
# What this script does
#   - Loads a Seurat object (RDS) and converts it to AnnData (.h5ad) via sceasy.
#   - This is typically run BEFORE Fig3/fig3_sctour_tnode.py.
#
# Recommended way to run (repo wrapper)
#   bash Fig3/run_fig3_combined.sh --only rds2h5ad -i <seurat.rds> -o <out_dir>
#
# Run this script directly (module-level YAML)
#   Rscript Fig3/fig3_rds_to_h5ad.R --config <fig3_module.yaml>
#
# Key config fields (module-level YAML)
#   rds2h5ad_input_rds, rds2h5ad_output_h5ad, rds2h5ad_main_layer
#
# Outputs
#   - rds2h5ad_output_h5ad
#
# Dependencies
#   Seurat, sceasy, dplyr, yaml
#
# Notes
#   - sceasy relies on a working Python/reticulate environment with `anndata`
#     installed. If conversion fails, check reticulate configuration.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(sceasy)
  library(dplyr)
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
"fig3_rds_to_h5ad.R

Required:
  --config <yaml>     YAML config path (reuse Fig3 config)

Example:
  Rscript fig3_rds_to_h5ad.R --config Fig3/configs/fig3_monocle3.yaml
", sep = "")
}

need_yaml <- function() {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Please run: install.packages('yaml')", call. = FALSE)
  }
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

need_yaml()
cfg <- yaml::read_yaml(args$config)

# Inputs/outputs
in_rds <- NULL
if (!is.null(cfg$rds2h5ad_input_rds) && nzchar(cfg$rds2h5ad_input_rds)) {
  in_rds <- cfg$rds2h5ad_input_rds
} else if (!is.null(cfg$out_obj_rds) && nzchar(cfg$out_obj_rds) && file.exists(cfg$out_obj_rds)) {
  in_rds <- cfg$out_obj_rds
} else if (!is.null(cfg$input_rds) && nzchar(cfg$input_rds)) {
  in_rds <- cfg$input_rds
} else {
  stop("Missing input rds: please provide rds2h5ad_input_rds or input_rds")
}

if (!file.exists(in_rds)) {
  stop("Missing input rds: ", in_rds)
}

out_h5ad <- NULL
if (!is.null(cfg$rds2h5ad_output_h5ad) && nzchar(cfg$rds2h5ad_output_h5ad)) {
  out_h5ad <- cfg$rds2h5ad_output_h5ad
} else if (!is.null(cfg$out_dir) && nzchar(cfg$out_dir)) {
  out_h5ad <- file.path(cfg$out_dir, "sctour", "seurat_to_h5ad.h5ad")
} else {
  out_h5ad <- "seurat_to_h5ad.h5ad"
}

dir.create(dirname(out_h5ad), recursive = TRUE, showWarnings = FALSE)
main_layer <- if (!is.null(cfg$rds2h5ad_main_layer)) cfg$rds2h5ad_main_layer else "data"

if (!file.exists(in_rds)) stop("Missing input rds: ", in_rds, call. = FALSE)

message("[fig3-rds2h5ad] Loading: ", in_rds)
obj <- readRDS(in_rds)

message("[fig3-rds2h5ad] Converting -> ", out_h5ad, " (main_layer=", main_layer, ")")
sceasy::convertFormat(
  obj,
  from = "seurat",
  to = "anndata",
  outFile = out_h5ad,
  main_layer = main_layer
)

message("[fig3-rds2h5ad] Done.")
