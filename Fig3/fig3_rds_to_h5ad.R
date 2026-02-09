#!/usr/bin/env Rscript
# Fig.3: Convert Seurat .rds -> AnnData .h5ad (for scTour input)
#
# Usage:
#   Rscript fig3_rds_to_h5ad.R --config Fig3/configs/fig3_monocle3.yaml
#
# This step is typically run BEFORE fig3_sctour_tnode.py.
#
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
in_rds <- if (!is.null(cfg$rds2h5ad_input_rds)) cfg$rds2h5ad_input_rds else "oocyte_some_c1_filter_1211.rds"
out_h5ad <- if (!is.null(cfg$rds2h5ad_output_h5ad)) cfg$rds2h5ad_output_h5ad else "oocyte_some_c1_filter_1211.h5ad"
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
