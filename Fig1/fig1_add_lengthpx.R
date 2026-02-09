#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Fig.1 utility: Add per-cell Length_px metadata to Seurat objects and filter missing cells
#
# This script:
#   1) Loads <sample>.rds (Seurat object)
#   2) Loads line length table (default: line_lengths_<sample>.csv) containing columns: Name, Length_px
#   3) Writes Length_px into obj@meta.data by matching cell names
#   4) Filters cells with missing Length_px
#   5) Saves updated Seurat object (overwrite or new file)
#
# Usage:
#   Rscript add_lengthpx.R --samples A1,A2,A3 --overwrite
#
# Options:
#   --samples        Comma-separated sample prefixes (required)
#   --rds-dir        Directory containing <sample>.rds (default: .)
#   --csv-dir        Directory containing CSVs (default: .)
#   --csv-prefix     CSV prefix before sample (default: line_lengths_)
#   --csv-suffix     CSV suffix after sample (default: .csv)
#   --overwrite      If set, overwrite <sample>.rds (default: FALSE)
#   --out-suffix     Suffix for output rds if not overwriting (default: _withLengthpx)
#
# Example:
#   Rscript add_lengthpx.R \
#     --samples A04086F2,A04087F4 \
#     --rds-dir . --csv-dir . \
#     --overwrite
#
suppressPackageStartupMessages({
  library(Seurat)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) return(list())
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (startsWith(key, "--")) {
      k <- sub("^--", "", key)
      if (k %in% c("overwrite", "help")) {
        out[[k]] <- TRUE
        i <- i + 1
        next
      }
      if (i == length(args)) stop("Missing value for ", key)
      out[[k]] <- args[[i + 1]]
      i <- i + 2
    } else {
      stop("Unknown argument: ", key)
    }
  }
  out
}

print_help <- function() {
  cat(
"add_lengthpx.R

Required:
  --samples        Comma-separated sample prefixes (e.g., A1,A2,A3)

Optional:
  --rds-dir        Directory containing <sample>.rds (default: .)
  --csv-dir        Directory containing CSVs (default: .)
  --csv-prefix     CSV prefix (default: line_lengths_)
  --csv-suffix     CSV suffix (default: .csv)
  --overwrite      Overwrite <sample>.rds (default: FALSE)
  --out-suffix     Output suffix if not overwriting (default: _withLengthpx)

Example:
  Rscript add_lengthpx.R --samples A04086F2,A04087F4 --overwrite
", sep = "")
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) {
  print_help()
  quit(status = 0)
}

if (is.null(args$samples) || nchar(args$samples) == 0) {
  print_help()
  stop("Missing required --samples")
}

samples <- strsplit(args$samples, ",", fixed = TRUE)[[1]]
samples <- trimws(samples)
samples <- samples[nzchar(samples)]
if (length(samples) == 0) stop("No valid samples parsed from --samples")

rds_dir <- if (!is.null(args$`rds-dir`)) args$`rds-dir` else "."
csv_dir <- if (!is.null(args$`csv-dir`)) args$`csv-dir` else "."
csv_prefix <- if (!is.null(args$`csv-prefix`)) args$`csv-prefix` else "line_lengths_"
csv_suffix <- if (!is.null(args$`csv-suffix`)) args$`csv-suffix` else ".csv"
overwrite <- !is.null(args$overwrite) && isTRUE(args$overwrite)
out_suffix <- if (!is.null(args$`out-suffix`)) args$`out-suffix` else "_withLengthpx"

for (s in samples) {
  rds_file <- file.path(rds_dir, paste0(s, ".rds"))
  csv_file <- file.path(csv_dir, paste0(csv_prefix, s, csv_suffix))

  message("=== Processing: ", s, " ===")
  if (!file.exists(rds_file)) stop("Missing rds: ", rds_file)
  if (!file.exists(csv_file)) stop("Missing csv: ", csv_file)

  obj <- readRDS(rds_file)
  length_df <- read.csv(csv_file, stringsAsFactors = FALSE)

  req_cols <- c("Name", "Length_px")
  miss <- setdiff(req_cols, colnames(length_df))
  if (length(miss) > 0) stop("CSV missing required columns: ", paste(miss, collapse = ", "), " | file: ", csv_file)

  idx <- match(rownames(obj@meta.data), length_df$Name)
  obj@meta.data$Length_px <- length_df$Length_px[idx]

  n_before <- ncol(obj)
  n_has <- sum(!is.na(obj@meta.data$Length_px))
  n_missing <- n_before - n_has

  cells_keep <- colnames(obj)[!is.na(obj@meta.data$Length_px)]
  obj <- subset(obj, cells = cells_keep)
  n_after <- ncol(obj)

  out_file <- if (overwrite) rds_file else file.path(rds_dir, paste0(s, out_suffix, ".rds"))
  saveRDS(obj, out_file, compress = TRUE)

  message("Saved: ", out_file,
          " | cells: ", n_before, " -> ", n_after,
          " (missing removed = ", n_missing, ")")
}

message("All done.")
