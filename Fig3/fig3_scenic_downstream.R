#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Fig3 | SCENIC downstream analysis (RSS heatmap across pseudotime bins)
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(Seurat)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(); i <- 1
  while (i <= length(args)) {
    k <- args[[i]]
    if (!startsWith(k, "--")) stop("Unknown arg: ", k)
    key <- sub("^--", "", k)
    if (key %in% c("help")) {
      out[[key]] <- TRUE; i <- i + 1; next
    }
    if (i == length(args)) stop("Missing value for ", k)
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

help_msg <- function() {
  cat("fig3_scenic_downstream.R\n\nRequired:\n  --config <yaml>\n", sep = "")
}

write_text <- function(path, txt) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

dir_for <- function(path) {
  d <- dirname(path)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

first_nonempty <- function(...) {
  xs <- list(...)
  for (x in xs) {
    if (!is.null(x) && nzchar(as.character(x))) return(as.character(x))
  }
  NULL
}

canonicalize_ids <- function(x) {
  x <- as.character(x)
  x <- sub("^X", "", x)
  x <- gsub("\\.", "-", x)
  x
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

cfg <- yaml::read_yaml(args$config)

# ---- Inputs ----
loom_path <- first_nonempty(cfg$scenic_result_loom)
if (is.null(loom_path)) {
  scenic_out_dir <- first_nonempty(cfg$scenic_out_dir, file.path("results", "Fig3", "scenic"))
  scenic_prefix <- first_nonempty(cfg$scenic_prefix, "scenic_withoutMII")
  loom_path <- file.path(scenic_out_dir, paste0(scenic_prefix, "_result.loom"))
}
if (!file.exists(loom_path)) stop("Missing scenic_result_loom: ", loom_path, call. = FALSE)

in_rds <- first_nonempty(
  cfg$scenic_downstream_input_rds,
  cfg$scenic_input_rds,
  cfg$input_rds,
  if (!is.null(cfg$out_obj_rds) && nzchar(as.character(cfg$out_obj_rds)) && file.exists(as.character(cfg$out_obj_rds))) cfg$out_obj_rds else NULL
)
if (is.null(in_rds)) stop("Missing scenic_downstream_input_rds/scenic_input_rds/input_rds", call. = FALSE)
if (!file.exists(in_rds)) stop("Missing scenic_downstream input rds: ", in_rds, call. = FALSE)

stage_col <- first_nonempty(cfg$scenic_stage_column, "stage")
exclude_stage <- first_nonempty(cfg$scenic_exclude_stage_value, "MII")
rss_var <- first_nonempty(cfg$scenic_rss_var, "pseudotime_bin_q")
rss_keep_levels <- cfg$scenic_rss_keep_levels
rss_zThreshold <- if (!is.null(cfg$scenic_rss_zThreshold)) as.numeric(cfg$scenic_rss_zThreshold) else 1
rss_thr <- if (!is.null(cfg$scenic_rss_thr)) as.numeric(cfg$scenic_rss_thr) else 0.1
rss_cluster_columns <- if (!is.null(cfg$scenic_rss_cluster_columns)) as.logical(cfg$scenic_rss_cluster_columns) else FALSE
rss_order_rows <- if (!is.null(cfg$scenic_rss_order_rows)) as.logical(cfg$scenic_rss_order_rows) else TRUE
col_low <- first_nonempty(cfg$scenic_rss_col_low, "#005887")
col_mid <- first_nonempty(cfg$scenic_rss_col_mid, "#926987")
col_high <- first_nonempty(cfg$scenic_rss_col_high, "#ED7B86")
revCol <- if (!is.null(cfg$scenic_rss_revCol)) as.logical(cfg$scenic_rss_revCol) else FALSE
verbose <- if (!is.null(cfg$scenic_rss_verbose)) as.logical(cfg$scenic_rss_verbose) else TRUE

# ---- Outputs ----
base_out_dir <- first_nonempty(cfg$scenic_out_dir, file.path("results", "Fig3", "scenic"))
out_dir <- first_nonempty(cfg$scenic_downstream_out_dir, file.path(base_out_dir, "downstream"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rss_csv <- first_nonempty(cfg$scenic_rss_csv, file.path(out_dir, "rss_matrix.csv"))
out_pdf <- first_nonempty(cfg$scenic_rss_plot_pdf, file.path(out_dir, "rss_plot.pdf"))
out_png <- first_nonempty(cfg$scenic_rss_plot_png, file.path(out_dir, "rss_plot.png"))
out_regulons_rds <- first_nonempty(cfg$scenic_regulons_rds, file.path(out_dir, "regulons_list.rds"))
out_thr_rds <- first_nonempty(cfg$scenic_auc_thresholds_rds, file.path(out_dir, "regulon_auc_thresholds.rds"))
for (p in c(out_rss_csv, out_pdf, out_png, out_regulons_rds, out_thr_rds)) dir_for(p)

message("[fig3-scenic-down] Open loom: ", loom_path)
loom <- open_loom(loom_path)
on.exit(try(close_loom(loom), silent = TRUE), add = TRUE)

regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
saveRDS(regulons, out_regulons_rds)

regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
saveRDS(regulonAucThresholds, out_thr_rds)

message("[fig3-scenic-down] Load Seurat: ", in_rds)
obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Expected Seurat object in scenic_downstream input RDS", call. = FALSE)

if (stage_col %in% colnames(obj@meta.data) && nzchar(exclude_stage)) {
  message("[fig3-scenic-down] Filter stage: ", stage_col, " != ", exclude_stage)
  keep_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[stage_col]]) != exclude_stage]
  obj <- subset(obj, cells = keep_cells)
} else {
  message("[fig3-scenic-down] stage filter skipped (column missing or empty exclude value)")
}

if (!rss_var %in% colnames(obj@meta.data)) {
  stop("meta.data missing RSS annotation column: ", rss_var, call. = FALSE)
}

# Exact overlap first, then a light canonicalization fallback.
common_cells <- intersect(colnames(obj), colnames(regulonAUC))
used_canonical <- FALSE
cell_map <- NULL
if (length(common_cells) == 0) {
  obj_can <- canonicalize_ids(colnames(obj))
  auc_can <- canonicalize_ids(colnames(regulonAUC))
  common_can <- intersect(obj_can, auc_can)
  if (length(common_can) > 0) {
    used_canonical <- TRUE
    obj_idx <- match(common_can, obj_can)
    auc_idx <- match(common_can, auc_can)
    cell_map <- data.frame(
      seurat_cell = colnames(obj)[obj_idx],
      loom_cell = colnames(regulonAUC)[auc_idx],
      canonical_cell = common_can,
      stringsAsFactors = FALSE
    )
    common_cells <- cell_map$seurat_cell
  }
}
if (length(common_cells) == 0) {
  stop("No overlapping cells between Seurat and loom regulonAUC.", call. = FALSE)
}

if (used_canonical) {
  message("[fig3-scenic-down] exact overlap=0; canonicalized overlap=", length(common_cells))
  obj <- subset(obj, cells = cell_map$seurat_cell)
  sub_regulonAUC <- regulonAUC[, cell_map$loom_cell, drop = FALSE]
  colnames(sub_regulonAUC) <- cell_map$seurat_cell
} else {
  obj <- subset(obj, cells = common_cells)
  sub_regulonAUC <- regulonAUC[, colnames(obj), drop = FALSE]
}
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ]

cell_ann <- data.frame(row.names = colnames(obj), tmp = obj@meta.data[[rss_var]], stringsAsFactors = FALSE)
colnames(cell_ann) <- rss_var
if (!is.factor(cell_ann[[rss_var]])) {
  vals <- as.character(cell_ann[[rss_var]])
  suppressWarnings(num <- as.numeric(vals))
  if (all(!is.na(num))) {
    levs <- as.character(sort(unique(num)))
    cell_ann[[rss_var]] <- factor(vals, levels = levs, ordered = TRUE)
  } else {
    cell_ann[[rss_var]] <- factor(vals)
  }
}

message("[fig3-scenic-down] calcRSS by: ", rss_var, " (cells=", ncol(sub_regulonAUC), ")")
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), cellAnnotation = cell_ann[colnames(sub_regulonAUC), rss_var, drop = FALSE])
rss <- na.omit(rss)

if (!is.null(rss_keep_levels)) {
  keep <- as.character(unlist(rss_keep_levels))
} else {
  cn <- colnames(rss)
  suppressWarnings(num <- as.numeric(cn))
  if (all(!is.na(num))) {
    keep <- as.character(sort(num))[1:min(10, length(num))]
  } else {
    keep <- cn[1:min(10, length(cn))]
  }
}
keep <- intersect(keep, colnames(rss))
if (length(keep) == 0) stop("No valid columns to keep in RSS. Check scenic_rss_keep_levels.", call. = FALSE)
rss <- rss[, keep, drop = FALSE]

rss_df <- as.data.frame(rss)
rss_df$regulon <- rownames(rss_df)
rss_df <- rss_df[, c("regulon", keep), drop = FALSE]
write.csv(rss_df, out_rss_csv, row.names = FALSE)

message("[fig3-scenic-down] plotRSS -> ", out_pdf)
rssPlot <- plotRSS(
  rss,
  labelsToDiscard = NULL,
  zThreshold = rss_zThreshold,
  cluster_columns = rss_cluster_columns,
  order_rows = rss_order_rows,
  thr = rss_thr,
  varName = rss_var,
  col.low = col_low,
  col.mid = col_mid,
  col.high = col_high,
  revCol = revCol,
  verbose = verbose
)

ggsave(out_pdf, plot = rssPlot, width = 10, height = 8, units = "in")
ggsave(out_png, plot = rssPlot, width = 10, height = 8, units = "in", dpi = 300)

params_used <- list(
  scenic_result_loom = loom_path,
  scenic_downstream_input_rds = in_rds,
  scenic_stage_column = stage_col,
  scenic_exclude_stage_value = exclude_stage,
  scenic_rss_var = rss_var,
  scenic_rss_keep_levels = keep,
  scenic_rss_zThreshold = rss_zThreshold,
  scenic_rss_thr = rss_thr,
  scenic_rss_cluster_columns = rss_cluster_columns,
  scenic_rss_order_rows = rss_order_rows,
  scenic_rss_col_low = col_low,
  scenic_rss_col_mid = col_mid,
  scenic_rss_col_high = col_high,
  scenic_rss_revCol = revCol,
  scenic_rss_verbose = verbose,
  used_canonical_cellname_fallback = used_canonical,
  n_cells_overlap = length(common_cells),
  outputs = list(
    rss_csv = out_rss_csv,
    rss_plot_pdf = out_pdf,
    rss_plot_png = out_png,
    regulons_rds = out_regulons_rds,
    thresholds_rds = out_thr_rds
  )
)
yaml::write_yaml(params_used, file.path(out_dir, "params_used_fig3_scenic_downstream.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig3_scenic_downstream.txt"), capture.output(sessionInfo()))
message("[fig3-scenic-down] Done.")
