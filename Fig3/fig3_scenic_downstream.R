#!/usr/bin/env Rscript
# Fig.3: SCENIC downstream analysis (pySCENIC result loom + Seurat pseudotime bins)
#
# Usage:
#   Rscript fig3_scenic_downstream.R --config Fig3/configs/fig3_monocle3.yaml
#
# Linked I/O behavior:
# - Uses scenic_result_loom (from SCENIC Step2) as the primary input.
# - If scenic_downstream_input_rds is not provided, uses:
#     out_obj_rds (preferred) -> input_rds
#
# Outputs (defaults under results/Fig3/scenic/downstream):
# - rss_matrix.csv
# - rss_plot.pdf / rss_plot.png
# - regulons_list.rds
# - regulon_auc_thresholds.rds
# - params_used_fig3_scenic_downstream.yaml
# - sessionInfo_fig3_scenic_downstream.txt
#
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
"fig3_scenic_downstream.R

Required:
  --config <yaml>     YAML config path (reuse Fig3 config)

Example:
  Rscript fig3_scenic_downstream.R --config Fig3/configs/fig3_monocle3.yaml
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

# ---- Inputs ----
loom_path <- cfg$scenic_result_loom
if (is.null(loom_path) || !nzchar(as.character(loom_path))) {
  loom_path <- "results/Fig3/scenic/oocyte_1211_withoutMII_result.loom"
}
if (!file.exists(loom_path)) stop("Missing scenic_result_loom: ", loom_path, call. = FALSE)

in_rds <- cfg$scenic_downstream_input_rds
if (is.null(in_rds) || !nzchar(as.character(in_rds))) {
  if (!is.null(cfg$out_obj_rds) && nzchar(as.character(cfg$out_obj_rds))) {
    in_rds <- cfg$out_obj_rds
  } else if (!is.null(cfg$input_rds) && nzchar(as.character(cfg$input_rds))) {
    in_rds <- cfg$input_rds
  } else {
    in_rds <- "obj_oo.rds"
  }
}
if (!file.exists(in_rds)) stop("Missing scenic_downstream_input_rds: ", in_rds, call. = FALSE)

stage_col <- if (!is.null(cfg$scenic_stage_column)) as.character(cfg$scenic_stage_column) else "stage"
exclude_stage <- if (!is.null(cfg$scenic_exclude_stage_value)) as.character(cfg$scenic_exclude_stage_value) else "MII"

rss_var <- if (!is.null(cfg$scenic_rss_var)) as.character(cfg$scenic_rss_var) else "pseudotime_bin_q"
rss_keep_levels <- cfg$scenic_rss_keep_levels  # optional list like ["1","2",...]
rss_zThreshold <- if (!is.null(cfg$scenic_rss_zThreshold)) as.numeric(cfg$scenic_rss_zThreshold) else 1
rss_thr <- if (!is.null(cfg$scenic_rss_thr)) as.numeric(cfg$scenic_rss_thr) else 0.1
rss_cluster_columns <- if (!is.null(cfg$scenic_rss_cluster_columns)) as.logical(cfg$scenic_rss_cluster_columns) else FALSE
rss_order_rows <- if (!is.null(cfg$scenic_rss_order_rows)) as.logical(cfg$scenic_rss_order_rows) else TRUE
col_low <- if (!is.null(cfg$scenic_rss_col_low)) as.character(cfg$scenic_rss_col_low) else "#005887"
col_mid <- if (!is.null(cfg$scenic_rss_col_mid)) as.character(cfg$scenic_rss_col_mid) else "#926987"
col_high <- if (!is.null(cfg$scenic_rss_col_high)) as.character(cfg$scenic_rss_col_high) else "#ED7B86"
revCol <- if (!is.null(cfg$scenic_rss_revCol)) as.logical(cfg$scenic_rss_revCol) else FALSE
verbose <- if (!is.null(cfg$scenic_rss_verbose)) as.logical(cfg$scenic_rss_verbose) else TRUE

# ---- Outputs ----
base_out_dir <- if (!is.null(cfg$scenic_out_dir)) cfg$scenic_out_dir else "results/Fig3/scenic"
out_dir <- if (!is.null(cfg$scenic_downstream_out_dir)) cfg$scenic_downstream_out_dir else file.path(base_out_dir, "downstream")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rss_csv <- if (!is.null(cfg$scenic_rss_csv)) cfg$scenic_rss_csv else file.path(out_dir, "rss_matrix.csv")
out_pdf <- if (!is.null(cfg$scenic_rss_plot_pdf)) cfg$scenic_rss_plot_pdf else file.path(out_dir, "rss_plot.pdf")
out_png <- if (!is.null(cfg$scenic_rss_plot_png)) cfg$scenic_rss_plot_png else file.path(out_dir, "rss_plot.png")
out_regulons_rds <- if (!is.null(cfg$scenic_regulons_rds)) cfg$scenic_regulons_rds else file.path(out_dir, "regulons_list.rds")
out_thr_rds <- if (!is.null(cfg$scenic_auc_thresholds_rds)) cfg$scenic_auc_thresholds_rds else file.path(out_dir, "regulon_auc_thresholds.rds")

message("[fig3-scenic-down] Open loom: ", loom_path)
loom <- open_loom(loom_path)

regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
saveRDS(regulons, out_regulons_rds)

regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
saveRDS(regulonAucThresholds, out_thr_rds)

# Load Seurat object and filter stage
message("[fig3-scenic-down] Load Seurat: ", in_rds)
obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Expected Seurat object in scenic_downstream_input_rds", call. = FALSE)
if (!stage_col %in% colnames(obj@meta.data)) stop("meta.data missing stage column: ", stage_col, call. = FALSE)

obj <- subset(obj, subset = get(stage_col) != exclude_stage)

# Match cells between Seurat and AUC
common_cells <- intersect(colnames(obj), colnames(regulonAUC))
if (length(common_cells) == 0) stop("No overlapping cells between Seurat and loom regulonAUC.", call. = FALSE)

obj <- subset(obj, cells = common_cells)

sub_regulonAUC <- regulonAUC[, match(colnames(obj), colnames(regulonAUC))]
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ]

if (!rss_var %in% colnames(obj@meta.data)) {
  stop("meta.data missing RSS annotation column: ", rss_var, call. = FALSE)
}

cell_ann <- data.frame(row.names = colnames(obj), tmp = obj@meta.data[[rss_var]], stringsAsFactors = FALSE)
colnames(cell_ann) <- rss_var

# Ensure the annotation is a factor with ordered levels if it looks numeric
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

message("[fig3-scenic-down] calcRSS by: ", rss_var)
rss <- calcRSS(AUC = getAUC(sub_regulonAUC),
               cellAnnotation = cell_ann[colnames(sub_regulonAUC), rss_var, drop = FALSE])
rss <- na.omit(rss)

# Keep selected bins/levels (default: first 10 numeric levels)
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

# Save RSS matrix
rss_df <- as.data.frame(rss)
rss_df$regulon <- rownames(rss_df)
rss_df <- rss_df[, c("regulon", keep), drop = FALSE]
write.csv(rss_df, out_rss_csv, row.names = FALSE)

# Plot
message("[fig3-scenic-down] plotRSS -> ", out_pdf)
rssPlot <- plotRSS(rss,
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
                   verbose = verbose)

ggsave(out_pdf, plot = rssPlot, width = 10, height = 8, units = "in")
ggsave(out_png, plot = rssPlot, width = 10, height = 8, units = "in", dpi = 300)

# Close loom
try({
  close_loom(loom)
}, silent = TRUE)

# Record params + session
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
