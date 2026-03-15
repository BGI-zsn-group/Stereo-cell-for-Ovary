#!/usr/bin/env Rscript
# Figs3: Somatic-cell QC for Seurat object (standalone figure)
#
# This script is intentionally standalone (NOT part of Fig3 pipeline).
# It performs:
#   - percent.mt calculation (mouse/human/Axolotl)
#   - mitochondrial filtering
#   - UMI upper-bound filtering (IQR rule)
#   - doublet detection via scDblFinder
# and saves a QC-processed Seurat object plus QC plots and logs.
#
# Usage:
#   Rscript figs3_somatic_qc.R --config Figs3/configs/figs3_qc.yaml
#
suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(Matrix)
  library(scDblFinder)
  library(SingleCellExperiment)
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
"figs3_somatic_qc.R

Required:
  --config <yaml>     YAML config path

Example:
  Rscript figs3_somatic_qc.R --config Figs3/configs/figs3_qc.yaml
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

enabled <- cfg$qc_enabled
if (is.null(enabled)) enabled <- TRUE
enabled <- as.logical(enabled)
if (!enabled) {
  message("[figs3-qc] qc_enabled=false; skip.")
  quit(status = 0)
}

# ---- Config ----
species <- cfg$qc_species
if (is.null(species) || !nzchar(as.character(species))) species <- "mouse"
species <- as.character(species)

in_rds <- cfg$qc_input_rds
if (is.null(in_rds) || !nzchar(as.character(in_rds))) stop("qc_input_rds is required for Figs3.", call. = FALSE)
in_rds <- as.character(in_rds)
if (!file.exists(in_rds)) stop("Missing qc_input_rds: ", in_rds, call. = FALSE)

out_dir <- cfg$qc_out_dir
if (is.null(out_dir) || !nzchar(as.character(out_dir))) out_dir <- "results/Figs3/qc"
out_dir <- as.character(out_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- cfg$qc_prefix
if (is.null(prefix) || !nzchar(as.character(prefix))) prefix <- "somatic"
prefix <- as.character(prefix)

group_by <- cfg$qc_group_by
if (is.null(group_by) || !nzchar(as.character(group_by))) group_by <- "orig.ident"
group_by <- as.character(group_by)

# Optional stage filtering
stage_col <- cfg$qc_stage_column
exclude_stage <- cfg$qc_exclude_stage_value
do_stage_filter <- (!is.null(stage_col) && nzchar(as.character(stage_col)) &&
                    !is.null(exclude_stage) && nzchar(as.character(exclude_stage)))
if (!is.null(stage_col)) stage_col <- as.character(stage_col)
if (!is.null(exclude_stage)) exclude_stage <- as.character(exclude_stage)

mt_max <- cfg$qc_percent_mt_max
if (is.null(mt_max)) mt_max <- 20
mt_max <- as.numeric(mt_max)

umi_iqr_mult <- cfg$qc_umi_iqr_multiplier
if (is.null(umi_iqr_mult)) umi_iqr_mult <- 1.5
umi_iqr_mult <- as.numeric(umi_iqr_mult)

dbl_rate <- cfg$qc_doublet_rate
if (is.null(dbl_rate)) dbl_rate <- 0.05
dbl_rate <- as.numeric(dbl_rate)

out_rds <- cfg$qc_out_rds
if (is.null(out_rds) || !nzchar(as.character(out_rds))) out_rds <- file.path(out_dir, paste0(prefix, ".qc_processed.rds"))
out_rds <- as.character(out_rds)

# ---- Logging ----
log_path <- file.path(out_dir, paste0(prefix, ".qc.log"))
log_file <- file(log_path, open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

message("[figs3-qc] out_dir: ", out_dir)
message("[figs3-qc] species: ", species)
message("[figs3-qc] input_rds: ", in_rds)
message("[figs3-qc] prefix: ", prefix)
message("[figs3-qc] start: ", Sys.time())

# ---- Load ----
scrna0 <- readRDS(in_rds)
if (!inherits(scrna0, "Seurat")) stop("Expected a Seurat object in qc_input_rds.", call. = FALSE)
message("[figs3-qc] initial cells: ", ncol(scrna0))

# Ensure orig.ident exists for grouping plots
if (!"orig.ident" %in% colnames(scrna0@meta.data)) {
  scrna0$orig.ident <- Idents(scrna0)
}

# Optional stage filter
if (do_stage_filter) {
  if (!stage_col %in% colnames(scrna0@meta.data)) stop("meta.data missing stage column: ", stage_col, call. = FALSE)
  scrna0 <- subset(scrna0, subset = get(stage_col) != exclude_stage)
  message("[figs3-qc] stage filter: excluded ", stage_col, " == '", exclude_stage, "'. cells: ", ncol(scrna0))
}

# ---- 1) percent.mt ----
message("[figs3-qc] compute percent.mt ...")
if (species == "mouse") {
  scrna0[["percent.mt"]] <- PercentageFeatureSet(scrna0, pattern = "^mt-")
} else if (species == "human") {
  scrna0[["percent.mt"]] <- PercentageFeatureSet(scrna0, pattern = "^MT-")
} else if (species == "Axolotl") {
  MT_gene <- c("AMEX60DD020311","AMEX60DD020312","AMEX60DD026192","AMEX60DD038501",
               "AMEX60DD040676","AMEX60DD041809","AMEX60DD041810","AMEX60DD041811",
               "AMEX60DD041812","AMEX60DD053672")
  counts_matrix <- GetAssayData(scrna0, assay = "RNA", slot = "counts")
  exp_MTgene <- unique(MT_gene[MT_gene %in% rownames(counts_matrix)])
  scrna0[["percent.mt"]] <- PercentageFeatureSet(scrna0, features = exp_MTgene)
} else {
  stop("Unsupported species. Use mouse/human/Axolotl", call. = FALSE)
}

raw_pdf <- file.path(out_dir, paste0(prefix, ".Raw.QC.VlnPlot.pdf"))
p1 <- VlnPlot(scrna0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = group_by, ncol = 3, pt.size = 0.01) +
  ggtitle("Raw QC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(raw_pdf, p1, width = 12, height = 5)

# ---- 2) MT filter ----
message("[figs3-qc] filter percent.mt < ", mt_max, " ...")
scrna1 <- subset(scrna0, subset = percent.mt < mt_max)
message("[figs3-qc] after MT filter cells: ", ncol(scrna1))

mt_pdf <- file.path(out_dir, paste0(prefix, ".MTQCed.VlnPlot.pdf"))
p2 <- VlnPlot(scrna1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = group_by, ncol = 3, pt.size = 0.1) +
  ggtitle("After MT QC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(mt_pdf, p2, width = 12, height = 5)

# ---- 3) UMI QC (IQR upper bound) ----
message("[figs3-qc] compute UMI upper bound (IQR* ", umi_iqr_mult, ") ...")
Q1 <- quantile(scrna1$nCount_RNA, probs = 0.25, na.rm = TRUE)
Q3 <- quantile(scrna1$nCount_RNA, probs = 0.75, na.rm = TRUE)
IQRv <- Q3 - Q1
upper <- Q3 + IQRv * umi_iqr_mult
message("[figs3-qc] UMI upper: ", as.numeric(upper))

scrna2 <- subset(scrna1, subset = nCount_RNA < upper)
message("[figs3-qc] after UMI QC cells: ", ncol(scrna2))

umi_pdf <- file.path(out_dir, paste0(prefix, ".UMIQCed.VlnPlot.pdf"))
p3 <- VlnPlot(scrna2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = group_by, ncol = 3, pt.size = 0.01) +
  ggtitle("After UMI QC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(umi_pdf, p3, width = 12, height = 5)

# ---- 4) Doublet detection ----
message("[figs3-qc] run scDblFinder (dbr=", dbl_rate, ") ...")
sce <- as.SingleCellExperiment(scrna2)
sce <- scDblFinder(sce, dbr = dbl_rate)

scrna2$scDblFinder_score <- sce$scDblFinder.score
scrna2$scDblFinder_class <- sce$scDblFinder.class

dbl_scores_path <- file.path(out_dir, paste0(prefix, ".doublet_scores.tsv"))
doublet_df <- data.frame(
  Cell = colnames(scrna2),
  Doublet_Score = scrna2$scDblFinder_score,
  Doublet_Class = scrna2$scDblFinder_class,
  Group = scrna2@meta.data[[group_by]],
  stringsAsFactors = FALSE
)
write.table(doublet_df, file = dbl_scores_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("[figs3-qc] doublet class table:")
print(table(scrna2$scDblFinder_class))

# ---- 5) Filter singlets ----
message("[figs3-qc] filter singlets ...")
scrna3 <- subset(scrna2, subset = scDblFinder_class == "singlet")
message("[figs3-qc] final cells: ", ncol(scrna3))

final_pdf <- file.path(out_dir, paste0(prefix, ".Final.QC.VlnPlot.pdf"))
p4 <- VlnPlot(scrna3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = group_by, ncol = 3, pt.size = 0.1) +
  ggtitle("After Doublet Filter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(final_pdf, p4, width = 12, height = 5)

# ---- Save ----
saveRDS(scrna3, file = out_rds)
message("[figs3-qc] saved: ", out_rds)

# ---- Stats ----
stats <- list(
  initial_cells = ncol(scrna0),
  after_mt_cells = ncol(scrna1),
  after_umi_cells = ncol(scrna2),
  final_cells = ncol(scrna3),
  umi_upper = as.numeric(upper)
)
yaml::write_yaml(stats, file.path(out_dir, paste0(prefix, ".qc_stats.yaml")))

# ---- Provenance ----
params_used <- list(
  qc_enabled = enabled,
  qc_species = species,
  qc_input_rds = in_rds,
  qc_out_dir = out_dir,
  qc_prefix = prefix,
  qc_group_by = group_by,
  qc_stage_column = if (do_stage_filter) stage_col else NULL,
  qc_exclude_stage_value = if (do_stage_filter) exclude_stage else NULL,
  qc_percent_mt_max = mt_max,
  qc_umi_iqr_multiplier = umi_iqr_mult,
  qc_umi_upper = as.numeric(upper),
  qc_doublet_rate = dbl_rate,
  outputs = list(
    out_rds = out_rds,
    log = log_path,
    raw_qc_pdf = raw_pdf,
    mt_qc_pdf = mt_pdf,
    umi_qc_pdf = umi_pdf,
    final_qc_pdf = final_pdf,
    doublet_scores = dbl_scores_path
  )
)
yaml::write_yaml(params_used, file.path(out_dir, "params_used_figs3_somatic_qc.yaml"))
write_text(file.path(out_dir, "sessionInfo_figs3_somatic_qc.txt"), capture.output(sessionInfo()))

message("[figs3-qc] done: ", Sys.time())

sink()
sink(type = "message")
close(log_file)
