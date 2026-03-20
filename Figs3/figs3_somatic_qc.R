#!/usr/bin/env Rscript
# Figs3: Somatic-cell QC with support for multiple input Seurat RDS files.
# Each input can specify a `sample` label. QC is performed independently per RDS,
# then all QC-passed objects are merged. The merged object retains a single
# metadata column `sample` carrying the user-provided labels.

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

Accepted config modes:
  1) Legacy single-input:
       qc_input_rds: /path/to/sample.rds
  2) Multi-input:
       qc_inputs:
         - path: /path/PD14_1.rds
           sample: PD14-1
         - path: /path/PD14_2.rds
           sample: PD14-2

Example:
  Rscript figs3_somatic_qc.R --config Figs3/configs/figs3_qc.yaml
", sep = "")
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

safe_tag <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "sample"
  x
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

normalize_qc_inputs <- function(cfg) {
  qci <- cfg$qc_inputs
  if (!is.null(qci) && length(qci) > 0) {
    if (!is.list(qci)) stop("qc_inputs must be a list.", call. = FALSE)

    # Support both YAML sequences and dotted CLI overrides such as
    # qc_inputs.0.path=... qc_inputs.0.sample=...
    if (!is.null(names(qci)) && all(grepl("^[0-9]+$", names(qci)))) {
      ord <- order(as.integer(names(qci)))
      qci <- qci[ord]
    }

    out <- list()
    idx <- 1
    for (it in qci) {
      if (is.null(it)) next
      path <- it$path %||% it$rds %||% it$input_rds
      sample <- it$sample %||% it$sample_label %||% it$name
      if (is.null(path) || !nzchar(as.character(path))) stop("Each qc_inputs item must provide `path`.", call. = FALSE)
      if (is.null(sample) || !nzchar(as.character(sample))) stop("Each qc_inputs item must provide `sample`.", call. = FALSE)
      out[[idx]] <- list(path = as.character(path), sample = as.character(sample))
      idx <- idx + 1
    }
    if (length(out) == 0) stop("qc_inputs is present but empty after normalization.", call. = FALSE)
    return(out)
  }

  in_rds <- cfg$qc_input_rds
  if (is.null(in_rds) || !nzchar(as.character(in_rds))) stop("Provide either qc_input_rds or qc_inputs.", call. = FALSE)
  sample <- cfg$qc_input_sample %||% cfg$qc_single_sample %||% tools::file_path_sans_ext(basename(as.character(in_rds)))
  list(list(path = as.character(in_rds), sample = as.character(sample)))
}

run_qc_one <- function(scrna0, sample_label, species, group_by, do_stage_filter, stage_col,
                       exclude_stage, mt_max, umi_iqr_mult, dbl_rate, per_sample_dir, prefix) {
  if (!inherits(scrna0, "Seurat")) stop("Expected a Seurat object.", call. = FALSE)

  if (!"orig.ident" %in% colnames(scrna0@meta.data)) {
    scrna0$orig.ident <- Idents(scrna0)
  }

  # Enforce one user-facing metadata key: sample
  scrna0$sample <- sample_label
  scrna0$orig.ident <- sample_label

  if (do_stage_filter) {
    if (!stage_col %in% colnames(scrna0@meta.data)) stop("meta.data missing stage column: ", stage_col, call. = FALSE)
    keep_cells <- rownames(scrna0@meta.data)[as.character(scrna0@meta.data[[stage_col]]) != exclude_stage]
    scrna0 <- subset(scrna0, cells = keep_cells)
    message("[figs3-qc] [", sample_label, "] stage filter: excluded ", stage_col, " == '", exclude_stage, "'. cells: ", ncol(scrna0))
  }

  message("[figs3-qc] [", sample_label, "] initial cells: ", ncol(scrna0))

  # percent.mt
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

  sample_tag <- safe_tag(sample_label)
  raw_pdf <- file.path(per_sample_dir, paste0(prefix, ".", sample_tag, ".Raw.QC.VlnPlot.pdf"))
  p1 <- VlnPlot(scrna0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = group_by, ncol = 3, pt.size = 0.01) +
    ggtitle(paste0("Raw QC | ", sample_label)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(raw_pdf, p1, width = 12, height = 5)

  # MT filter
  scrna1 <- subset(scrna0, subset = percent.mt < mt_max)
  message("[figs3-qc] [", sample_label, "] after MT filter cells: ", ncol(scrna1))

  mt_pdf <- file.path(per_sample_dir, paste0(prefix, ".", sample_tag, ".MTQCed.VlnPlot.pdf"))
  p2 <- VlnPlot(scrna1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = group_by, ncol = 3, pt.size = 0.1) +
    ggtitle(paste0("After MT QC | ", sample_label)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(mt_pdf, p2, width = 12, height = 5)

  # UMI upper bound
  Q1 <- quantile(scrna1$nCount_RNA, probs = 0.25, na.rm = TRUE)
  Q3 <- quantile(scrna1$nCount_RNA, probs = 0.75, na.rm = TRUE)
  IQRv <- Q3 - Q1
  upper <- Q3 + IQRv * umi_iqr_mult
  message("[figs3-qc] [", sample_label, "] UMI upper: ", as.numeric(upper))

  scrna2 <- subset(scrna1, subset = nCount_RNA < upper)
  message("[figs3-qc] [", sample_label, "] after UMI QC cells: ", ncol(scrna2))

  umi_pdf <- file.path(per_sample_dir, paste0(prefix, ".", sample_tag, ".UMIQCed.VlnPlot.pdf"))
  p3 <- VlnPlot(scrna2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = group_by, ncol = 3, pt.size = 0.01) +
    ggtitle(paste0("After UMI QC | ", sample_label)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(umi_pdf, p3, width = 12, height = 5)

  # Doublet detection
  message("[figs3-qc] [", sample_label, "] run scDblFinder (dbr=", dbl_rate, ") ...")
  sce <- as.SingleCellExperiment(scrna2)
  sce <- scDblFinder(sce, dbr = dbl_rate)

  scrna2$scDblFinder_score <- sce$scDblFinder.score
  scrna2$scDblFinder_class <- sce$scDblFinder.class

  dbl_scores_path <- file.path(per_sample_dir, paste0(prefix, ".", sample_tag, ".doublet_scores.tsv"))
  doublet_df <- data.frame(
    Cell = colnames(scrna2),
    Doublet_Score = scrna2$scDblFinder_score,
    Doublet_Class = scrna2$scDblFinder_class,
    Group = scrna2@meta.data[[group_by]],
    stringsAsFactors = FALSE
  )
  write.table(doublet_df, file = dbl_scores_path, sep = "\t", quote = FALSE, row.names = FALSE)

  message("[figs3-qc] [", sample_label, "] doublet class table:")
  print(table(scrna2$scDblFinder_class))

  # Filter singlets
  scrna3 <- subset(scrna2, subset = scDblFinder_class == "singlet")
  scrna3$sample <- sample_label
  scrna3$orig.ident <- sample_label
  message("[figs3-qc] [", sample_label, "] final cells: ", ncol(scrna3))

  final_pdf <- file.path(per_sample_dir, paste0(prefix, ".", sample_tag, ".Final.QC.VlnPlot.pdf"))
  p4 <- VlnPlot(scrna3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = group_by, ncol = 3, pt.size = 0.1) +
    ggtitle(paste0("After Doublet Filter | ", sample_label)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(final_pdf, p4, width = 12, height = 5)

  stats <- list(
    sample = sample_label,
    initial_cells = ncol(scrna0),
    after_mt_cells = ncol(scrna1),
    after_umi_cells = ncol(scrna2),
    final_cells = ncol(scrna3),
    umi_upper = as.numeric(upper),
    outputs = list(
      raw_qc_pdf = raw_pdf,
      mt_qc_pdf = mt_pdf,
      umi_qc_pdf = umi_pdf,
      final_qc_pdf = final_pdf,
      doublet_scores = dbl_scores_path
    )
  )

  list(obj = scrna3, stats = stats)
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

species <- cfg$qc_species %||% "mouse"
species <- as.character(species)

out_dir <- as.character(cfg$qc_out_dir %||% "results/Figs3/qc")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- as.character(cfg$qc_prefix %||% "somatic")
group_by <- as.character(cfg$qc_group_by %||% "sample")

stage_col <- cfg$qc_stage_column
exclude_stage <- cfg$qc_exclude_stage_value
do_stage_filter <- (!is.null(stage_col) && nzchar(as.character(stage_col)) &&
                    !is.null(exclude_stage) && nzchar(as.character(exclude_stage)))
if (!is.null(stage_col)) stage_col <- as.character(stage_col)
if (!is.null(exclude_stage)) exclude_stage <- as.character(exclude_stage)

mt_max <- as.numeric(cfg$qc_percent_mt_max %||% 20)
umi_iqr_mult <- as.numeric(cfg$qc_umi_iqr_multiplier %||% 1.5)
dbl_rate <- as.numeric(cfg$qc_doublet_rate %||% 0.05)

out_rds <- as.character(cfg$qc_out_rds %||% file.path(out_dir, paste0(prefix, ".qc_processed.rds")))
if (!dir.exists(dirname(out_rds))) dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, paste0(prefix, ".qc.log"))
log_file <- file(log_path, open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

message("[figs3-qc] out_dir: ", out_dir)
message("[figs3-qc] species: ", species)
message("[figs3-qc] prefix: ", prefix)
message("[figs3-qc] start: ", Sys.time())

qc_inputs <- normalize_qc_inputs(cfg)
message("[figs3-qc] n_inputs: ", length(qc_inputs))

per_sample_dir <- file.path(out_dir, "per_sample_qc")
if (!dir.exists(per_sample_dir)) dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)

qc_objs <- list()
qc_stats <- list()

for (i in seq_along(qc_inputs)) {
  path_i <- qc_inputs[[i]]$path
  sample_i <- qc_inputs[[i]]$sample
  if (!file.exists(path_i)) stop("Missing qc input: ", path_i, call. = FALSE)
  message("[figs3-qc] [", sample_i, "] input_rds: ", path_i)
  obj_i <- readRDS(path_i)
  res_i <- run_qc_one(
    scrna0 = obj_i,
    sample_label = sample_i,
    species = species,
    group_by = group_by,
    do_stage_filter = do_stage_filter,
    stage_col = stage_col,
    exclude_stage = exclude_stage,
    mt_max = mt_max,
    umi_iqr_mult = umi_iqr_mult,
    dbl_rate = dbl_rate,
    per_sample_dir = per_sample_dir,
    prefix = prefix
  )
  qc_objs[[sample_i]] <- res_i$obj
  qc_stats[[sample_i]] <- res_i$stats
}

merged <- NULL
if (length(qc_objs) == 1) {
  merged <- qc_objs[[1]]
} else {
  samples <- names(qc_objs)
  merged <- qc_objs[[samples[[1]]]]
  merged <- merge(merged, y = qc_objs[samples[-1]], add.cell.ids = samples)
}
merged$sample <- factor(as.character(merged$sample))
merged$orig.ident <- factor(as.character(merged$sample))

saveRDS(merged, file = out_rds)
message("[figs3-qc] saved merged QC object: ", out_rds)

merged_pdf <- file.path(out_dir, paste0(prefix, ".Merged.Final.QC.VlnPlot.pdf"))
p_merged <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    group.by = "sample", ncol = 3, pt.size = 0.1) +
  ggtitle("Merged QC object") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(merged_pdf, p_merged, width = 12, height = 5)

stats <- list(
  n_inputs = length(qc_inputs),
  samples = lapply(qc_stats, identity),
  merged_final_cells = ncol(merged)
)
yaml::write_yaml(stats, file.path(out_dir, paste0(prefix, ".qc_stats.yaml")))

params_used <- list(
  qc_enabled = enabled,
  qc_species = species,
  qc_inputs = qc_inputs,
  qc_out_dir = out_dir,
  qc_prefix = prefix,
  qc_group_by = group_by,
  qc_stage_column = if (do_stage_filter) stage_col else NULL,
  qc_exclude_stage_value = if (do_stage_filter) exclude_stage else NULL,
  qc_percent_mt_max = mt_max,
  qc_umi_iqr_multiplier = umi_iqr_mult,
  qc_doublet_rate = dbl_rate,
  outputs = list(
    out_rds = out_rds,
    log = log_path,
    per_sample_dir = per_sample_dir,
    merged_final_qc_pdf = merged_pdf
  )
)
yaml::write_yaml(params_used, file.path(out_dir, "params_used_figs3_somatic_qc.yaml"))
write_text(file.path(out_dir, "sessionInfo_figs3_somatic_qc.txt"), capture.output(sessionInfo()))

message("[figs3-qc] done: ", Sys.time())

sink()
sink(type = "message")
close(log_file)
