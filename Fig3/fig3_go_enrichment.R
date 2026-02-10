#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Fig3 | GO enrichment for merged Monocle3 gene modules
#
# What this script does
#   - Reads gene module assignments produced by fig3_monocle3_modules.R.
#   - Merges modules according to `go_merge_map` and runs GO enrichment (clusterProfiler).
#   - Optionally simplifies redundant GO terms and writes a merged result table.
#
# Recommended way to run (repo wrapper)
#   bash Fig3/run_fig3_combined.sh --only go -o <out_dir>
#
# Run this script directly (module-level YAML)
#   Rscript Fig3/fig3_go_enrichment.R --config <fig3_module.yaml>
#
# Key config fields (module-level YAML)
#   out_dir, out_gene_modules_csv, go_out_csv,
#   go_keyType, go_ont, go_pAdjustMethod, go_pvalueCutoff, go_qvalueCutoff,
#   go_simplify_cutoff, go_simplify_by, go_simplify_select_fun, go_merge_map
#
# Outputs
#   - GO_BP_merged_modules_simplified.csv (or go_out_csv)
#   - params_used_fig3_go.yaml, sessionInfo_fig3_go.txt
#
# Dependencies
#   clusterProfiler, org.Mm.eg.db (Bioconductor), dplyr, yaml
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
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
"fig3_go_enrichment.R

Required:
  --config <yaml>     YAML config path (reuse Fig3 config)

Example:
  Rscript fig3_go_enrichment.R --config Fig3/configs/fig3_monocle3.yaml
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

# Convert YAML list mapping into merge_map list of character vectors
# Expected YAML structure:
# go_merge_map:
#   M1: ["3"]
#   M2: ["1"]
yaml_merge_map_to_list <- function(x) {
  if (is.null(x)) return(NULL)
  out <- list()
  for (nm in names(x)) {
    out[[nm]] <- as.character(unlist(x[[nm]]))
  }
  out
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

need_yaml()
cfg <- yaml::read_yaml(args$config)

out_dir <- if (!is.null(cfg$out_dir)) cfg$out_dir else "results/Fig3"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

modules_csv <- if (!is.null(cfg$out_gene_modules_csv)) cfg$out_gene_modules_csv else file.path(out_dir, "gene_modules.csv")
if (!file.exists(modules_csv)) stop("Missing gene_modules.csv: ", modules_csv, call. = FALSE)

go_out_csv <- if (!is.null(cfg$go_out_csv)) cfg$go_out_csv else file.path(out_dir, "GO_BP_merged_modules_simplified.csv")

# Enrichment params
keyType <- if (!is.null(cfg$go_keyType)) as.character(cfg$go_keyType) else "SYMBOL"
ont <- if (!is.null(cfg$go_ont)) as.character(cfg$go_ont) else "BP"
pAdjustMethod <- if (!is.null(cfg$go_pAdjustMethod)) as.character(cfg$go_pAdjustMethod) else "BH"
pvalueCutoff <- if (!is.null(cfg$go_pvalueCutoff)) as.numeric(cfg$go_pvalueCutoff) else 0.05
qvalueCutoff <- if (!is.null(cfg$go_qvalueCutoff)) as.numeric(cfg$go_qvalueCutoff) else 0.2
simplify_cutoff <- if (!is.null(cfg$go_simplify_cutoff)) as.numeric(cfg$go_simplify_cutoff) else 0.7
simplify_by <- if (!is.null(cfg$go_simplify_by)) as.character(cfg$go_simplify_by) else "p.adjust"
simplify_select_fun <- if (!is.null(cfg$go_simplify_select_fun)) as.character(cfg$go_simplify_select_fun) else "min"

merge_map <- yaml_merge_map_to_list(cfg$go_merge_map)

# Default merge map (your latest)
if (is.null(merge_map)) {
  merge_map <- list(
    M1 = c("3"),
    M2 = c("1"),
    M3 = c("4"),
    M4 = c("5"),
    M5 = c("2")
  )
}

message("[fig3-go] Loading modules: ", modules_csv)
gene_module_df <- read.csv(modules_csv, stringsAsFactors = FALSE)

required_cols <- c("id", "module")
if (!all(required_cols %in% colnames(gene_module_df))) {
  stop("gene_modules.csv must contain columns: id, module", call. = FALSE)
}

# Extract genes per module
modules <- split(gene_module_df$id, as.character(gene_module_df$module))

# Merge modules by merge_map
merged_modules <- lapply(merge_map, function(group) {
  genes <- unlist(modules[group], use.names = FALSE)
  genes <- genes[!is.na(genes)]
  unique(as.character(genes))
})

# Run GO enrichment
ego_all <- list()

for (m in names(merged_modules)) {
  message("[fig3-go] Processing merged module ", m, " ...")
  genes <- merged_modules[[m]]
  if (length(genes) == 0) next

  ego <- tryCatch({
    enrichGO(
      gene          = genes,
      OrgDb         = org.Mm.eg.db,
      keyType       = keyType,
      ont           = ont,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff,
      readable      = TRUE
    )
  }, error = function(e) NULL)

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_simple <- tryCatch({
      # simplify_select_fun: "min" or "max"
      sf <- if (simplify_select_fun == "max") max else min
      simplify(ego, cutoff = simplify_cutoff, by = simplify_by, select_fun = sf)
    }, error = function(e) ego)

    df <- as.data.frame(ego_simple)
    df$merged_module <- m
    ego_all[[m]] <- df
  }
}

ego_all_df <- bind_rows(ego_all)
write.csv(ego_all_df, file = go_out_csv, row.names = FALSE)

message("[fig3-go] ✅ Enrichment done. Results: ", nrow(ego_all_df), " rows -> ", go_out_csv)

# Reproducibility artifacts
params_used <- list(
  out_dir = out_dir,
  modules_csv = modules_csv,
  go_out_csv = go_out_csv,
  go_keyType = keyType,
  go_ont = ont,
  go_pAdjustMethod = pAdjustMethod,
  go_pvalueCutoff = pvalueCutoff,
  go_qvalueCutoff = qvalueCutoff,
  go_simplify_cutoff = simplify_cutoff,
  go_simplify_by = simplify_by,
  go_simplify_select_fun = simplify_select_fun,
  go_merge_map = merge_map
)

yaml::write_yaml(params_used, file.path(out_dir, "params_used_fig3_go.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig3_go.txt"), capture.output(sessionInfo()))
