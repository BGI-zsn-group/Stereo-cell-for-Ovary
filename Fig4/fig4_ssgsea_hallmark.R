\
#!/usr/bin/env Rscript
# Fig4: Hallmark ssGSEA on GC celltype-average expression (standalone module)
#
# This script:
#   1) Loads a Seurat object (.rds) that contains GC cells with a `celltype` column.
#   2) Computes average expression per group (default: group_by = "celltype") using log-normalized data.
#   3) Runs GSVA::gsva(method = "ssgsea") using MSigDB Hallmark gene sets for mouse.
#   4) Keeps a curated subset of Hallmark pathways (default list embedded; configurable).
#   5) Reorders columns by a desired celltype order and writes matrix + heatmap.
#
# Inputs (YAML; extracted from Fig4/configs/fig4_combined.yaml by run_fig4_combined.sh)
#   - input_rds: Seurat object path
#   - out_dir: output directory
#   - prefix: output prefix
#   - assay: default "RNA"
#   - group_by: default "celltype"
#   - slot: default "data" (log-normalized expression)
#   - species: default "Mus musculus"
#   - category: default "H" (Hallmark)
#   - method: default "ssgsea"
#   - kcdf: default "Gaussian" (for log-expression)
#   - min_sz / max_sz: gene set size filters
#   - parallel_sz: threads for GSVA
#   - keep_pathways: optional list of Hallmark pathway names (without "HALLMARK_" prefix)
#   - celltype_order: optional list to reorder columns
#
# Outputs (under out_dir)
#   - <prefix>.hallmark_ssgsea.all.csv            (all Hallmark pathways)
#   - <prefix>.hallmark_ssgsea.selected.csv       (selected pathways)
#   - <prefix>.hallmark_ssgsea.selected.pdf       (heatmap)
#   - <prefix>.hallmark_ssgsea.selected.rds       (matrix as R object)
#
# Usage:
#   Rscript Fig4/fig4_ssgsea_hallmark.R --config <module_config.yaml>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(GSVA)
  library(msigdbr)
  library(pheatmap)
})

help_msg <- function() {
  cat(
"Usage:
  Rscript Fig4/fig4_ssgsea_hallmark.R --config <yaml>

Options:
  --config <yaml>   Module config YAML extracted from combined YAML.

"
  )
}

# simple arg parsing
args <- commandArgs(trailingOnly = TRUE)
cfg_path <- NULL
if (length(args) == 0) { help_msg(); stop("Missing --config", call. = FALSE) }
i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") {
    if (i + 1 > length(args)) stop("Missing value for --config", call. = FALSE)
    cfg_path <- args[i + 1]
    i <- i + 2
  } else if (args[i] %in% c("-h", "--help")) {
    help_msg(); quit(status = 0)
  } else {
    stop(paste0("Unknown arg: ", args[i]), call. = FALSE)
  }
}

if (is.null(cfg_path) || !file.exists(cfg_path)) stop(paste0("Config not found: ", cfg_path), call. = FALSE)
cfg <- yaml::read_yaml(cfg_path)

`%||%` <- function(a, b) if (!is.null(a)) a else b

input_rds <- cfg$input_rds %||% stop("Missing config: input_rds", call. = FALSE)
out_dir   <- cfg$out_dir %||% "results/Fig4/ssgsea_hallmark"
prefix    <- cfg$prefix %||% "gc"
assay     <- cfg$assay %||% "RNA"
group_by  <- cfg$group_by %||% "celltype"
slot_use  <- cfg$slot %||% "data"

species   <- cfg$species %||% "Mus musculus"
category  <- cfg$category %||% "H"

method    <- cfg$method %||% "ssgsea"
kcdf      <- cfg$kcdf %||% "Gaussian"
min_sz    <- cfg$min_sz %||% 10
max_sz    <- cfg$max_sz %||% Inf
parallel_sz <- cfg$parallel_sz %||% 10

keep_pathways <- cfg$keep_pathways
celltype_order <- cfg$celltype_order

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

logf <- file.path(out_dir, paste0(prefix, ".hallmark_ssgsea.log"))
sink(logf, split = TRUE)
cat("[Info] input_rds:", input_rds, "\n")
cat("[Info] out_dir:", out_dir, "\n")
cat("[Info] prefix:", prefix, "\n")
cat("[Info] assay:", assay, " group_by:", group_by, " slot:", slot_use, "\n")
cat("[Info] GSVA method:", method, " kcdf:", kcdf, " min_sz:", min_sz, " max_sz:", max_sz, " threads:", parallel_sz, "\n")

if (!file.exists(input_rds)) stop(paste0("Input RDS not found: ", input_rds), call. = FALSE)
obj <- readRDS(input_rds)

if (!(assay %in% Assays(obj))) stop(paste0("Assay not found in object: ", assay), call. = FALSE)
DefaultAssay(obj) <- assay

if (!(group_by %in% colnames(obj@meta.data))) {
  stop(paste0("group_by column not found in meta.data: ", group_by), call. = FALSE)
}

# Average expression per group
avg_expr_list <- AverageExpression(
  obj,
  group.by = group_by,
  assays = assay,
  slot = slot_use
)
avg_expr_mat <- avg_expr_list[[assay]]
cat("[Info] avg_expr_mat:", nrow(avg_expr_mat), "genes x", ncol(avg_expr_mat), "groups\n")

# Hallmark gene sets (mouse)
hallmark_msig <- msigdbr(species = species, category = category)
hallmark_list <- split(hallmark_msig$gene_symbol, hallmark_msig$gs_name)

# Run GSVA/ssGSEA with compatibility for different GSVA versions
run_gsva <- function(expr_mat, gsets, method, kcdf, min_sz, max_sz, parallel_sz) {
  # Old API
  out <- tryCatch({
    GSVA::gsva(
      expr = expr_mat,
      gset.idx.list = gsets,
      method = method,
      kcdf = kcdf,
      min.sz = min_sz,
      max.sz = max_sz,
      parallel.sz = parallel_sz
    )
  }, error = function(e) {
    msg <- conditionMessage(e)
    # Newer GSVA APIs may require param objects
    if (grepl("Param|param", msg, ignore.case = TRUE)) {
      if (method == "ssgsea" && exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
        param <- get("ssgseaParam", envir = asNamespace("GSVA"))(expr_mat, gsets, kcdf = kcdf, minSize = min_sz, maxSize = max_sz)
        return(GSVA::gsva(param, verbose = FALSE))
      }
      if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
        param <- get("gsvaParam", envir = asNamespace("GSVA"))(expr_mat, gsets, method = method, kcdf = kcdf, minSize = min_sz, maxSize = max_sz)
        return(GSVA::gsva(param, verbose = FALSE))
      }
    }
    stop(e)
  })
  out
}

gsva_res <- run_gsva(avg_expr_mat, hallmark_list, method, kcdf, min_sz, max_sz, parallel_sz)

# Remove HALLMARK_ prefix in rownames
rownames(gsva_res) <- sub("^HALLMARK_", "", rownames(gsva_res))

# Save all pathways
all_csv <- file.path(out_dir, paste0(prefix, ".hallmark_ssgsea.all.csv"))
write.csv(gsva_res, all_csv, quote = FALSE)
cat("[Info] wrote:", all_csv, "\n")

# Default curated list (26 pathways; names without HALLMARK_)
default_keep <- c(
  "E2F_TARGETS","G2M_CHECKPOINT","MITOTIC_SPINDLE","MYC_TARGETS_V1","MYC_TARGETS_V2","DNA_REPAIR",
  "WNT_BETA_CATENIN_SIGNALING","NOTCH_SIGNALING","HEDGEHOG_SIGNALING","TGF_BETA_SIGNALING","ANGIOGENESIS",
  "OXIDATIVE_PHOSPHORYLATION","GLYCOLYSIS","FATTY_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS","PEROXISOME",
  "MTORC1_SIGNALING","PI3K_AKT_MTOR_SIGNALING",
  "P53_PATHWAY","HYPOXIA","REACTIVE_OXYGEN_SPECIES_PATHWAY","UNFOLDED_PROTEIN_RESPONSE",
  "APOPTOSIS","EPITHELIAL_MESENCHYMAL_TRANSITION",
  "TNFA_SIGNALING_VIA_NFKB","IL6_JAK_STAT3_SIGNALING"
)

keep <- keep_pathways %||% default_keep
keep_in <- intersect(keep, rownames(gsva_res))
miss <- setdiff(keep, rownames(gsva_res))
if (length(miss) > 0) cat("[Warn] pathways not found:", paste(miss, collapse = ", "), "\n")

mat_sel <- gsva_res[keep_in, , drop = FALSE]

# Reorder columns if requested
if (!is.null(celltype_order)) {
  celltype_order <- intersect(celltype_order, colnames(mat_sel))
  if (length(celltype_order) > 0) {
    mat_sel <- mat_sel[, celltype_order, drop = FALSE]
  }
}

sel_csv <- file.path(out_dir, paste0(prefix, ".hallmark_ssgsea.selected.csv"))
write.csv(mat_sel, sel_csv, quote = FALSE)
cat("[Info] wrote:", sel_csv, "\n")

sel_rds <- file.path(out_dir, paste0(prefix, ".hallmark_ssgsea.selected.rds"))
saveRDS(mat_sel, sel_rds)
cat("[Info] wrote:", sel_rds, "\n")

# Heatmap
pdf_path <- file.path(out_dir, paste0(prefix, ".hallmark_ssgsea.selected.pdf"))
pdf(pdf_path, width = 9, height = max(4, 0.25 * nrow(mat_sel) + 2))
pheatmap::pheatmap(mat_sel)
dev.off()
cat("[Info] wrote:", pdf_path, "\n")

cat("[OK] Done.\n")
sink()
