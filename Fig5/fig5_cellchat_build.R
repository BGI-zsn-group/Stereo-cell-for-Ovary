\
#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Fig5 | Build CellChat objects across oocyte stages + granulosa cell states
#
# This script merges:
#   - a GC Seurat object (obj_gr) and an oocyte Seurat object (obj_oo),
# then constructs CellChat objects for multiple stage-specific subsets:
#   EGO, GO1, GO2, GO3, FGO
#
# It supports two "modes" (same code path):
#   1) mode = "all"              -> CellChatDB.mouse (all interaction categories)
#   2) mode = "cellcell_contact" -> subsetDB(CellChatDB.mouse, search="Cell-Cell Contact")
#
# Key features:
#   - Optionally re-annotate GC `celltype` from cluster IDs (configurable mapping).
#   - Build stage-specific merged subsets, run full CellChat pipeline, and (optionally)
#     export ligand-receptor tables with interaction "annotation".
#   - liftCellChat() to a unified group level ordering, then save stage-wise CellChat RDS.
#
# Usage:
#   Rscript Fig5/fig5_cellchat_build.R --config <module_config.yaml>
#
# Outputs (under out_dir):
#   - cellchat_<STAGE>_<tag>.rds
#   - comm_<STAGE>_<tag>.csv          (optional; default TRUE)
#   - <prefix>.cellchat_all_<tag>.rds (combined helper)
#   - <prefix>.cellchat_build.log
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(reshape2)
  library(ggplot2)
  library(grid)
  library(magrittr)
})

help_msg <- function() {
  cat(
"Usage:
  Rscript Fig5/fig5_cellchat_build.R --config <yaml>

Options:
  --config <yaml>   Module config YAML extracted from combined YAML.
"
  )
}

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") {
    if (i + 1 > length(args)) stop("Missing value for --config", call. = FALSE)
    cfg_path <- args[i + 1]; i <- i + 2
  } else if (args[i] %in% c("-h", "--help")) {
    help_msg(); quit(status = 0)
  } else {
    stop(paste0("Unknown arg: ", args[i]), call. = FALSE)
  }
}
if (is.null(cfg_path) || !file.exists(cfg_path)) { help_msg(); stop("Config not found.", call. = FALSE) }

cfg <- yaml::read_yaml(cfg_path)
`%||%` <- function(a, b) if (!is.null(a)) a else b

obj_gr_rds <- cfg$obj_gr_rds %||% stop("Missing config: obj_gr_rds", call. = FALSE)
obj_oo_rds <- cfg$obj_oo_rds %||% stop("Missing config: obj_oo_rds", call. = FALSE)

out_dir <- cfg$out_dir %||% "results/Fig5/cellchat"
prefix  <- cfg$prefix %||% "fig5"
mode <- cfg$mode %||% "all"  # "all" or "cellcell_contact"
export_comm_tables <- cfg$export_comm_tables %||% TRUE

gc_cluster_col <- cfg$gc_cluster_col %||% "seurat_clusters"
gc_sample_col  <- cfg$gc_sample_col %||% "sample"
oo_stage_col   <- cfg$oo_stage_col %||% "stage"

gc_celltype_col <- cfg$gc_celltype_col %||% "celltype"
gc_celltype_levels <- cfg$gc_celltype_levels %||% c(
  "GC_Progenitor","GC_Preantral_1","GC_Preantral_2",
  "GC_Cumulus_1","GC_Cumulus_2","GC_Antral_Mural","GC_Atretic"
)
cluster_to_celltype <- cfg$cluster_to_celltype

oocyte_exclude_stages <- cfg$oocyte_exclude_stages %||% c("MII")
gc_exclude_celltypes  <- cfg$gc_exclude_celltypes %||% c("GC_Atretic")

stage_subsets <- cfg$stage_subsets
if (is.null(stage_subsets)) {
  if (mode == "cellcell_contact") {
    stage_subsets <- list(
      list(name="EGO", keep_labels=c("EGO","GC_Progenitor")),
      list(name="GO1", keep_labels=c("GO1","GC_Preantral_1","GC_Preantral_2","GC_Cumulus_1")),
      list(name="GO2", keep_labels=c("GO2","GC_Cumulus_1","GC_Cumulus_2")),
      list(name="GO3", keep_labels=c("GO3","GC_Cumulus_2")),
      list(name="FGO", keep_labels=c("FGO","GC_Cumulus_2"))
    )
  } else {
    stage_subsets <- list(
      list(name="EGO", keep_labels=c("EGO","GC_Progenitor")),
      list(name="GO1", keep_labels=c("GO1","GC_Preantral_1","GC_Preantral_2","GC_Cumulus_1")),
      list(name="GO2", keep_labels=c("GO2","GC_Cumulus_1","GC_Cumulus_2")),
      list(name="GO3", keep_labels=c("GO3","GC_Cumulus_2","GC_Antral_Mural")),
      list(name="FGO", keep_labels=c("FGO","GC_Cumulus_2","GC_Antral_Mural"))
    )
  }
}

group_new_levels <- cfg$group_new_levels %||% c(
  "EGO","GO1","GO2","GO3","FGO",
  "GC_Progenitor","GC_Preantral_1","GC_Preantral_2",
  "GC_Cumulus_1","GC_Cumulus_2","GC_Antral_Mural"
)

compute_type <- cfg$compute_type %||% "triMean"
population_size <- cfg$population_size %||% FALSE
seed_use <- cfg$seed_use %||% 888
min_cells <- cfg$min_cells %||% 10

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, paste0(prefix, ".cellchat_build.log"))
sink(log_file, split = TRUE)

cat("[Info] mode:", mode, "\n")
cat("[Info] obj_gr_rds:", obj_gr_rds, "\n")
cat("[Info] obj_oo_rds:", obj_oo_rds, "\n")
cat("[Info] out_dir:", out_dir, "\n")
cat("[Info] prefix:", prefix, "\n")

if (!file.exists(obj_gr_rds)) stop("obj_gr_rds not found: ", obj_gr_rds, call. = FALSE)
if (!file.exists(obj_oo_rds)) stop("obj_oo_rds not found: ", obj_oo_rds, call. = FALSE)
obj_gr <- readRDS(obj_gr_rds)
obj_oo <- readRDS(obj_oo_rds)

# Optional: re-annotate GC celltype from clusters
if (!is.null(cluster_to_celltype)) {
  if (!(gc_cluster_col %in% colnames(obj_gr@meta.data))) {
    stop("GC cluster column not found: ", gc_cluster_col, call. = FALSE)
  }
  obj_gr[[gc_celltype_col]] <- NA_character_
  clust_vec <- as.character(obj_gr@meta.data[[gc_cluster_col]])
  for (ct in names(cluster_to_celltype)) {
    cls <- as.character(unlist(cluster_to_celltype[[ct]]))
    obj_gr@meta.data[[gc_celltype_col]][clust_vec %in% cls] <- ct
  }
}

if (!(gc_celltype_col %in% colnames(obj_gr@meta.data))) {
  stop("GC celltype column not found: ", gc_celltype_col, call. = FALSE)
}
obj_gr@meta.data[[gc_celltype_col]] <- factor(
  as.character(obj_gr@meta.data[[gc_celltype_col]]),
  levels = gc_celltype_levels
)

if (!(gc_sample_col %in% colnames(obj_gr@meta.data))) {
  stop("GC sample column not found: ", gc_sample_col, call. = FALSE)
}
obj_gr$labels <- obj_gr@meta.data[[gc_celltype_col]]
obj_gr$T <- obj_gr@meta.data[[gc_sample_col]]

if (!is.null(gc_exclude_celltypes) && length(gc_exclude_celltypes) > 0) {
  obj_gr <- subset(obj_gr, subset = !(labels %in% gc_exclude_celltypes))
}

if (!(oo_stage_col %in% colnames(obj_oo@meta.data))) {
  stop("Oocyte stage column not found: ", oo_stage_col, call. = FALSE)
}
obj_oo$labels <- obj_oo@meta.data[[oo_stage_col]]

if (!is.null(oocyte_exclude_stages) && length(oocyte_exclude_stages) > 0) {
  obj_oo <- subset(obj_oo, subset = !(labels %in% oocyte_exclude_stages))
}

obj_merge <- merge(obj_gr, obj_oo)
obj_merge <- NormalizeData(obj_merge)

CellChatDB <- CellChatDB.mouse
tag <- "triMean"
if (mode == "cellcell_contact") {
  CellChatDB <- subsetDB(CellChatDB.mouse, search = "Cell-Cell Contact")
  tag <- "triMean_CCC"
}

process_cellchat <- function(cc_obj) {
  cc_obj_processed <- cc_obj %>%
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    smoothData(adj = PPI.mouse) %>%
    computeCommunProb(type = compute_type, population.size = population_size, seed.use = seed_use) %>%
    filterCommunication(min.cells = min_cells) %>%
    computeCommunProbPathway() %>%
    aggregateNet()

  df.net <- subsetCommunication(cc_obj_processed)

  db.inter <- cc_obj_processed@DB$interaction
  db.inter.small <- unique(db.inter[, c("interaction_name", "annotation")])
  df.net.anno <- merge(df.net, db.inter.small, by = "interaction_name", all.x = TRUE)

  list(cellchat = cc_obj_processed, communication = df.net.anno)
}

group.new <- factor(group_new_levels, levels = group_new_levels)

results <- list()

for (st in stage_subsets) {
  st_name <- st$name
  keep_labels <- st$keep_labels
  if (is.null(st_name) || is.null(keep_labels)) stop("stage_subsets must have name + keep_labels", call. = FALSE)

  obj_sub <- subset(obj_merge, subset = labels %in% keep_labels)
  obj_sub$labels <- factor(as.character(obj_sub$labels), levels = keep_labels)

  data.input <- GetAssayData(obj_sub, slot = "data")
  Idents(obj_sub) <- "labels"
  idents <- droplevels(Idents(obj_sub))

  cat("[Stage]", st_name, "cells:", ncol(obj_sub), "groups:", paste(levels(idents), collapse = ","), "\n")

  cc <- createCellChat(
    object = data.input,
    meta = data.frame(group = idents, row.names = names(idents)),
    group.by = "group"
  )
  cc@DB <- CellChatDB

  processed <- process_cellchat(cc)
  cc_lift <- liftCellChat(processed$cellchat, group.new)

  out_rds <- file.path(out_dir, paste0("cellchat_", st_name, "_", tag, ".rds"))
  saveRDS(cc_lift, out_rds)
  cat("[Info] wrote:", out_rds, "\n")

  if (isTRUE(export_comm_tables)) {
    out_csv <- file.path(out_dir, paste0("comm_", st_name, "_", tag, ".csv"))
    write.csv(processed$communication, out_csv, quote = FALSE, row.names = FALSE)
    cat("[Info] wrote:", out_csv, "\n")
  }

  results[[st_name]] <- list(cellchat = cc_lift, communication = processed$communication)
}

out_all <- file.path(out_dir, paste0(prefix, ".cellchat_all_", tag, ".rds"))
saveRDS(results, out_all)
cat("[Info] wrote:", out_all, "\n")

cat("[OK] Done.\n")
sink()
