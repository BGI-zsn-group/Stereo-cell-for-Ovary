#!/usr/bin/env Rscript
# Fig4: Granulosa cell (GC) processing (standalone module for Nature-style reproducibility)
#
# This script:
#   1) Loads an annotated Seurat object (obj_all_withanno.rds)
#   2) Subsets GC-related celltypes (default: GC_Antral_Mural, GC_Preantral, GC_Column)
#   3) Harmony integration + clustering + UMAP (Round 1)
#   4) Removes low-quality clusters (Round 1)
#   5) Re-runs Harmony pipeline (Round 2) and removes additional low-quality cluster(s)
#   6) Final clustering/UMAP and assigns refined GC subtypes by cluster->celltype mapping
#
# Usage:
#   Rscript fig4_gc_processing.R --config Fig4/configs/fig4_gc_processing.yaml
#
suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

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
"fig4_gc_processing.R

Required:
  --config <yaml>     YAML config path

Example:
  Rscript fig4_gc_processing.R --config Fig4/configs/fig4_gc_processing.yaml
", sep = "")
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

cell_cycle_lists_mouse <- function() {
  s.genes <- c(
    "Mcm5","Pcna","Tyms","Fen1","Mcm7","Mcm4","Rrm1","Ung","Gins2","Mcm6",
    "Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Polr1b","Nasp",
    "Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2",
    "Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2",
    "Usp1","Clspn","Pola1","Chaf1b","Mrpl36","E2f8"
  )
  g2m.genes <- c(
    "Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2",
    "Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Pimreg","Smc4","Ccnb2",
    "Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1",
    "Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk","Cdc25c","Kif2c",
    "Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr",
    "Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3",
    "Gas2l3","Cbx5","Cenpa"
  )
  list(s.genes = s.genes, g2m.genes = g2m.genes)
}

prep_split_list <- function(obj, nfeatures = 2000) {
  obj.list <- SplitObject(obj)
  obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    x
  })
  obj.list
}

select_features_excluding <- function(obj, obj.list, nfeatures = 2000, exclude_regex = "^(Rp|mt)", manual_exclude = NULL) {
  feats <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = nfeatures)
  fl_genes <- grep(exclude_regex, rownames(obj), value = TRUE)
  feats <- setdiff(feats, fl_genes)
  if (!is.null(manual_exclude) && length(manual_exclude) > 0) {
    manual_exclude <- as.character(unlist(manual_exclude))
    feats <- setdiff(feats, manual_exclude)
  }
  feats
}

run_harmony_pipeline <- function(
  obj,
  features,
  regress_vars = c("S.Score", "G2M.Score"),
  pca_npcs = 50,
  harmony_group = "orig.ident",
  theta = NULL,
  lambda = NULL,
  sigma = NULL,
  neighbors_dim_max = 30,
  cluster_resolution = 1.0,
  umap_dim_max = 30,
  verbose = TRUE
) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  VariableFeatures(obj) <- features

  obj <- ScaleData(obj, vars.to.regress = regress_vars, verbose = FALSE)
  obj <- RunPCA(obj, features = features, npcs = pca_npcs, verbose = FALSE)

  # Harmony params are optional; if NULL, Harmony defaults are used.
  args <- list(
    object = obj,
    group.by.vars = harmony_group,
    reduction = "pca",
    reduction.save = "harmony",
    verbose = verbose
  )
  if (!is.null(theta))  args$theta  <- theta
  if (!is.null(lambda)) args$lambda <- lambda
  if (!is.null(sigma))  args$sigma  <- sigma
  obj <- do.call(RunHarmony, args)

  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:neighbors_dim_max, verbose = FALSE)
  obj <- FindClusters(obj, reduction = "harmony", resolution = cluster_resolution, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:umap_dim_max, verbose = FALSE)
  obj
}

apply_celltype_map <- function(obj, celltype_map) {
  if (is.null(celltype_map) || length(celltype_map) == 0) return(obj)
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("seurat_clusters missing in object.", call. = FALSE)

  obj$celltype_refined <- NA_character_
  for (ct in names(celltype_map)) {
    cls <- as.character(unlist(celltype_map[[ct]]))
    obj$celltype_refined[obj$seurat_clusters %in% cls] <- ct
  }
  obj$celltype_refined <- factor(obj$celltype_refined, levels = names(celltype_map))
  obj
}

# ---------------- Main ----------------
args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

cfg <- yaml::read_yaml(args$config)

out_dir <- as.character(cfg$out_dir %||% "results/Fig4/gc_processing")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- as.character(cfg$prefix %||% "gc")
log_path <- file.path(out_dir, paste0(prefix, ".log"))

log_file <- file(log_path, open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

message("[fig4-gc] start: ", Sys.time())
message("[fig4-gc] out_dir: ", out_dir)

in_rds <- as.character(cfg$input_rds %||% "obj_all_withanno.rds")
if (!file.exists(in_rds)) stop("Missing input_rds: ", in_rds, call. = FALSE)

obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Expected Seurat object in input_rds.", call. = FALSE)

subset_celltypes <- as.character(unlist(cfg$subset_celltypes %||% c("GC_Antral_Mural", "GC_Preantral", "GC_Column")))
if (!"celltype" %in% colnames(obj@meta.data)) stop("meta.data missing `celltype` column.", call. = FALSE)

obj_gr <- subset(obj, subset = celltype %in% subset_celltypes)
message("[fig4-gc] subset celltypes: ", paste(subset_celltypes, collapse = ", "))
message("[fig4-gc] subset cells: ", ncol(obj_gr))

# Make sure identities exist
if (!"orig.ident" %in% colnames(obj_gr@meta.data)) obj_gr$orig.ident <- Idents(obj_gr)
obj_gr$orig.ident <- factor(obj_gr$orig.ident)
Idents(obj_gr) <- "orig.ident"

# ---- Feature selection (Round 1) ----
nfeatures1 <- as.integer(cfg$nfeatures_round1 %||% 2000)
exclude_regex <- as.character(cfg$exclude_regex %||% "^(Rp|mt)")
manual_exclude1 <- cfg$manual_exclude_genes_round1

obj.list1 <- prep_split_list(obj_gr, nfeatures = nfeatures1)
features1 <- select_features_excluding(obj_gr, obj.list1, nfeatures = nfeatures1,
                                       exclude_regex = exclude_regex, manual_exclude = manual_exclude1)

# Cell cycle scoring (mouse lists)
cc <- cell_cycle_lists_mouse()
cc$s.genes <- intersect(cc$s.genes, rownames(obj_gr))
cc$g2m.genes <- intersect(cc$g2m.genes, rownames(obj_gr))
obj_gr <- CellCycleScoring(obj_gr, s.features = cc$s.genes, g2m.features = cc$g2m.genes, set.ident = FALSE)

# ---- Round 1 pipeline ----
p1 <- cfg$pipeline_round1 %||% list()
obj_r1 <- run_harmony_pipeline(
  obj = obj_gr,
  features = features1,
  regress_vars = as.character(unlist(p1$vars_to_regress %||% c("S.Score","G2M.Score"))),
  pca_npcs = as.integer(p1$pca_npcs %||% 50),
  harmony_group = as.character(p1$harmony_group %||% "orig.ident"),
  theta  = if (!is.null(p1$theta))  as.numeric(p1$theta)  else NULL,
  lambda = if (!is.null(p1$lambda)) as.numeric(p1$lambda) else NULL,
  sigma  = if (!is.null(p1$sigma))  as.numeric(p1$sigma)  else NULL,
  neighbors_dim_max = as.integer(p1$neighbors_dim_max %||% 40),
  cluster_resolution = as.numeric(p1$cluster_resolution %||% 1.5),
  umap_dim_max = as.integer(p1$umap_dim_max %||% 40),
  verbose = TRUE
)

out_round1_rds <- as.character(cfg$out_round1_rds %||% file.path(out_dir, paste0(prefix, ".round1.rds")))
saveRDS(obj_r1, out_round1_rds)
message("[fig4-gc] saved round1: ", out_round1_rds)

# Round1 plots
p_umap1_clust <- DimPlot(obj_r1, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("Round1 clusters")
p_umap1_orig  <- DimPlot(obj_r1, reduction = "umap", group.by = "orig.ident", label = FALSE) + ggtitle("Round1 orig.ident")
ggsave(file.path(out_dir, paste0(prefix, ".round1.umap.pdf")), p_umap1_orig + p_umap1_clust, width = 12, height = 5)

# ---- Remove low-quality clusters (Round 1) ----
exclude1 <- as.character(unlist(cfg$exclude_clusters_round1 %||% c("6","7","15","16","19")))
obj_r1f <- subset(obj_r1, subset = !(seurat_clusters %in% exclude1))
message("[fig4-gc] exclude round1 clusters: ", paste(exclude1, collapse = ", "))
message("[fig4-gc] cells after round1 filter: ", ncol(obj_r1f))

# ---- Round 2 pipeline ----
# Feature selection (Round 2) can be either:
#   - reuse integration features from split list (recommended) OR
#   - use default FindVariableFeatures() (as in your original second stage)
round2_mode <- as.character(cfg$round2_feature_mode %||% "default")  # "default" or "integration"

p2 <- cfg$pipeline_round2 %||% list()

if (round2_mode == "integration") {
  Idents(obj_r1f) <- "orig.ident"
  nfeatures2 <- as.integer(cfg$nfeatures_round2 %||% 2000)
  manual_exclude2 <- cfg$manual_exclude_genes_round2
  obj.list2 <- prep_split_list(obj_r1f, nfeatures = nfeatures2)
  features2 <- select_features_excluding(obj_r1f, obj.list2, nfeatures = nfeatures2,
                                         exclude_regex = exclude_regex, manual_exclude = manual_exclude2)
} else {
  # default: Normalize + FindVariableFeatures without forcing a precomputed list
  obj_r1f <- NormalizeData(obj_r1f, verbose = FALSE)
  obj_r1f <- FindVariableFeatures(obj_r1f, verbose = FALSE)
  features2 <- VariableFeatures(obj_r1f)
}

obj_r2 <- run_harmony_pipeline(
  obj = obj_r1f,
  features = features2,
  regress_vars = as.character(unlist(p2$vars_to_regress %||% c("S.Score","G2M.Score"))),
  pca_npcs = as.integer(p2$pca_npcs %||% 50),
  harmony_group = as.character(p2$harmony_group %||% "orig.ident"),
  theta  = if (!is.null(p2$theta))  as.numeric(p2$theta)  else NULL,
  lambda = if (!is.null(p2$lambda)) as.numeric(p2$lambda) else NULL,
  sigma  = if (!is.null(p2$sigma))  as.numeric(p2$sigma)  else NULL,
  neighbors_dim_max = as.integer(p2$neighbors_dim_max %||% 30),
  cluster_resolution = as.numeric(p2$cluster_resolution %||% 1.0),
  umap_dim_max = as.integer(p2$umap_dim_max %||% 30),
  verbose = TRUE
)

# Remove additional low-quality cluster(s) after round2
exclude2 <- as.character(unlist(cfg$exclude_clusters_round2 %||% c("8")))
obj_r2f <- subset(obj_r2, subset = !(seurat_clusters %in% exclude2))
message("[fig4-gc] exclude round2 clusters: ", paste(exclude2, collapse = ", "))
message("[fig4-gc] cells after round2 filter: ", ncol(obj_r2f))

# Final clustering pass (as in your script)
p3 <- cfg$pipeline_final %||% list()
obj_final <- obj_r2f
obj_final <- FindNeighbors(obj_final, reduction = "harmony", dims = 1:as.integer(p3$neighbors_dim_max %||% 30), verbose = FALSE)
obj_final <- FindClusters(obj_final, resolution = as.numeric(p3$cluster_resolution %||% 1.5), verbose = FALSE)
obj_final <- RunUMAP(obj_final, reduction = "harmony", dims = 1:as.integer(p3$umap_dim_max %||% 30), verbose = FALSE)

# Assign refined GC subtypes by cluster mapping
celltype_map <- cfg$celltype_map %||% list()
obj_final <- apply_celltype_map(obj_final, celltype_map)

# Save final
out_final_rds <- as.character(cfg$out_final_rds %||% file.path(out_dir, paste0(prefix, ".final.rds")))
saveRDS(obj_final, out_final_rds)
message("[fig4-gc] saved final: ", out_final_rds)

# Final plots
p_umap_final_ct <- if ("celltype_refined" %in% colnames(obj_final@meta.data)) {
  DimPlot(obj_final, reduction = "umap", group.by = "celltype_refined", label = TRUE, repel = TRUE) + ggtitle("Final refined celltypes")
} else {
  DimPlot(obj_final, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("Final clusters")
}
p_umap_final_orig <- DimPlot(obj_final, reduction = "umap", group.by = "orig.ident", label = FALSE) + ggtitle("Final orig.ident")
ggsave(file.path(out_dir, paste0(prefix, ".final.umap.pdf")), p_umap_final_orig + p_umap_final_ct, width = 14, height = 5)

# Cell counts
count_df <- obj_final@meta.data %>%
  dplyr::count(orig.ident, sample = if ("sample" %in% colnames(obj_final@meta.data)) sample else NA,
               celltype_refined = if ("celltype_refined" %in% colnames(obj_final@meta.data)) celltype_refined else NA,
               name = "n_cells")
write.csv(count_df, file.path(out_dir, paste0(prefix, ".cell_counts.csv")), row.names = FALSE)

# Provenance
params_used <- cfg
params_used$computed_features_round1_n <- length(features1)
params_used$computed_features_round2_n <- length(features2)
params_used$outputs <- list(
  out_round1_rds = out_round1_rds,
  out_final_rds = out_final_rds,
  log_path = log_path
)
yaml::write_yaml(params_used, file.path(out_dir, "params_used_fig4_gc_processing.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig4_gc_processing.txt"), capture.output(sessionInfo()))

message("[fig4-gc] done: ", Sys.time())

sink()
sink(type = "message")
close(log_file)
