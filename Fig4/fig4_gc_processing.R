#!/usr/bin/env Rscript
# Fig4: updated Granulosa cell (GC) processing
#
# This version follows the user's newer GC workflow more closely:
#   1) Optionally relabel one or more source clusters in obj_total before subsetting GC
#   2) Subset GC-related celltypes
#   3) Cell-cycle scoring
#   4) Round 1: Scale (optionally regress), PCA, Harmony, neighbors/clusters/UMAP
#      - supports two resolutions run sequentially on the same harmony embedding
#   5) Remove specified low-quality clusters
#   6) Intermediate re-clustering on existing harmony embedding and optional stage-1/2 celltype labels
#   7) Final rerun with cell-cycle regression, Harmony, neighbors/clusters/UMAP, final labels
#   8) Optional removal of one or more final clusters, then re-cluster again while preserving celltype labels
#
# Output object stores refined labels in meta.data$celltype so downstream ssGSEA / c2l_sc can use it directly.

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
  Rscript fig4_gc_processing.R --config Fig4/configs/fig4_combined.yaml
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
    feats <- setdiff(feats, as.character(unlist(manual_exclude)))
  }
  feats
}

apply_cluster_to_celltype <- function(obj, mapping, target_col = "celltype", levels_vec = NULL, init_col = NULL) {
  if (is.null(mapping) || length(mapping) == 0) return(obj)
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("seurat_clusters missing in object.", call. = FALSE)

  md <- obj@meta.data
  if (!is.null(init_col) && init_col %in% colnames(md)) {
    new_vals <- as.character(md[[init_col]])
  } else if (target_col %in% colnames(md)) {
    new_vals <- as.character(md[[target_col]])
  } else {
    new_vals <- rep(NA_character_, nrow(md))
  }

  clus <- as.character(md$seurat_clusters)
  for (ct in names(mapping)) {
    cls <- as.character(unlist(mapping[[ct]]))
    new_vals[clus %in% cls] <- ct
  }

  if (is.null(levels_vec) || length(levels_vec) == 0) {
    levels_vec <- names(mapping)
  }
  md[[target_col]] <- factor(new_vals, levels = as.character(levels_vec))
  obj@meta.data <- md
  obj
}

subset_excluding_clusters <- function(obj, exclude_clusters) {
  exclude_clusters <- as.character(unlist(exclude_clusters))
  if (length(exclude_clusters) == 0) return(obj)
  keep_cells <- rownames(obj@meta.data)[!(as.character(obj@meta.data$seurat_clusters) %in% exclude_clusters)]
  subset(obj, cells = keep_cells)
}

run_scale_pca_harmony <- function(obj, features, regress_vars = NULL, pca_npcs = 50, harmony_group = "orig.ident",
                                  theta = NULL, lambda = NULL, sigma = NULL, verbose = TRUE) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  VariableFeatures(obj) <- features

  if (!is.null(regress_vars) && length(regress_vars) > 0) {
    obj <- ScaleData(obj, vars.to.regress = regress_vars, verbose = FALSE)
  } else {
    obj <- ScaleData(obj, verbose = FALSE)
  }

  obj <- RunPCA(obj, features = features, npcs = pca_npcs, verbose = FALSE)

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
  obj
}

run_neighbors_clusters_umap <- function(obj, dims = 1:30, resolution = 1.0, reduction = "harmony") {
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, reduction = reduction, resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = reduction, dims = dims, verbose = FALSE)
  obj
}

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
if (!"celltype" %in% colnames(obj@meta.data)) stop("meta.data missing `celltype` column.", call. = FALSE)

# Optional pre-subset relabel in full object (e.g. cluster 9 -> GC_Column)
pre_map <- cfg$pre_subset_cluster_to_celltype
if (!is.null(pre_map) && length(pre_map) > 0) {
  obj <- apply_cluster_to_celltype(obj, pre_map, target_col = "celltype", init_col = "celltype")
  message("[fig4-gc] applied pre_subset_cluster_to_celltype to full object.")
}

subset_celltypes <- as.character(unlist(cfg$subset_celltypes %||% c("GC_Antral_Mural", "GC_Preantral", "GC_Atretic", "GC_Column")))
obj_gr <- subset(obj, subset = celltype %in% subset_celltypes)
message("[fig4-gc] subset celltypes: ", paste(subset_celltypes, collapse = ", "))
message("[fig4-gc] subset cells: ", ncol(obj_gr))

if (!"orig.ident" %in% colnames(obj_gr@meta.data)) obj_gr$orig.ident <- Idents(obj_gr)
obj_gr$orig.ident <- factor(obj_gr$orig.ident)
Idents(obj_gr) <- "orig.ident"

# features + cell cycle
nfeatures1 <- as.integer(cfg$nfeatures_round1 %||% 2000)
exclude_regex <- as.character(cfg$exclude_regex %||% "^(Rp|mt)")
obj.list1 <- prep_split_list(obj_gr, nfeatures = nfeatures1)
features1 <- select_features_excluding(obj_gr, obj.list1, nfeatures = nfeatures1,
                                       exclude_regex = exclude_regex,
                                       manual_exclude = cfg$manual_exclude_genes_round1)
message("[fig4-gc] features round1: ", length(features1))

cc <- cell_cycle_lists_mouse()
cc$s.genes <- intersect(cc$s.genes, rownames(obj_gr))
cc$g2m.genes <- intersect(cc$g2m.genes, rownames(obj_gr))
obj_gr <- CellCycleScoring(obj_gr, s.features = cc$s.genes, g2m.features = cc$g2m.genes, set.ident = FALSE)
message("[fig4-gc] cell-cycle genes used: S=", length(cc$s.genes), ", G2M=", length(cc$g2m.genes))

# Round 1: match user's first pipeline (no regression by default)
p1 <- cfg$pipeline_round1 %||% list()
obj_r1 <- run_scale_pca_harmony(
  obj_gr,
  features = features1,
  regress_vars = if (!is.null(p1$vars_to_regress)) as.character(unlist(p1$vars_to_regress)) else NULL,
  pca_npcs = as.integer(p1$pca_npcs %||% 50),
  harmony_group = as.character(p1$harmony_group %||% "orig.ident"),
  theta = if (!is.null(p1$theta)) as.numeric(p1$theta) else NULL,
  lambda = if (!is.null(p1$lambda)) as.numeric(p1$lambda) else NULL,
  sigma = if (!is.null(p1$sigma)) as.numeric(p1$sigma) else NULL,
  verbose = TRUE
)

dims_r1 <- seq_len(as.integer(p1$dims %||% 40))
res1a <- as.numeric(p1$cluster_resolution_stage1 %||% 0.5)
res1b <- as.numeric(p1$cluster_resolution_stage2 %||% 1.5)
obj_r1 <- run_neighbors_clusters_umap(obj_r1, dims = dims_r1, resolution = res1a)
message("[fig4-gc] round1 stage1 resolution: ", res1a)
out_round1_stage1_rds <- as.character(cfg$out_round1_stage1_rds %||% file.path(out_dir, paste0(prefix, ".round1_stage1.rds")))
saveRDS(obj_r1, out_round1_stage1_rds)
obj_r1 <- run_neighbors_clusters_umap(obj_r1, dims = dims_r1, resolution = res1b)
message("[fig4-gc] round1 stage2 resolution: ", res1b)
out_round1_rds <- as.character(cfg$out_round1_rds %||% file.path(out_dir, paste0(prefix, ".round1.rds")))
saveRDS(obj_r1, out_round1_rds)

p_umap1 <- DimPlot(obj_r1, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("GC round1 clusters")
p_orig1 <- DimPlot(obj_r1, reduction = "umap", group.by = "orig.ident", label = FALSE) + ggtitle("GC round1 orig.ident")
ggsave(file.path(out_dir, paste0(prefix, ".round1.umap.pdf")), p_orig1 + p_umap1, width = 12, height = 5)

exclude1 <- as.character(unlist(cfg$exclude_clusters_round1 %||% c("6","7","15","16")))
obj_r1f <- subset_excluding_clusters(obj_r1, exclude1)
message("[fig4-gc] exclude round1 clusters: ", paste(exclude1, collapse = ", "))
message("[fig4-gc] cells after round1 filter: ", ncol(obj_r1f))

# Stage 2: recluster on existing harmony, optional stage labels
p2 <- cfg$pipeline_stage2 %||% list()
dims_stage2 <- seq_len(as.integer(p2$dims %||% 40))
res_stage2 <- as.numeric(p2$cluster_resolution %||% 1.0)
obj_stage2 <- run_neighbors_clusters_umap(obj_r1f, dims = dims_stage2, resolution = res_stage2)
message("[fig4-gc] stage2 recluster resolution: ", res_stage2)

if (!is.null(cfg$stage2_label_levels) && !is.null(cfg$stage2_cluster_to_celltype)) {
  obj_stage2 <- apply_cluster_to_celltype(
    obj_stage2,
    mapping = cfg$stage2_cluster_to_celltype,
    target_col = as.character(cfg$celltype_col %||% "celltype"),
    levels_vec = as.character(unlist(cfg$stage2_label_levels))
  )
  message("[fig4-gc] applied stage2 cluster->celltype mapping.")
}

# Final rerun with cell-cycle regression
p3 <- cfg$pipeline_final %||% list()
obj_final <- run_scale_pca_harmony(
  obj_stage2,
  features = features1,
  regress_vars = as.character(unlist(p3$vars_to_regress %||% c("S.Score", "G2M.Score"))),
  pca_npcs = as.integer(p3$pca_npcs %||% 50),
  harmony_group = as.character(p3$harmony_group %||% "orig.ident"),
  theta = if (!is.null(p3$theta)) as.numeric(p3$theta) else NULL,
  lambda = if (!is.null(p3$lambda)) as.numeric(p3$lambda) else NULL,
  sigma = if (!is.null(p3$sigma)) as.numeric(p3$sigma) else NULL,
  verbose = TRUE
)

dims_final <- seq_len(as.integer(p3$dims %||% 30))
res_final <- as.numeric(p3$cluster_resolution %||% 1.0)
obj_final <- run_neighbors_clusters_umap(obj_final, dims = dims_final, resolution = res_final)
message("[fig4-gc] final resolution before extra exclusion: ", res_final)

celltype_col <- as.character(cfg$celltype_col %||% "celltype")
if (!is.null(cfg$final_label_levels) && !is.null(cfg$final_cluster_to_celltype)) {
  obj_final <- apply_cluster_to_celltype(
    obj_final,
    mapping = cfg$final_cluster_to_celltype,
    target_col = celltype_col,
    levels_vec = as.character(unlist(cfg$final_label_levels))
  )
  message("[fig4-gc] applied final cluster->celltype mapping.")
}

# optional extra exclusion after final mapping (preserve celltype labels)
exclude_final <- as.character(unlist(cfg$exclude_clusters_after_final_map %||% character(0)))
if (length(exclude_final) > 0) {
  obj_final <- subset_excluding_clusters(obj_final, exclude_final)
  message("[fig4-gc] excluded clusters after final mapping: ", paste(exclude_final, collapse = ", "))
  obj_final <- run_neighbors_clusters_umap(obj_final, dims = dims_final, resolution = res_final)
  message("[fig4-gc] reran neighbors/clusters/UMAP after post-final exclusion.")
}

# save outputs
out_stage2_rds <- as.character(cfg$out_stage2_rds %||% file.path(out_dir, paste0(prefix, ".stage2.rds")))
saveRDS(obj_stage2, out_stage2_rds)
out_final_rds <- as.character(cfg$out_final_rds %||% file.path(out_dir, paste0(prefix, ".final.rds")))
saveRDS(obj_final, out_final_rds)
message("[fig4-gc] saved final: ", out_final_rds)

plot_group <- if (celltype_col %in% colnames(obj_final@meta.data)) celltype_col else "seurat_clusters"
p_final_ct <- DimPlot(obj_final, reduction = "umap", group.by = plot_group, label = TRUE, repel = TRUE) + ggtitle("Final GC labels")
p_final_cl <- DimPlot(obj_final, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("Final GC clusters")
p_final_orig <- DimPlot(obj_final, reduction = "umap", group.by = "orig.ident") + ggtitle("Final orig.ident")
ggsave(file.path(out_dir, paste0(prefix, ".final.umap.pdf")), (p_final_orig | p_final_cl) / p_final_ct, width = 14, height = 9)

count_df <- obj_final@meta.data %>%
  dplyr::count(orig.ident,
               sample = if ("sample" %in% colnames(obj_final@meta.data)) sample else NA,
               seurat_clusters,
               label = if (celltype_col %in% colnames(obj_final@meta.data)) .data[[celltype_col]] else NA,
               name = "n_cells")
write.csv(count_df, file.path(out_dir, paste0(prefix, ".cell_counts.csv")), row.names = FALSE)

params_used <- cfg
params_used$resolved <- list(
  input_rds = in_rds,
  out_round1_stage1_rds = out_round1_stage1_rds,
  out_round1_rds = out_round1_rds,
  out_stage2_rds = out_stage2_rds,
  out_final_rds = out_final_rds,
  log_path = log_path,
  features_round1_n = length(features1)
)
yaml::write_yaml(params_used, file.path(out_dir, "params_used_fig4_gc_processing.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig4_gc_processing.txt"), capture.output(sessionInfo()))

message("[fig4-gc] done: ", Sys.time())

sink()
sink(type = "message")
close(log_file)
