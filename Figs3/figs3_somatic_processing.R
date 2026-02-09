#!/usr/bin/env Rscript
# Figs3: Somatic-cell integration / clustering / annotation (standalone figure)
#
# This script merges multiple QC-processed Seurat objects across timepoints,
# performs variable-feature selection (excluding ribosomal/mitochondrial genes),
# regresses out cell-cycle scores, runs PCA + Harmony integration, clustering, UMAP,
# removes specified low-quality clusters in two rounds, and annotates cell types
# from a cluster->celltype mapping.
#
# Usage:
#   Rscript figs3_somatic_processing.R --config Figs3/configs/figs3_somatic_processing.yaml
#
suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(harmony)
  library(ggplot2)
  library(patchwork)
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
"figs3_somatic_processing.R

Required:
  --config <yaml>     YAML config path

Example:
  Rscript figs3_somatic_processing.R --config Figs3/configs/figs3_somatic_processing.yaml
", sep = "")
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

%||% <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

# ---------- Helpers ----------
load_group_object <- function(group_cfg) {
  # group_cfg:
  #   name: "PD14"
  #   sample_label: "P14"
  #   items: [{id:"PD14-1", rds:"..."}, ...]
  if (is.null(group_cfg$items) || length(group_cfg$items) == 0) stop("Each group must have non-empty items.", call. = FALSE)
  group_name <- as.character(group_cfg$name)
  if (!nzchar(group_name)) stop("Group name missing.", call. = FALSE)
  sample_label <- as.character(group_cfg$sample_label)
  if (!nzchar(sample_label)) sample_label <- group_name

  objs <- list()
  ids <- c()
  for (it in group_cfg$items) {
    id <- as.character(it$id)
    rds <- as.character(it$rds)
    if (!nzchar(id) || !nzchar(rds)) stop("Each item needs id and rds.", call. = FALSE)
    if (!file.exists(rds)) stop("Missing RDS: ", rds, call. = FALSE)
    message("[figs3] read: ", id, " <- ", rds)
    o <- readRDS(rds)
    if (!inherits(o, "Seurat")) stop("Expected Seurat object: ", rds, call. = FALSE)
    objs[[id]] <- o
    ids <- c(ids, id)
  }

  # Merge within group
  if (length(objs) == 1) {
    obj <- objs[[1]]
    # Add sample label
    obj$sample <- sample_label
    return(obj)
  }

  first_id <- ids[[1]]
  obj <- objs[[first_id]]
  rest <- objs[ids[-1]]
  obj <- merge(obj, y = rest, add.cell.ids = ids)
  obj$sample <- sample_label
  return(obj)
}

prep_split_list <- function(obj, nfeatures = 2000) {
  # Split by orig.ident (Idents should be set before calling)
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

cell_cycle_lists_mouse <- function() {
  s.genes <- c(
    "Mcm5","Pcna","Tyms","Fen1","Mcm7","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl",
    "Prim1","Uhrf1","Cenpu","Hells","Rfc2","Polr1b","Nasp","Rad51ap1","Gmnn","Wdr76",
    "Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1",
    "Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Mrpl36","E2f8"
  )
  g2m.genes <- c(
    "Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67",
    "Tmpo","Cenpf","Tacc3","Pimreg","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
    "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2",
    "Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf",
    "Nek2","G2e3","Gas2l3","Cbx5","Cenpa"
  )
  list(s.genes = s.genes, g2m.genes = g2m.genes)
}

run_harmony_pipeline <- function(
  obj,
  features,
  s.genes,
  g2m.genes,
  harmony_group = "orig.ident",
  pca_dims = 30,
  harmony_theta = 1.5,
  harmony_lambda = 1.5,
  harmony_sigma = 0.05,
  neighbors_dims = 1:30,
  cluster_resolution = 1.0,
  umap_dims = 1:30,
  verbose = TRUE
) {
  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(obj, verbose = FALSE)
  VariableFeatures(obj) <- features

  # cell cycle scoring
  s.genes <- intersect(s.genes, rownames(obj))
  g2m.genes <- intersect(g2m.genes, rownames(obj))
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

  # scale + regress
  obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)

  # PCA
  obj <- RunPCA(obj, features = features, verbose = FALSE)

  # Harmony
  obj <- RunHarmony(
    object = obj,
    group.by.vars = harmony_group,
    reduction = "pca",
    reduction.save = "harmony",
    theta = harmony_theta,
    lambda = harmony_lambda,
    sigma = harmony_sigma,
    verbose = verbose
  )

  # graph + clustering + UMAP
  obj <- FindNeighbors(obj, reduction = "harmony", dims = neighbors_dims, verbose = FALSE)
  obj <- FindClusters(obj, reduction = "harmony", resolution = cluster_resolution, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = umap_dims, verbose = FALSE)

  obj
}

apply_celltype_map <- function(obj, celltype_map) {
  # celltype_map is a named list: {Celltype: [cluster_ids]}
  if (is.null(celltype_map) || length(celltype_map) == 0) return(obj)
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("seurat_clusters missing in object.", call. = FALSE)

  obj$celltype <- NA_character_
  for (ct in names(celltype_map)) {
    cls <- as.character(unlist(celltype_map[[ct]]))
    obj$celltype[obj$seurat_clusters %in% cls] <- ct
  }
  obj$celltype <- factor(obj$celltype)
  obj
}

# ---------- Main ----------
args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

cfg <- yaml::read_yaml(args$config)

out_dir <- cfg$out_dir
if (is.null(out_dir) || !nzchar(as.character(out_dir))) out_dir <- "results/Figs3/somatic_processing"
out_dir <- as.character(out_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- cfg$prefix
if (is.null(prefix) || !nzchar(as.character(prefix))) prefix <- "somatic"
prefix <- as.character(prefix)

log_path <- file.path(out_dir, paste0(prefix, ".processing.log"))
log_file <- file(log_path, open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")

message("[figs3] out_dir: ", out_dir)
message("[figs3] start: ", Sys.time())

# ---- Load inputs ----
if (is.null(cfg$groups) || length(cfg$groups) == 0) {
  stop("Config must contain `groups` (list).", call. = FALSE)
}

group_objs <- list()
for (g in cfg$groups) {
  obj_g <- load_group_object(g)
  group_name <- as.character(g$name)
  group_objs[[group_name]] <- obj_g
}

# Merge across groups
group_names <- names(group_objs)
if (length(group_names) == 1) {
  obj_merge <- group_objs[[1]]
} else {
  obj_merge <- group_objs[[group_names[[1]]]]
  rest <- group_objs[group_names[-1]]
  obj_merge <- merge(obj_merge, y = rest)
}
message("[figs3] merged cells: ", ncol(obj_merge))

# Ensure factors
obj_merge$orig.ident <- factor(obj_merge$orig.ident)
obj_merge$sample <- factor(obj_merge$sample)

# ---- Round 1 integration ----
Idents(obj_merge) <- "orig.ident"

nfeatures_round1 <- cfg$nfeatures_round1
if (is.null(nfeatures_round1)) nfeatures_round1 <- 2000
nfeatures_round1 <- as.integer(nfeatures_round1)

exclude_regex <- cfg$exclude_regex
if (is.null(exclude_regex) || !nzchar(as.character(exclude_regex))) exclude_regex <- "^(Rp|mt)"
exclude_regex <- as.character(exclude_regex)

manual_exclude <- cfg$manual_exclude_genes_round1

obj.list <- prep_split_list(obj_merge, nfeatures = nfeatures_round1)
features1 <- select_features_excluding(obj_merge, obj.list, nfeatures = nfeatures_round1,
                                      exclude_regex = exclude_regex, manual_exclude = manual_exclude)

cc <- cell_cycle_lists_mouse()
p1 <- cfg$pipeline_round1
if (is.null(p1)) p1 <- list()

obj_merge_proc <- run_harmony_pipeline(
  obj = obj_merge,
  features = features1,
  s.genes = cc$s.genes,
  g2m.genes = cc$g2m.genes,
  harmony_group = as.character(p1$harmony_group %||% "orig.ident"),
  pca_dims = as.integer(p1$pca_dims %||% 30),
  harmony_theta = as.numeric(p1$theta %||% 1.5),
  harmony_lambda = as.numeric(p1$lambda %||% 1.5),
  harmony_sigma = as.numeric(p1$sigma %||% 0.05),
  neighbors_dims = 1:as.integer(p1$neighbors_dim_max %||% 30),
  cluster_resolution = as.numeric(p1$cluster_resolution %||% 1.0),
  umap_dims = 1:as.integer(p1$umap_dim_max %||% 30),
  verbose = TRUE
)

# Save round1 object (optional)
out_round1_rds <- cfg$out_round1_rds
if (is.null(out_round1_rds) || !nzchar(as.character(out_round1_rds))) {
  out_round1_rds <- file.path(out_dir, paste0(prefix, ".round1.rds"))
}
saveRDS(obj_merge_proc, out_round1_rds)
message("[figs3] saved round1: ", out_round1_rds)

# Round1 plots
p_umap1_sample <- DimPlot(obj_merge_proc, reduction = "umap", group.by = "sample", label = FALSE) + ggtitle("UMAP by sample (round1)")
p_umap1_orig   <- DimPlot(obj_merge_proc, reduction = "umap", group.by = "orig.ident", label = FALSE) + ggtitle("UMAP by orig.ident (round1)")
p_umap1_clust  <- DimPlot(obj_merge_proc, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("UMAP by clusters (round1)")
ggsave(file.path(out_dir, paste0(prefix, ".round1.umap.pdf")), p_umap1_sample + p_umap1_orig + p_umap1_clust, width = 14, height = 5)

# ---- Round 1 exclusion ----
exclude1 <- cfg$exclude_clusters_round1
if (is.null(exclude1)) exclude1 <- c("3", "14")
exclude1 <- as.character(unlist(exclude1))

obj_subset <- subset(obj_merge_proc, subset = !(seurat_clusters %in% exclude1))
message("[figs3] after exclude round1 clusters {", paste(exclude1, collapse = ","), "} cells: ", ncol(obj_subset))

# ---- Round 2 integration ----
Idents(obj_subset) <- "orig.ident"

nfeatures_round2 <- cfg$nfeatures_round2
if (is.null(nfeatures_round2)) nfeatures_round2 <- 2000
nfeatures_round2 <- as.integer(nfeatures_round2)

manual_exclude2 <- cfg$manual_exclude_genes_round2

obj.list2 <- prep_split_list(obj_subset, nfeatures = nfeatures_round2)
features2 <- select_features_excluding(obj_subset, obj.list2, nfeatures = nfeatures_round2,
                                      exclude_regex = exclude_regex, manual_exclude = manual_exclude2)

p2 <- cfg$pipeline_round2
if (is.null(p2)) p2 <- list()

obj_subset_proc <- run_harmony_pipeline(
  obj = obj_subset,
  features = features2,
  s.genes = cc$s.genes,
  g2m.genes = cc$g2m.genes,
  harmony_group = as.character(p2$harmony_group %||% "orig.ident"),
  pca_dims = as.integer(p2$pca_dims %||% 40),
  harmony_theta = as.numeric(p2$theta %||% 1.5),
  harmony_lambda = as.numeric(p2$lambda %||% 1.5),
  harmony_sigma = as.numeric(p2$sigma %||% 0.03),
  neighbors_dims = 1:as.integer(p2$neighbors_dim_max %||% 40),
  cluster_resolution = as.numeric(p2$cluster_resolution %||% 1.0),
  umap_dims = 1:as.integer(p2$umap_dim_max %||% 40),
  verbose = TRUE
)

# Round2 plots
p_umap2_sample <- DimPlot(obj_subset_proc, reduction = "umap", group.by = "sample", label = FALSE) + ggtitle("UMAP by sample (round22)")
p_umap2_orig   <- DimPlot(obj_subset_proc, reduction = "umap", group.by = "orig.ident", label = FALSE) + ggtitle("UMAP by orig.ident (round2)")
p_umap2_clust  <- DimPlot(obj_subset_proc, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("UMAP by clusters (round2)")
ggsave(file.path(out_dir, paste0(prefix, ".round2.umap.pdf")), p_umap2_sample + p_umap2_orig + p_umap2_clust, width = 14, height = 5)

# ---- Round 2 exclusion ----
exclude2 <- cfg$exclude_clusters_round2
if (is.null(exclude2)) exclude2 <- c("25", "12")
exclude2 <- as.character(unlist(exclude2))

obj_final <- subset(obj_subset_proc, subset = !(seurat_clusters %in% exclude2))
message("[figs3] after exclude round2 clusters {", paste(exclude2, collapse = ","), "} cells: ", ncol(obj_final))

# ---- Cell type annotation ----
celltype_map <- cfg$celltype_map
obj_final <- apply_celltype_map(obj_final, celltype_map)

# Save final object
out_final_rds <- cfg$out_final_rds
if (is.null(out_final_rds) || !nzchar(as.character(out_final_rds))) {
  out_final_rds <- file.path(out_dir, paste0(prefix, ".final.rds"))
}
saveRDS(obj_final, out_final_rds)
message("[figs3] saved final: ", out_final_rds)

# Final plots
p_umap_ct <- DimPlot(obj_final, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("UMAP by celltype (final)")
p_umap_sample <- DimPlot(obj_final, reduction = "umap", group.by = "sample", label = FALSE) + ggtitle("UMAP by sample (final)")
p_umap_orig <- DimPlot(obj_final, reduction = "umap", group.by = "orig.ident", label = FALSE) + ggtitle("UMAP by orig.ident (final)")
ggsave(file.path(out_dir, paste0(prefix, ".final.umap.pdf")), p_umap_sample + p_umap_orig + p_umap_ct, width = 16, height = 5)

# Counts table
ct_tbl <- obj_final@meta.data %>%
  dplyr::count(sample, celltype, name = "n_cells") %>%
  dplyr::arrange(sample, dplyr::desc(n_cells))
write.csv(ct_tbl, file.path(out_dir, paste0(prefix, ".cell_counts_by_sample_celltype.csv")), row.names = FALSE)

# Provenance
params_used <- cfg
params_used$out_dir <- out_dir
params_used$out_round1_rds <- out_round1_rds
params_used$out_final_rds <- out_final_rds
params_used$log_path <- log_path
yaml::write_yaml(params_used, file.path(out_dir, "params_used_figs3_somatic_processing.yaml"))
write_text(file.path(out_dir, "sessionInfo_figs3_somatic_processing.txt"), capture.output(sessionInfo()))

message("[figs3] done: ", Sys.time())

sink()
sink(type = "message")
close(log_file)
