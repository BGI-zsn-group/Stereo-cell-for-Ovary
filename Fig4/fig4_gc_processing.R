#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    k <- args[[i]]
    if (!startsWith(k, "--")) stop("Unknown arg: ", k, call. = FALSE)
    key <- sub("^--", "", k)
    if (key %in% c("help")) {
      out[[key]] <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) stop("Missing value for ", k, call. = FALSE)
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

help_msg <- function() {
  cat(
"fig4_gc_processing.R\n\nRequired:\n  --config <yaml>\n", sep = "")
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

normalize_map <- function(x) {
  if (is.null(x) || length(x) == 0) return(list())
  out <- list()
  for (nm in names(x)) out[[nm]] <- as.character(unlist(x[[nm]], use.names = FALSE))
  out
}

apply_celltype_map <- function(obj, cmap, levels_vec = NULL, overwrite_col = "celltype") {
  cmap <- normalize_map(cmap)
  if (length(cmap) == 0) return(obj)
  md <- obj@meta.data
  if (!"seurat_clusters" %in% colnames(md)) stop("Missing seurat_clusters in metadata.", call. = FALSE)
  if (overwrite_col %in% colnames(md)) md[[overwrite_col]] <- as.character(md[[overwrite_col]]) else md[[overwrite_col]] <- NA_character_
  sc <- as.character(md$seurat_clusters)
  for (ct in names(cmap)) {
    idx <- which(sc %in% cmap[[ct]])
    if (length(idx) > 0) md[idx, overwrite_col] <- ct
  }
  if (!is.null(levels_vec) && length(levels_vec) > 0) {
    md[[overwrite_col]] <- factor(md[[overwrite_col]], levels = as.character(levels_vec))
  } else {
    md[[overwrite_col]] <- factor(md[[overwrite_col]])
  }
  obj@meta.data <- md
  obj
}

subset_excluding_clusters <- function(obj, clusters) {
  if (is.null(clusters) || length(clusters) == 0) return(obj)
  md <- obj@meta.data
  keep <- rownames(md)[!(as.character(md$seurat_clusters) %in% as.character(clusters))]
  subset(obj, cells = keep)
}

run_harmony_chain <- function(obj, dims = 1:40, res = 1.0, npcs = 50, group_by = "orig.ident",
                              regress_vars = NULL, do_scale = TRUE, max_value = 10) {
  if (do_scale) {
    if (!is.null(regress_vars) && length(regress_vars) > 0) {
      obj <- ScaleData(obj, vars.to.regress = regress_vars, max.value = max_value)
    } else {
      obj <- ScaleData(obj, max.value = max_value)
    }
  }
  obj <- RunPCA(obj, npcs = npcs)
  obj <- RunHarmony(obj, group.by.vars = group_by, reduction = "pca")
  obj <- FindNeighbors(obj, reduction = "harmony", dims = dims)
  obj <- FindClusters(obj, reduction = "harmony", resolution = res)
  obj <- RunUMAP(obj, reduction = "harmony", dims = dims)
  obj
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

cfg <- yaml::read_yaml(args$config)
in_rds <- as.character(cfg$input_rds %||% "results/Fig4/obj_all_withanno.rds")
out_dir <- as.character(cfg$out_dir %||% "results/Fig4/gc_processing")
prefix <- as.character(cfg$prefix %||% "gc")
if (!file.exists(in_rds)) stop("Missing input_rds: ", in_rds, call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, paste0(prefix, ".fig4_gc_processing.log"))
log_file <- file(log_path, open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log_file), silent = TRUE)
}, add = TRUE)

message("[fig4-gc] start: ", Sys.time())
message("[fig4-gc] input_rds: ", in_rds)
message("[fig4-gc] out_dir: ", out_dir)

obj_total <- readRDS(in_rds)
if (!inherits(obj_total, "Seurat")) stop("input_rds must be a Seurat object.", call. = FALSE)
if (!"celltype" %in% colnames(obj_total@meta.data)) stop("input_rds is missing metadata column `celltype`.", call. = FALSE)
if (!"seurat_clusters" %in% colnames(obj_total@meta.data)) stop("input_rds is missing metadata column `seurat_clusters`.", call. = FALSE)

obj_total$celltype <- as.character(obj_total$celltype)
cluster_relabel <- cfg$cluster_relabel %||% list("9" = "GC_Column")
for (cl in names(cluster_relabel)) {
  idx <- which(as.character(obj_total@meta.data$seurat_clusters) == as.character(cl))
  if (length(idx) > 0) obj_total@meta.data$celltype[idx] <- as.character(cluster_relabel[[cl]])
}

subset_celltypes <- as.character(unlist(cfg$subset_celltypes %||% c("GC_Antral_Mural", "GC_Preantral", "GC_Atretic", "GC_Column")))
keep_cells <- rownames(obj_total@meta.data)[obj_total@meta.data$celltype %in% subset_celltypes]
obj_gr <- subset(obj_total, cells = keep_cells)
message("[fig4-gc] subset cells: ", ncol(obj_gr))

split_by <- as.character(cfg$split_by %||% "orig.ident")
if (!split_by %in% colnames(obj_gr@meta.data)) stop("Missing split_by column: ", split_by, call. = FALSE)
Idents(obj_gr) <- split_by
obj.list <- SplitObject(obj_gr)
obj.list <- lapply(obj.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = as.integer(cfg$nfeatures_select %||% 2000))
  x
})
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = as.integer(cfg$nfeatures_select %||% 2000))
exclude_regex <- as.character(cfg$exclude_regex %||% "^(Rp|mt)")
if (!is.na(exclude_regex) && nzchar(exclude_regex)) {
  fl_genes <- grep(exclude_regex, rownames(obj_gr), value = TRUE)
  features <- setdiff(features, fl_genes)
}
obj_gr <- NormalizeData(obj_gr)
obj_gr@assays$RNA@var.features <- features

s.genes <- as.character(unlist(cfg$s_genes %||% character(0)))
g2m.genes <- as.character(unlist(cfg$g2m_genes %||% character(0)))
s.genes <- intersect(s.genes, rownames(obj_gr))
g2m.genes <- intersect(g2m.genes, rownames(obj_gr))
if (length(s.genes) > 0 && length(g2m.genes) > 0) {
  obj_gr <- CellCycleScoring(obj_gr, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
}

# Round 1: no cell-cycle regression, follow notebook logic
obj_gr <- ScaleData(obj_gr)
obj_gr <- RunPCA(obj_gr, npcs = as.integer(cfg$round1_npcs %||% 50), features = features)
obj_gr <- RunHarmony(obj_gr, group.by.vars = split_by, reduction = "pca")
obj_gr <- FindNeighbors(obj_gr, reduction = "harmony", dims = seq_len(as.integer(cfg$round1_dim_max %||% 40)))
obj_gr <- FindClusters(obj_gr, reduction = "harmony", resolution = as.numeric(cfg$round1_resolution_1 %||% 0.5))
obj_gr <- RunUMAP(obj_gr, reduction = "harmony", dims = seq_len(as.integer(cfg$round1_dim_max %||% 40)))
obj_gr <- FindNeighbors(obj_gr, reduction = "harmony", dims = seq_len(as.integer(cfg$round1_dim_max %||% 40)))
obj_gr <- FindClusters(obj_gr, reduction = "harmony", resolution = as.numeric(cfg$round1_resolution_2 %||% 1.5))
obj_gr <- RunUMAP(obj_gr, reduction = "harmony", dims = seq_len(as.integer(cfg$round1_dim_max %||% 40)))

obj_round1 <- obj_gr
out_round1_rds <- as.character(cfg$out_round1_rds %||% file.path(out_dir, paste0(prefix, ".round1.rds")))
saveRDS(obj_round1, out_round1_rds)

# Filter bad clusters after round1
obj_gr <- subset_excluding_clusters(obj_gr, cfg$exclude_clusters_round1 %||% c("6", "7", "15", "16"))
message("[fig4-gc] after exclude round1 clusters: ", ncol(obj_gr))

# Intermediate annotation (for traceability; downstream uses final labels)
obj_gr <- apply_celltype_map(
  obj_gr,
  cfg$celltype_map_intermediate %||% list(
    GC_Progenitor = c("11"),
    GC_Preantral = c("17", "5", "2", "0", "14"),
    GC_Cumulus = c("18", "10", "9", "12"),
    GC_Antral_Mural_Inhba = c("8"),
    GC_Antral_Mural_Id2 = c("4"),
    GC_Antral_Mural_Slc38a5 = c("1"),
    GC_Atretic = c("3", "13")
  ),
  levels_vec = as.character(unlist(cfg$celltype_levels_intermediate %||% c(
    "GC_Progenitor", "GC_Preantral", "GC_Cumulus",
    "GC_Antral_Mural_Inhba", "GC_Antral_Mural_Id2", "GC_Antral_Mural_Slc38a5", "GC_Atretic"
  )))
)

# Recluster after filtering
obj_gr <- FindNeighbors(obj_gr, reduction = "harmony", dims = seq_len(as.integer(cfg$round1_dim_max %||% 40)))
obj_gr <- FindClusters(obj_gr, reduction = "harmony", resolution = as.numeric(cfg$intermediate_resolution %||% 1.0))
obj_gr <- RunUMAP(obj_gr, reduction = "harmony", dims = seq_len(as.integer(cfg$round1_dim_max %||% 40)))

# Pre-final label set
obj_gr <- apply_celltype_map(
  obj_gr,
  cfg$celltype_map_prefinal %||% list(
    GC_Preantral = c("9", "10", "1", "2"),
    GC_Mitotic_S = c("5"),
    GC_Mitotic_G2M = c("6", "8"),
    GC_Cumulus = c("7"),
    GC_Antral_Mural = c("4"),
    GC_Atretic = c("0", "3", "11")
  ),
  levels_vec = as.character(unlist(cfg$celltype_levels_prefinal %||% c(
    "GC_Preantral", "GC_Mitotic_S", "GC_Mitotic_G2M", "GC_Cumulus", "GC_Antral_Mural", "GC_Atretic"
  )))
)

# Final rerun with cell-cycle regression
regress_vars <- as.character(unlist(cfg$final_regress_vars %||% c("S.Score", "G2M.Score")))
obj_gr <- run_harmony_chain(
  obj_gr,
  dims = seq_len(as.integer(cfg$final_dim_max %||% 30)),
  res = as.numeric(cfg$final_resolution %||% 1.0),
  npcs = as.integer(cfg$final_npcs %||% 50),
  group_by = split_by,
  regress_vars = regress_vars,
  do_scale = TRUE,
  max_value = as.numeric(cfg$scale_max_value %||% 10)
)

obj_gr <- apply_celltype_map(
  obj_gr,
  cfg$celltype_map_final %||% list(
    GC_Progenitor = c("7"),
    GC_Preantral = c("11", "0", "1", "8"),
    GC_Mitotic = c("9"),
    GC_Cumulus = c("5", "10"),
    GC_Antral_Mural = c("6", "3"),
    GC_Atretic = c("2", "4")
  ),
  levels_vec = as.character(unlist(cfg$celltype_levels_final %||% c(
    "GC_Progenitor", "GC_Preantral", "GC_Mitotic", "GC_Cumulus", "GC_Antral_Mural", "GC_Atretic"
  )))
)

# Remove final low-quality cluster and re-embed/recluster, keeping final celltype labels
obj_gr <- subset_excluding_clusters(obj_gr, cfg$exclude_clusters_final %||% c("8"))
obj_gr <- RunUMAP(obj_gr, reduction = "harmony", dims = seq_len(as.integer(cfg$final_dim_max %||% 30)))
obj_gr <- FindNeighbors(obj_gr, reduction = "harmony", dims = seq_len(as.integer(cfg$final_dim_max %||% 30)))
obj_gr <- FindClusters(obj_gr, reduction = "harmony", resolution = as.numeric(cfg$postfilter_resolution %||% 1.0))

out_final_rds <- as.character(cfg$out_final_rds %||% file.path(out_dir, "obj_gr_newanno.rds"))
saveRDS(obj_gr, out_final_rds)

# lightweight plots
pdf(file.path(out_dir, paste0(prefix, ".umap_by_celltype.pdf")), width = 6.5, height = 5.5)
print(DimPlot(obj_gr, reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend())
dev.off()

pdf(file.path(out_dir, paste0(prefix, ".umap_by_cluster.pdf")), width = 6.5, height = 5.5)
print(DimPlot(obj_gr, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
dev.off()

yaml::write_yaml(cfg, file.path(out_dir, "params_used_fig4_gc_processing.yaml"))
write_text(file.path(out_dir, "sessionInfo_fig4_gc_processing.txt"), capture.output(sessionInfo()))
message("[fig4-gc] final cells: ", ncol(obj_gr))
message("[fig4-gc] done: ", Sys.time())
