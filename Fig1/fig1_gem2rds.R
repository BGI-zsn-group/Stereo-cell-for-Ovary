#!/usr/bin/env Rscript

# gem2rds (Fig.1 part-2): Convert cut GEM to Seurat RDS
#
# Input:
#   1) <prefix>.QuPath.cut.gem.gz (output of cutgem.py)
# Output:
#   Seurat object saved as .rds
#
# Usage:
#   Rscript gem2rds.R <cut_gem_path> <out_rds_path> <project_name>
#
# Example:
#   Rscript gem2rds.R gem_cut/sample1.QuPath.cut.gem.gz fig1_sample1.rds sample1

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript gem2rds.R <cut_gem_path> <out_rds_path> <project_name>\n",
       "Example: Rscript gem2rds.R gem_cut/sample1.QuPath.cut.gem.gz sample1.rds sample1")
}

in_gem <- args[1]
out_rds <- args[2]
project <- args[3]

message("[gem2rds] reading: ", in_gem)

# Robust read for .gz
if (grepl("\\.gz$", in_gem, ignore.case = TRUE)) {
  gem <- fread(cmd = paste("zcat", shQuote(in_gem)), header = TRUE)
} else {
  gem <- fread(in_gem, header = TRUE)
}

required_cols <- c("geneID", "MIDCount", "Label")
missing_cols <- setdiff(required_cols, colnames(gem))
if (length(missing_cols) > 0) {
  stop("[gem2rds] Missing required columns: ", paste(missing_cols, collapse = ", "))
}

gem$label <- as.character(gem$Label)

genes <- unique(gem$geneID)
cells <- unique(gem$label)

gene_index <- seq_along(genes); names(gene_index) <- genes
cell_index <- seq_along(cells); names(cell_index) <- cells

# Build sparse matrix (duplicates are summed by sparseMatrix)
mat1 <- sparseMatrix(
  i = gene_index[gem$geneID],
  j = cell_index[gem$label],
  x = gem$MIDCount,
  dims = c(length(genes), length(cells))
)
rownames(mat1) <- genes
colnames(mat1) <- cells

message("[gem2rds] matrix built: genes=", nrow(mat1), " cells=", ncol(mat1), " nnzero=", length(mat1@x))

obj <- CreateSeuratObject(counts = mat1, project = project)
saveRDS(obj, file = out_rds, compress = TRUE)

message("[gem2rds] saved: ", out_rds)
message("[gem2rds] done")
