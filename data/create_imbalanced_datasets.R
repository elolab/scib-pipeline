#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to create imbalanced PBMC cell types for benchmarking 
#integration through the scib-pipeline. 
# Date: 18/07/2025
# Last update: 18/07/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Packages
library("SingleCellExperiment")
library("Seurat")
library("anndata")

## Process data 
# Import data
data.dir <- file.path("data", "08_extdata_figure_6")
input <- list(
  batch1 = file.path(data.dir, "tran_exp5_pbmc_batch1_balanced.h5ad"), 
  batch2 = file.path(data.dir, "tran_exp5_pbmc_batch2_balanced.h5ad")
)
sce <- lapply(input, function(x) zellkonverter::readH5AD(file = x, X_name = "counts"))

# Concatenate SCE batch datasets
stopifnot(all(row.names(sce$batch1) == row.names(sce$batch2)))
sce <- cbind(sce$batch1, sce$batch2)

# Parse variables
sce$batch <- as.character(sce$batch)
sce$celltype <- as.character(sce$celltype)

## Create imbalanced datasets by iteratively down sampling cell types from batch_2 to
#0, 5, 10, and 100% of the cells, i.e., from being completely absent in one 
#of the batches to completely balanced. 
cell.types <- c("B cell", "Monocyte_CD14", "Monocyte_FCGR3A")
down <- c(0, 0.05, 0.1) 
batch.down <- "batch_2"
i <- 1
for (cell in cell.types) {
  for (d in down) {
    if (i == 1) {
      cat("\n Creating balanced dataset\n")
      tmp.sce <- sce
      # Remove lowly & non-expressed genes & normalize
      pick.genes <- (rowSums(counts(tmp.sce))>10) # remove non-expressing genes
      tmp.sce <- tmp.sce[pick.genes,sort(colnames(tmp.sce))] # pick genes & sort cell names
      tmp.sce <- as.SingleCellExperiment(NormalizeData(object = as.Seurat(tmp.sce, data = NULL)))
      # Format data for anndata 
      counts(tmp.sce) <- as(as(counts(tmp.sce), "RsparseMatrix"), "dgRMatrix")
      logcounts(tmp.sce) <- as(as(logcounts(tmp.sce), "RsparseMatrix"), "dgRMatrix")
      # Create anndata object
      cat("\n Selecting:", nrow(tmp.sce), "genes x", ncol(tmp.sce), "cells\n")
      adata <- AnnData(
        X = t(logcounts(tmp.sce)),
        obs = as.data.frame(colData(tmp.sce)),
        layers = list(counts = t(counts(tmp.sce)))
      )
      # Export anndata object as h5ad
      adata$write_h5ad(filename = file.path("datasets", "balanced_pbmc.h5ad"))
      rm(list = c("adata", "tmp.sce")); gc();
    }
    cat("\n Creating imbalanced dataset for cell type:", cell, paste0("(downsampled to: ", d*100, "%)\n"))
    tmp.sce <- sce
    name <- paste(gsub(" ", "_", cell), d, sep = "_")
    # Downsample
    down.label <- colnames(tmp.sce)[ (tmp.sce$batch==batch.down & tmp.sce$celltype == cell) ]
    other.labels <- colnames(tmp.sce)[ ! (colnames(tmp.sce) %in% down.label) ]
    set.seed(289)
    pick.cells <- sample(x = down.label, size = length(down.label)*d)
    pick.cells <- sort(c(pick.cells, other.labels))
    tmp.sce <- tmp.sce[,pick.cells]
    # Remove lowly & non-expressed genes & normalize
    pick.genes <- (rowSums(counts(tmp.sce))>10) # remove non-expressing genes
    tmp.sce <- tmp.sce[pick.genes,]
    tmp.sce <- as.SingleCellExperiment(NormalizeData(object = as.Seurat(tmp.sce, data = NULL)))
    # Format data for anndata 
    counts(tmp.sce) <- as(as(counts(tmp.sce), "RsparseMatrix"), "dgRMatrix")
    logcounts(tmp.sce) <- as(as(logcounts(tmp.sce), "RsparseMatrix"), "dgRMatrix")
    # Create anndata object
    cat("\n Selecting:", nrow(tmp.sce), "genes x", ncol(tmp.sce), "cells\n")
    adata <- AnnData(
      X = t(logcounts(tmp.sce)),
      obs = as.data.frame(colData(tmp.sce)),
      layers = list(counts = t(counts(tmp.sce)))
    )
    # Export anndata object as h5ad
    adata$write_h5ad(filename = file.path("datasets", paste0("imbalanced_pbmc_", name, ".h5ad")))
    rm(list = c("adata", "tmp.sce")); gc();
    i <- i + 1
  }
}
#------------------------------------------------------------------------------#
