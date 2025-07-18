#------------------------------------------------------------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to create ifnb h5ad file for scib-pipeline
# Date: 07/06/2024
# Last update: 21/01/2025
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Import packages
library("SingleCellExperiment")
library("Seurat")
library("anndata")

# Import data
ifnb <- SeuratData::LoadData("ifnb")

# Normalize & remove non-expressed genes
ifnb <- NormalizeData(ifnb)
sce.ifnb <- as.SingleCellExperiment(ifnb)
sce.ifnb <- Coralysis::PrepareData(sce.ifnb)

# Remove cells from one of the batches
pick.cells <- ((sce.ifnb$stim=="CTRL" & !(sce.ifnb$seurat_annotations %in% c("CD14 Mono", "CD4 Naive T"))) | 
	                        (sce.ifnb$stim=="STIM" & !(sce.ifnb$seurat_annotations %in% c("CD16 Mono", "CD4 Memory T"))))
sce.ifnb <- sce.ifnb[,pick.cells]

# Format data for anndata 
counts(sce.ifnb) <- as(as(counts(sce.ifnb), "RsparseMatrix"), "dgRMatrix")
logcounts(sce.ifnb) <- as(as(logcounts(sce.ifnb), "RsparseMatrix"), "dgRMatrix")
sce.ifnb$seurat_annotations <- as.character(sce.ifnb$seurat_annotations)
sce.ifnb$ident <- as.character(sce.ifnb$ident)
sce.ifnb$orig.ident <- as.character(sce.ifnb$orig.ident)

# Create anndata object
adata <- AnnData(
		 X = t(logcounts(sce.ifnb)),
		 obs = as.data.frame(colData(sce.ifnb)),
		 layers = list(counts = t(counts(sce.ifnb)))
)

# Export anndata object as h5ad 
adata$write <- h5ad(filename = "datasets/ifnb.h5ad")
#
#------------------------------------------------------------------------------#
