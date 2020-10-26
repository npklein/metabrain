# Title     : create_seurat_object.R
# Objective :
# Created by: mvochteloo
# Created on: 22/10/2020

# install.packages('Seurat')
library(Seurat)

data_dir <- '/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/2020-10-22-MVochteloo-Copy'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = expression_matrix)


