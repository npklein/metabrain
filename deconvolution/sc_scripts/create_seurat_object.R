# Title     : create_seurat_object.R
# Objective :
# Created by: mvochteloo
# Created on: 22/10/2020

# https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
# https://satijalab.org/seurat/v3.2/sctransform_vignette.html

# install.packages('Seurat')
library(Seurat)

work_dir <- '/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/2020-10-22-MVochteloo-Copy/'
list.files(work_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = work_dir)

projid_table <- read.table('/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/filtered_column_metadata.txt.gz', header = T, stringsAsFactors = F)
projid_table$TAG <- gsub('\\.', '-', projid_table$TAG)

projid_rosmapid_translate_table <- read.csv("/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/meta/ROSMAP_IDkey.csv", sep = ",")
phenotype_table <- read.csv("/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt", sep="\t")
meta_table <- merge(projid_table, merge(projid_rosmapid_translate_table, phenotype_table, by.x="rnaseq_id", by.y="rnaseq_id"), by.x="projid", by.y="projid")
meta_table <- meta_table[, colSums(is.na(meta_table)) < nrow(meta_table)]

rownames(meta_table) <- meta_table$TAG
meta_table$TAG <- NULL

seurat_object <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                            min.cells = 0,
                                            min.features = 0,
                                            project = "scRNAseq",
                                            meta.data = meta_table)

filtered_gene_row_file <- read.csv("/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/filtered_gene_row_names.txt.gz", sep = "\t")
filtered_gene_row_tags <- filtered_gene_row_file[, 1]
length(filtered_gene_row_tags)
# [1] 17925

seurat_object <- seurat_object[filtered_gene_row_tags, !is.na(seurat_object@meta.data$projid)]
seurat_object
# An object of class Seurat
# 17352 features across 70634 samples within 1 assay
# Active assay: RNA (17352 features, 0 variable features)

seurat_object <- SCTransform(seurat_object)


seurat_object <- RunPCA(seurat_object)
seurat_object <- RunTSNE(object = seurat_object, reduction = "pca",  dims = 1:20)
seurat_object <- RunUMAP(object = seurat_object, reduction = "pca",  dims = 1:20)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:50)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

library(ggplot2)
for (name in c("pre.cluster", "broad.cell.type", "Subcluster",
                 "SpecificBrainRegion", "BroadBrainRegion", "cohort",
                 "MetaCohort", "predicted.brain.region",
                 "reannotated_diangosis", "seurat_clusters")) {
  DimPlot(seurat_object, reduction = 'tsne', group.by = name)
  ggsave(paste0("ROSMAP_scRNAseq_tsne_groupedBy", name ,".pdf"), width = 10, height = 10)
}

seurat_object@meta.data
visualisation_table <- merge(merge(merge(seurat_object@meta.data, seurat_object[['pca']]@cell.embeddings, by="row.names", all=TRUE), seurat_object[['tsne']]@cell.embeddings, by="row.names", all=TRUE), seurat_object[['umap']]@cell.embeddings, by="row.names", all=TRUE)
visualisation_table$Row.names <- NULL
visualisation_table$Row.names <- NULL
visualisation_table$Row.names <- NULL
write.table(visualisation_table,
            file = paste0(work_dir, "2020-11-20-metadata-withDimReduction.txt"),
            quote = F, sep = "\t", col.names = NA)

length(rownames(seurat_object))
# [1] 16718 = genes
length(colnames(seurat_object))
# [1] 38383 = cells

write.table(cell_counts,
            file = paste0(work_dir, "expression.txt"),
            quote = F, sep = "\t", col.names = NA)

Idents(object=seurat_object) <- "broad.cell.type"

library(Matrix)

gene_trans_df <- read.csv("/groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz", sep = "\t")[, c("ArrayAddress", "Symbol")]
gene_trans_df <- gene_trans_df[complete.cases(gene_trans_df),]

mean_expression_per_cell_type <- function(seurat, sample.id.column.name = "projid", gene_trans_df = NULL) {
  individuals <- unique(seurat@meta.data[, sample.id.column.name])
  idents <- unique(seurat@active.ident)

  cell_counts <- matrix(nrow = length(idents), ncol = length(individuals),
                        dimnames = list(idents, individuals))

  for (ident in idents) {

    cells_cell_type <- seurat[, seurat@active.ident == ident]

    mean_expression_matrix <- matrix(nrow = nrow(cells_cell_type), ncol = length(individuals),
                                     dimnames = list(rownames(cells_cell_type), individuals))

    for (individual in individuals) {
      if (sum(cells_cell_type@meta.data[, sample.id.column.name] == individual) == 0) {
        mean_expression_matrix[, toString(individual)] <- 0
        cell_counts[ident, toString(individual)] <- 0
      } else {
        cells_cell_type_individual <- cells_cell_type[, cells_cell_type@meta.data[, sample.id.column.name] == individual]

        mean_expression_matrix[, toString(individual)] <- rowMeans(cells_cell_type_individual)
        cell_counts[ident, toString(individual)] <- ncol(cells_cell_type_individual)
      }
      print(cell_counts[ident, toString(individual)])
    }

    if (!is.null(gene_trans_df)) {
      rownames(mean_expression_matrix) <- gene_trans_df[match(rownames(mean_expression_matrix), gene_trans_df[,2]), 1]
      mean_expression_matrix <- mean_expression_matrix[!is.na(rownames(mean_expression_matrix)), ]
    }

    write.table(mean_expression_matrix,
                file = paste0(work_dir, ident, "_expression", ".tsv"),
                quote = F, sep = "\t", col.names = NA)

  }

  write.table(cell_counts,
              file = paste0(work_dir, "cell_counts.txt"),
              quote = F, sep = "\t", col.names = NA)
}

mean_expression_per_cell_type(seurat_object, gene_trans_df = gene_trans_df)