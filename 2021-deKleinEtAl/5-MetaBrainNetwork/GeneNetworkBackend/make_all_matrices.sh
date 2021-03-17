

ml Python
python make_go_matrices.py 2020-03-27-eigenvector-genes.txt 2020-03-25-uniprot_ensembl_map.txt.gz &
python make_hpo_matrix.py 2020-03-27-eigenvector-genes.txt 2020-03-25-ncbi-to-ensembl.txt.gz &
python make_kegg_matrix.py 2020-03-27-eigenvector-genes.txt 7.0 2020-03-25-ncbi-to-ensembl.txt.gz &
python make_reactome_matrix.py 2020-03-27-eigenvector-genes.txt &
