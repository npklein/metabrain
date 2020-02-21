hardcoded_path = '/groups/umcg-biogen/tmp04/output/2019-02-25-FreezeTwo/output/2019-07-29-MetaBrainNetwork/GeneNetworkBackend/'

pathway_files = [hardcoded_path+'PathwayMatrix/goa_human.gaf_F_matrix.txt',hardcoded_path+'PathwayMatrix/goa_human.gaf_F_genesInPathways.txt',hardcoded_path+'PathwayMatrix/goa_human.gaf_P_matrix.txt',
                 hardcoded_path+'PathwayMatrix/goa_human.gaf_P_genesInPathways.txt',hardcoded_path+'PathwayMatrix/c2.cp.kegg.v6.1.entrez.gmt_matrix.txt',
                 hardcoded_path+'PathwayMatrix/c2.cp.kegg.v6.1.entrez.gmt_genesInPathways.txt',hardcoded_path+'PathwayMatrix/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt',
                 hardcoded_path+'PathwayMatrix/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_genesInPathways.txt',hardcoded_path+'PathwayMatrix/Ensembl2Reactome_All_Levels.txt_matrix.txt',
                 hardcoded_path+'PathwayMatrix/Ensembl2Reactome_All_Levels.txt_genesInPathways.txt',hardcoded_path+'PathwayMatrix/goa_human.gaf_C_matrix.txt',hardcoded_path+'PathwayMatrix/goa_human.gaf_C_genesInPathways.txt']

vector_per_gene = {}
eigenvectors='/groups/umcg-biogen/tmp04/umcg-ndeklein/GeneNetwork/output/7_PCA_on_correlation_matrix/MetaBrain.eigenvectors.cronbach_0.9.txt'
print('parse eigenvectors')
with open(eigenvectors) as input_file:
    header = input_file.readline()
    for line in input_file:
        gene =line.split('\t')[0]
        vector_per_gene[gene] = line

print('done parsing eigenvectors')

gene_order = []
for index, f in enumerate(pathway_files):
    if index == 0:
        eigen_out = open(eigenvectors.replace('.txt','.filteredOnPathwayGenes.txt'),'w')
        eigen_out.write(header)
    k = f.rfind(".txt")
    out_name = f[:k] + "_filteredOnEigenVectorGenes.txt"
    print(f)
    print(out_name)
    print('-'*20)
    
    with open(f) as input_file, open(out_name,'w') as out:
        if not f.endswith('genesInPathways.txt'):
            out.write(input_file.readline())
        gene_index = 0
        for line in input_file:
            gene = line.strip().split('\t')[0]
            if gene in vector_per_gene:
                out.write(line)
                if index == 0:
                    eigen_out.write(vector_per_gene[gene])
                    gene_order.append(gene)
                else:
                    if f.endswith('matrix.txt'):
                        assert gene_index == gene_order.index(gene)
                gene_index += 1
    if index == 0:
        eigen_out.close()
