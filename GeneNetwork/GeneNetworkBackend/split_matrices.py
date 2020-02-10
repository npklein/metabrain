# Split GeneNetwork input matrices up into smaller chunks so that they can be processed in parallel

import os
import math
from multiprocessing import Pool
hardcoded_path = '/groups/umcg-biogen/tmp04/output/2019-02-25-FreezeTwo/output/2019-07-29-MetaBrainNetwork/GeneNetworkBackend/'



def split_files(i,f):
    initial_start = 1
    initial_end = 101
    
    start = initial_start + (initial_end*i)
    end = initial_end + (initial_end*i)
    k = f.rfind(".txt")
    with open(f) as input_file:
        n_columns = len(input_file.readline().split('\t'))
    out_name = f.split('/')[-1][:k]+'.start'+str(start)+'end'+str(end)+'.txt'
    print(out_name)
    with open(f) as input_file, open(hardcoded_path+'PathwayMatrix/split_matrices/'+type+'/'+out_name,'w') as out:
        for line in input_file:
            line = line.strip().split('\t')
            out.write(line[0])
            for y in range(start, end+1, 1):
                if y >= n_columns:
                    break
                out.write('\t'+line[y])
            out.write('\n')


def run(type):
    pathway_files = {'hpo':hardcoded_path+'PathwayMatrix/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix_filteredOnEigenVectorGenes.txt',
                 'reactome':hardcoded_path+'PathwayMatrix/Ensembl2Reactome_All_Levels.txt_matrix_filteredOnEigenVectorGenes.txt',
                 'go_F':hardcoded_path+'PathwayMatrix/goa_human.gaf_F_matrix_filteredOnEigenVectorGenes.txt',
                 'go_C':hardcoded_path+'PathwayMatrix/goa_human.gaf_C_matrix_filteredOnEigenVectorGenes.txt',
                 'go_P':hardcoded_path+'PathwayMatrix/goa_human.gaf_P_matrix_filteredOnEigenVectorGenes.txt',
                 'kegg':hardcoded_path+'PathwayMatrix/c2.cp.kegg.v6.1.entrez.gmt_matrix_filteredOnEigenVectorGenes.txt'}
    if not os.path.exists(hardcoded_path+'PathwayMatrix/split_matrices/'+type):
        os.makedirs(hardcoded_path+'PathwayMatrix/split_matrices/'+type)
    
    f = pathway_files[type]
    with open(f) as input_file:
        n_columns = len(input_file.readline().split('\t'))
    batches = math.ceil(n_columns/100)

    input_data = []
    for i in range(0, batches, 1):
        input_data.append((i, f))
    p = Pool(23)

    p.starmap(split_files, input_data)

types = ['kegg','hpo','reactome','go_F','go_C','go_P']
for type in types:
    run(type)
