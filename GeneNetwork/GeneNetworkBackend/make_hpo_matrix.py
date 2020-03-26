from datetime import datetime
import wget
from time import strftime
import os
from pathlib import Path
import gzip
import argparse

parser = argparse.ArgumentParser(description='Make HPO matrix.')
parser.add_argument('ordered_gene_list',
                    help='List with ordered gene IDs')
parser.add_argument('ncbi_to_ensembl_file',
                    help='Gzipped file with in first column ensembl IDs, second column NCBI IDs')
args = parser.parse_args()
Path("PathwayMatrix/").mkdir(parents=True, exist_ok=True)

ncbi_to_ensembl = {}
with gzip.open(args.ncbi_to_ensembl_file,'rt') as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip('\n').split('\t')
        ncbi_to_ensembl[line[1]] = line[0]

today = datetime.now().strftime("%Y-%m-%d") 
input_file_name = today+'-HPO-phenotype-to-genes.txt.gz'

if not os.path.exists(input_file_name):
    print('Getting lates hpo build')
    url = "http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/util/annotation/phenotype_to_genes.txt"
    wget.download(url, out=input_file_name.rstrip('.gz'))
    with open(input_file_name.rstrip('.gz'), 'rb') as f_in, gzip.open(input_file_name, 'wb') as f_out:
        f_out.writelines(f_in)

    os.remove(input_file.rstrip('.gz'))

pathway_genes = {}
pathways = set([])
print('Start reading '+input_file_name)
mapped = 0
not_mapped = 0
with  gzip.open(input_file_name,'rt') as input_file:
    input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        pathway = line[0]
        ncbi_gene = line[2]
        if pathway not in pathway_genes:
            pathway_genes[pathway] = set([])
        if ncbi_gene in ncbi_to_ensembl:
            mapped += 1
            pathway_genes[pathway].add(ncbi_to_ensembl[ncbi_gene])
        else:
            not_mapped += 1
        pathways.add(pathway)
print('Could not map ID',not_mapped,'times')
print('Could map ID',mapped,'times')

print('done')

pathways = sorted(pathways)
outfile = 'PathwayMatrix/'+input_file_name.replace('.txt.gz','')+'.matrix.txt'
print('start writing matrix')
with open(args.ordered_gene_list) as input_file, open(outfile,'w') as out, open(outfile.replace('matrix.txt','genesInPathways.txt'),'w') as out2:
    out.write(today)
    for pathway in pathways:
        out.write('\t'+pathway)
    out.write('\n')
    gene_in_atleast_1_hpo = set([])
    for gene in input_file:
        gene = gene.strip()
        out.write(gene)
        for pathway in pathways:
            out.write('\t')
            if gene in pathway_genes[pathway]:
                out.write('1')
                gene_in_atleast_1_hpo.add(gene)
            else:
                out.write('0')
        out.write('\n')
    for gene in gene_in_atleast_1_hpo:
        out2.write(gene+'\n')
print('Output written to '+outfile)
