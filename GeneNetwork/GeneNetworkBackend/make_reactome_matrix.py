from datetime import datetime
import wget
from time import strftime
import os
from pathlib import Path
import gzip
import argparse

parser = argparse.ArgumentParser(description='Make REACTOME matrix.')
parser.add_argument('ordered_gene_list',
                    help='List with ordered gene IDs')
args = parser.parse_args()
Path("PathwayMatrix/").mkdir(parents=True, exist_ok=True)


today = datetime.now().strftime("%Y-%m-%d") 
input_file = today+'-Ensembl2Reactome_All_Levels.txt.gz'
if not os.path.exists(input_file):
    print('Getting latest reactome annotation')
    url = "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt"
    wget.download(url, out=today+'-Ensembl2Reactome_All_Levels.txt')
    with open(today+'-Ensembl2Reactome_All_Levels.txt', 'rb') as f_in, gzip.open(input_file, 'wb') as f_out:
        f_out.writelines(f_in)

    os.remove(today+'-Ensembl2Reactome_All_Levels.txt')


pathway_genes = {}
pathways = set([])

print('Start reading '+input_file)
with  gzip.open(input_file,'rt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        if line[5] != 'Homo sapiens':
            continue
        ensembl_id = line[0]
        pathway = line[1]
        if pathway not in pathway_genes:
            pathway_genes[pathway] = set([])
        pathway_genes[pathway].add(ensembl_id)
        pathways.add(pathway)
print('done')

pathways = sorted(pathways)
outfile = 'PathwayMatrix/'+today+'-Ensembl2Reactome_All_Levels.matrix.txt'
print('start writing matrix')
with open(args.ordered_gene_list) as input_file, open(outfile,'w') as out, open(outfile.replace('matrix.txt','genesInPathways.txt'),'w') as out2:
    out.write(today)
    for pathway in pathways:
        out.write('\t'+pathway)
    out.write('\n')
    gene_in_atleast_1_reactome = set([])
    for gene in input_file:
        gene = gene.strip()
        out.write(gene)
        for pathway in pathways:
            out.write('\t')
            if gene in pathway_genes[pathway]:
                out.write('1')
                gene_in_atleast_1_reactome.add(gene)
            else:
                out.write('0')
        out.write('\n')
    for gene in gene_in_atleast_1_reactome:
        out2.write(gene+'\n')
print('Output written to '+outfile)
