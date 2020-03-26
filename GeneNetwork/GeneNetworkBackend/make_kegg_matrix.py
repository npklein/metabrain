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
parser.add_argument('kegg_version',
                    help='Kegg version to use (e.g.: 7.0)')
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
input_file = today+'-c2.cp.kegg.v'+args.kegg_version+'.entrez.gmt.gz'

if not os.path.exists(input_file):
    print('Getting kegg version', args.kegg_version)
    url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"+args.kegg_version+"/c2.cp.kegg.v"+args.kegg_version+".entrez.gmt"
    wget.download(url, out=input_file.rstrip('.gz'))
    with open(input_file.rstrip('.gz'), 'rb') as f_in, gzip.open(input_file, 'wb') as f_out:
        f_out.writelines(f_in)

    os.remove(input_file.rstrip('.gz'))

pathway_genes = {}
pathways = set([])

print('Start reading '+input_file)
with  gzip.open(input_file,'rt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        pathway = ' '.join(line[0].replace('KEGG_','').lower().split('_')).capitalize()
        ncbi_genes = line[2:]
        pathway_genes[pathway] = set([ncbi_to_ensembl[x] for x in ncbi_genes if x in ncbi_to_ensembl])
        pathways.add(pathway)
print('done')

pathways = sorted(pathways)
outfile = 'PathwayMatrix/'+today+'c2.cp.kegg.v'+args.kegg_version+'.matrix.txt'
print('start writing matrix')
with open(args.ordered_gene_list) as input_file, open(outfile,'w') as out, open(outfile.replace('matrix.txt','genesInPathways.txt'),'w') as out2:
    out.write(today)
    for pathway in pathways:
        out.write('\t'+pathway)
    out.write('\n')
    gene_in_at_least_1_kegg = set([])
    for gene in input_file:
        gene = gene.strip()
        out.write(gene)
        for pathway in pathways:
            out.write('\t')
            if gene in pathway_genes[pathway]:
                out.write('1')
                gene_in_at_least_1_kegg.add(gene)
            else:
                out.write('0')
        out.write('\n')
    for gene in gene_in_at_least_1_kegg:
        out2.write(gene+'\n')
print('Output written to '+outfile)
