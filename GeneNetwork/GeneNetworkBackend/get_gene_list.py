import os
from datetime import datetime
import argparse
parser = argparse.ArgumentParser(description='Make gene list based on GTF file that was used for expression quantification.')
parser.add_argument('gtf_file',
                    help='GTF file used for mapping expression')
args = parser.parse_args()


today = datetime.now().strftime("%Y-%m-%d")
genes = set([])
with open(args.gtf_file) as input_file:
    for line in input_file:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] == 'gene':
            ensembl_id = line[8].split('gene_id "')[1].split('.')[0]
            genes.add(ensembl_id)

genes_ordered = sorted(genes)

with open(today+'-'+os.path.basename(args.gtf_file).replace('.gtf','.orderedGeneList.txt'),'w') as out:
    for gene in genes_ordered:
        out.write(gene+'\n')
