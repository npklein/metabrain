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
args = parser.parse_args()
Path("PathwayMatrix/").mkdir(parents=True, exist_ok=True)


today = datetime.now().strftime("%Y-%m-%d") 
input_file = today+'-c2.cp.kegg.v'+args.kegg_version+'.entrez.gmt.gz'

if not os.path.exists(input_file):
    print('Getting kegg version', args.kegg_version)
    url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"+args.kegg_version+"/c2.cp.kegg.v"+args.kegg_version+".entrez.gmt"
    wget.download(url, out=input_file.rstrip('.gz'))
    with open(input_file.rstrip('.gz'), 'rb') as f_in, gzip.open(input_file, 'wb') as f_out:
        f_out.writelines(f_in)

    os.remove(input_file.rstrip('.gz'))
exit()

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
with open(args.ordered_gene_list) as input_file, open(outfile,'w') as out:
    out.write(today)
    for pathway in pathways:
        out.write('\t'+pathway)
    out.write('\n')
    for gene in input_file:
        gene = gene.strip()
        out.write(gene)
        for pathway in pathways:
            out.write('\t')
            if gene in pathway_genes[pathway]:
                out.write('1')
            else:
                out.write('0')
        out.write('\n')
print('Output written to '+outfile)
