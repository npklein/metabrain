import sys
import argparse
import gzip

parser = argparse.ArgumentParser(description='Annotate Sieberts eQTLs with TSS distance after liftover. The file format is quite specific so this is not reusable')
parser.add_argument('eqtl_file', help='File with columns as in description')
parser.add_argument('gtf', help='GTF file')
parser.add_argument('outfile', help='Name of outfile')

args = parser.parse_args()

#Loading data
gene_tss_pos = {}
with gzip.open(args.gtf,'rt') as gtf_file:
    for line in gtf_file:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] != 'gene':
            continue
        start = line[3]
        end = line[4]
        strand = line[6]
        gene = line[8].split('gene_id "')[1].split('.')[0]
        if strand == '-':
            gene_tss_pos[gene] = int(end)
        elif strand == '+':
            gene_tss_pos[gene] = int(start)
        else:
            raise RuntimeError('strand should be - or +, was: '+strand) 

top_eQTL = {}
lowest_fdr = {}
with open(args.eqtl_file) as input_file:
    for line in input_file:
        line = line.rstrip().split('\t')
        SNP_pos = int(line[2])
        gene = line[3].split(',')[4]
        FDR = float(line[3].split(',')[8])
        tss_dist = abs(SNP_pos)-gene_tss_pos[gene]
        line = '\t'.join(line)+','+str(tss_dist)+'\n'
        if gene not in top_eQTL:
            top_eQTL[gene] = line
            lowest_fdr = FDR
        elif FDR < lowest_fdr:
            top_eQTL[gene] = line
            lowest_fdr = FDR


with open(args.outfile,'w') as out:
    for gene in top_eQTL:
        out.write(top_eQTL[gene])

