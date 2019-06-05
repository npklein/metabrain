import pandas as pd;
import argparse
import gzip
parser = argparse.ArgumentParser(description='Compare eQTLgen to MetaBrain')
parser.add_argument('eQTL_file',
                    help='File with eQTL data')
parser.add_argument('gwas_file',
                    help='File with gwas data')
parser.add_argument('outfile',
                    help='outfile name')
args = parser.parse_args()

traits = {}
with open(args.gwas_file) as input_file:
    header = input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        trait = line[7]
        chr = line[11]
        pos = line[12]
        if chr+'_'+pos not in traits:
            traits[chr+'_'+pos] = set([])
        traits[chr+'_'+pos].add(trait)

with gzip.open(args.eQTL_file,'rt') as input_file, gzip.open(args.outfile,'wt') as out:
    header = input_file.readline().strip('\n')
    split_head = header.split('\t')
    fdr_index = split_head.index('FDR')
    out.write(header+'\tGWAS_TRAIT\n')
    for line in input_file:
        line = line.strip('\n')
        if float(line.split('\t')[fdr_index]) >= 0.05:
            continue
        chr = line.split('\t')[2]
        pos = line.split('\t')[3]
        if chr+'_'+pos in traits:
            out.write(line+'\t'+';;'.join(traits[chr+'_'+pos])+'\n')
        else:
            out.write(line+'\t-\n')

print('written to '+args.outfile)
