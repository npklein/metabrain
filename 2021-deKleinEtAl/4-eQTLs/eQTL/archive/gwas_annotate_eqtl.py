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
    header = input_file.readline().strip().split('\t')
    snps_index = header.index('SNPS')
    for line in input_file:
        line = line.strip().split('\t')
        trait = line[7]
        for snp in line[snps_index].split('; '):
            if 'rs' not in snp:
                continue
            if snp not in traits:
                traits[snp] = set([])
            traits[snp].add(trait)

with gzip.open(args.eQTL_file,'rt') as input_file, gzip.open(args.outfile,'wt') as out:
    header = input_file.readline().strip('\n')
    split_head = header.split('\t')
    fdr_index = split_head.index('FDR')
    snp_index = split_head.index('SNPName')
    out.write(header+'\tGWAS_TRAIT\n')
    for line in input_file:
        line = line.strip('\n')
        if float(line.split('\t')[fdr_index]) >= 0.05:
            continue
        snp = line.split('\t')[snp_index]
        if 'rs' in snp and ':' in snp:
            snp = snp.split(':')[2]
        if snp in traits:
            out.write(line+'\t'+';;'.join(traits[snp])+'\n')
#        else:
#            out.write(line+'\t-\n')

print('written to '+args.outfile)
