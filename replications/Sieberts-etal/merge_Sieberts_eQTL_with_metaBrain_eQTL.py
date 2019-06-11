import argparse
import gzip
parser = argparse.ArgumentParser(description='Merge MetaBrain with eQTLs from Sieberts et al.')
parser.add_argument('Sieberts_eqtls', help='eQTLs from Sieberts et al. study')
parser.add_argument('MetaBrain_eqtls', help='eQTLs from MetaBrain study')

args = parser.parse_args()

sieberts_eQTLs = {}
with open(args.Sieberts_eqtls) as input_file:
    header = input_file.readline().strip().split('\t')
    print(header)
    for line in input_file:
        line = line.strip().split('\t')
        chr = line[0]
        pos = line[1]
        snp_name = line[2]
        gene = line[4]
        FDR = line[8]
        expr_increasin_allele = line[13]
        beta = line[9]
        print(line)
        break

with gzip.open(args.MetaBrain_eqtls, 'rt') as input_file:
    header = input_file.readline().strip().split('\t')
    print(header)
    for line in input_file:
        line = line.strip().split('\t')
        print(line)
        exit()
