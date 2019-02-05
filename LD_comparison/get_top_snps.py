import gzip
import argparse

parser = argparse.ArgumentParser(description='Get top SNPs for metabrain and eqtlgen')
parser.add_argument('metaBrain', help='File with the metaBrain results')
parser.add_argument('eqtlGen', help='File with the eqtlGen results')

args = parser.parse_args()

def parse_line(line):
    '''Decode binary line and select snp, gene and fdr'''
    line = line.decode('utf-8').strip().split('\t')
    fdr = float(line[-1])
    snp = line[1]
    # because eqtlGen does not have version number, make sure that it is only the ensID
    gene = line[4]
    return(snp, gene, fdr)

top_snp_per_gene_metaBrain = {}
top_snp_fdr = {}
print('parse metaBrain')
with gzip.open(args.metaBrain) as input_file:
    header = input_file.readline().decode('utf-8').strip().split('\t')
    if header[-1] != 'FDR':
        raise RuntimeError('Header not in expected order')
    x = 0
    for line in input_file:
        x += 1
        if x % 100000 == 0:
            print(x)
        snp, gene, fdr = parse_line(line)
        # because eqtlGen does not have version number, make sure that it is only the ensID
        gene = gene.split('.')[0]

        # Only interested in those that are significant for this analysis
        if fdr > 0.05:
            break

        # only want the top SNP. Since can't be sure that file is sorted by pvalues, check if pval is lower than previous
        # first occurence of the gene is also the top SNP. Continue if gene already seen
        if gene not in top_snp_per_gene_metaBrain or fdr < top_snp_fdr[gene]:
            top_snp_per_gene_metaBrain[gene] = snp
            top_snp_fdr[gene] = fdr


top_snp_per_gene_eqtlGen = {}
top_snp_fdr = {}
print('parsing eqtlGen')
with gzip.open(args.eqtlGen) as input_file:
    header = input_file.readline().decode('utf-8').strip().split('\t')
    if header[-1] != 'FDR':
        raise RuntimeError('Header not in expected order')

    pref_fdr = 0
    x = 0
    for line in input_file:
        x += 1
        if x % 100000 == 0:
            print(x)

        snp, gene, fdr = parse_line(line)

        if fdr > 0.05:
            break

        # only use genes that are also found in metaBrain
        if gene not in top_snp_per_gene_metaBrain:
            continue

        # then select topSNPs
        if gene not in top_snp_per_gene_eqtlGen or fdr < top_snp_fdr[gene]:
            top_snp_per_gene_eqtlGen[gene] = snp
            top_snp_fdr[gene] = fdr


with open('top_snp_metaBrain_and_eQTLgen.txt','w') as out:
    out.write('gene\ttopSNP_metaBrain\ttopSNP_eqtlGen\n')
    for gene in top_snp_per_gene_metaBrain:
        if gene in top_snp_per_gene_eqtlGen:
            out.write(gene+'\t'+top_snp_per_gene_metaBrain[gene]+'\t'+top_snp_per_gene_eqtlGen[gene]+'\n')
