raise RuntimeError('This script is wrong, need to take in consideration that the genotypes can be swapped')
import argparse
import gzip
parser = argparse.ArgumentParser(description='Merge MetaBrain with eQTLs from eqtlGen et al.')
parser.add_argument('eqtlGen_eqtls', help='eQTLs from eqtlGen')
parser.add_argument('name', help='name that is given in output columns')
parser.add_argument('MetaBrain_eqtls', help='eQTLs from MetaBrain study')
parser.add_argument('output_file', help='output file name')

args = parser.parse_args()

n_in_eqtlGen = 0
eqtlGen_eQTLs = {}
print('start reading '+args.name+' data')
top_snp_per_gene = {}
lowest_pval_gene = {}

with gzip.open(args.eqtlGen_eqtls,'rt') as input_file:
    header = input_file.readline().strip().split('\t')
    for index, line in enumerate(input_file):
        if index == 10000:
            print(args.eqtlGen_eqtls,'-',index)
        line = line.strip().split('\t')

        FDR = line[-1]
        rs = line[1]
        chr = line[2]
        pos = line[3]
        gene = line[4].split('.')[0]
        isTopSNP = False
        if gene not in lowest_pval_gene:
            lowest_pval_gene[gene] = float(FDR)
            isTopSNP = True
        elif lowest_pval_gene[gene] < float(FDR):
            raise RuntimeError('not FDR sorted')
        alleles = line[8].split('/')
        allele_assessed  = line[9]
        zscore = line[10]
        eqtlGen_eQTLs[rs+'_'+gene] = [zscore, allele_assessed, FDR, alleles, isTopSNP]
        n_in_eqtlGen += 1
lowest_pval_gene = {}
print(args.eqtlGen_eqtls,'- done, start reading metabrain and start writing')
with gzip.open(args.MetaBrain_eqtls, 'rt') as input_file, open(args.output_file,'w') as out:
    header = input_file.readline().strip().split('\t')
    out.write('chr\tpos\tSNP\tgene\tmetaBrain_FDR\tmetaBrain_zScore\tmetaBrain_zScore_swapped\tmetaBrain_alleleAssessed')
    out.write('\t'+args.name+'_FDR\t'+args.name+'_zscore\t'+args.name+'_alleleAssessed\tsame_direction\tis_metabrain_topSNP\tis_eqtlGen_topSNP\n')
    for index, line in enumerate(input_file):
        if index == 10000:
            print(args.eqtlGen_eqtls,'-',index)
        line = line.strip().split('\t')
        FDR = line[-1]
        SNP = line[1]
        rs = SNP.split(':')[2]
        chr = line[2]
        pos = line[3]
        gene = line[4]
        snpType = line[8]
        alleles = snpType.split('/')
        allele_assessed  = line[9]
        gene = line[4].split('.')[0]
        zscore = line[10]

        isTopSNP = False
        if gene not in lowest_pval_gene:
            lowest_pval_gene[gene] = float(FDR)
            isTopSNP = True
        elif lowest_pval_gene[gene] < float(FDR):
            raise RuntimeError('not FDR sorted')

        if rs+'_'+gene not in eqtlGen_eQTLs:
            print(rs+'_'+gene+' not in')
            #raise RuntimeError(rs+'_'+gene+' not in')
            continue
        eqtlGen = eqtlGen_eQTLs[rs+'_'+gene]
        if alleles != eqtlGen[3] and [alleles[1], alleles[0]] != eqtlGen[3]:
           print('/'.join(alleles)+' '+ '/'.join(eqtlGen[3])+'. Genotype not the same, skip: '+SNP+'_'+gene)
           continue
        if eqtlGen[1] != allele_assessed:
            zscore_swapped = -1*float(zscore)
        else:
            zscore_swapped = zscore
        same_direction = float(zscore_swapped) * float(eqtlGen[0]) >= 0
        out.write(chr+'\t'+pos+'\t'+SNP+'\t'+gene+'\t'+FDR+'\t'+zscore+'\t'+str(zscore_swapped)+'\t'+allele_assessed+'\t'+eqtlGen[2]+'\t')
        out.write(eqtlGen[0]+'\t'+eqtlGen[1]+'\t'+str(same_direction)+'\t'+str(isTopSNP)+'\t'+str(eqtlGen[4])+'\n')

