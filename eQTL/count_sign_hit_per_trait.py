import pandas as pd;
import argparse
import gzip
from multiprocessing import Process,Pool

parser = argparse.ArgumentParser(description='Compare eQTLgen to MetaBrain')
parser.add_argument('meataBrain_eQTL_file',
                    help='File with eQTL data')
parser.add_argument('eQTLgen_eQTL_file',
                    help='File with gwas data')
parser.add_argument('outfile',
                    help='outfile name')
parser.add_argument('SNPs_to_include_1',
                    help='file with SNPs to include')
parser.add_argument('SNPs_to_include_2',
                    help='File with SNPs to include')
args = parser.parse_args()

snps_to_include_1 = set([])
with gzip.open(args.SNPs_to_include_1,'rt') as input_file:
    for line in input_file:
        snp = line.strip()
        if snp.startswith('rs'):
            snps_to_include_1.add(snp)

snps_to_include_2 = set([])
with gzip.open(args.SNPs_to_include_2,'rt') as input_file:
    for line in input_file:
        snp = line.strip()
        if snp.startswith('rs'):
            snps_to_include_2.add(snp)

snps_to_use = snps_to_include_1 & snps_to_include_2
print(len(snps_to_use),'overlapping snps')

def parse_traits(eqtl_file, snps_to_use):
    counts_per_trait = {'total':0}
    print('start reading '+eqtl_file)
    with gzip.open(eqtl_file,'rt') as input_file:
        header = input_file.readline().strip().split('\t')
        gwas_index = header.index('GWAS_TRAIT')
        snp_name_index = header.index('SNPName')
        x = 0
        for line in input_file:
            x += 1
            if x % 1000000 == 0:
                print(eqtl_file+': '+str(x))
            line = line.strip().split('\t')
            snp = line[snp_name_index]
            if 'rs' not in snp:
                continue
            if ':' in snp:
                snp = snp.split(':')[2]
            gwas_trait = line[gwas_index]
            if gwas_trait == '-':
                continue
            for trait in gwas_trait.split(';;'):
                if trait == '-':
                    raise RuntimeError('should not happen')
                if trait not in counts_per_trait:
                    counts_per_trait[trait] = 0

            if snp not in snps_to_use:
                continue

            for trait in gwas_trait.split(';;'):
                if trait == '-':
                    continue
                counts_per_trait[trait] += 1
            counts_per_trait['total']  += 1
    return(counts_per_trait)


if __name__ == '__main__':
    pool = Pool(3)
    workers = [pool.apply_async(parse_traits, args=(args.meataBrain_eQTL_file,snps_to_use,)),
               pool.apply_async(parse_traits, args=(args.eQTLgen_eQTL_file,snps_to_use,))] 
    final_result = [worker.get() for worker in workers]
    metaBrain = final_result[0]
    eQTLgen = final_result[1]
    unique_traits = set([])
    for trait in metaBrain:
        unique_traits.add(trait)
    for trait in eQTLgen:
        unique_traits.add(trait)
    with open(args.outfile,'w') as out:
        out.write('Trait\tn_in_MetaBrain\tn_in_eQTLgen\n')
        for trait in unique_traits:
            out.write(trait+'\t')
            if trait in metaBrain:
                out.write(str(metaBrain[trait]))
            else:
                out.write('0')
            if trait in eQTLgen:
                out.write('\t'+str(eQTLgen[trait]))
            else:
                out.write('\t0')
            out.write('\n')

print('written to '+args.outfile)

