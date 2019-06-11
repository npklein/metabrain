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
snps_to_use = snps_to_include_1.intersection(snps_to_include_2)

def parse_traits(eqtl_file, snps_to_use):
    counts_per_trait = {'total':0}
    print(eqtl_file+' - '+str(len(snps_to_use))+' overlapping snps')
    print('start reading '+eqtl_file)
    snps_no_gwas = 0
    snps_no_rs = 0
    snps_not_overlapping = 0
    snps_double = 0
    snps_seen = set([])
    snps_seen_all = set([])
    with gzip.open(eqtl_file,'rt') as input_file:
        header = input_file.readline().strip().split('\t')
        gwas_index = header.index('GWAS_TRAIT')
        snp_name_index = header.index('SNPName')
        x = 0
        for line in input_file:
            x += 1
            line = line.strip().split('\t')
            snp = line[snp_name_index]
            if 'rs' in snp and ':' in snp:
                snp = snp.split(':')[2]

            if snp not in snps_to_use:
                snps_not_overlapping += 1
                continue

            snps_seen_all.add(snp)
            gwas_trait = line[gwas_index]

            if gwas_trait == '-':
                snps_no_gwas += 1
                continue

            if snp in snps_seen:
                snps_double += 1
                continue
            snps_seen.add(snp)

            for trait in gwas_trait.split(';;'):
                if trait not in counts_per_trait:
                    counts_per_trait[trait] = 0
                counts_per_trait[trait] += 1
            counts_per_trait['total']  += 1
    print(eqtl_file+' - total: '+str(x))
    print(eqtl_file+' - SNPs not overlapping: '+str(snps_not_overlapping))
    print(eqtl_file +' SNPs FDR > 0.05:'+str(x-len(snps_seen_all)))
    print(eqtl_file+' - SNPs no GWAS: '+str(snps_no_gwas))
    print(eqtl_file+' - SNPs double: '+str(snps_not_overlapping))
    print(eqtl_file+' - total: '+str(counts_per_trait['total']))
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

