import argparse
import gzip
parser = argparse.ArgumentParser(description='Compare eQTLgen to MetaBrain')
parser.add_argument('eQTL_file',
                    help='File with eQTL data')
parser.add_argument('gwas_file',
                    help='File with gwas data')
parser.add_argument('gwas_annotation_file',
                    help='File with gwas annotation')
parser.add_argument('tested_snps',
                    help='File with tested SNPs')
parser.add_argument('outfile',
                    help='outfile name')
args = parser.parse_args()

print('start reading tested SNPs')
with gzip.open(args.tested_snps,'rt') as input_file:
    SNPs_tested = set(input_file.read().split('\n'))

trait_info = {}
print('start reading GWAS annotation file')
with open(args.gwas_annotation_file) as input_file:
    input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        trait_info[line[0].lower()] = [line[1], line[2]]

traits = {}
SNPs_in_gwas = set([])
print('start reading GWAS file')
snp_gwas_info = {}
with open(args.gwas_file) as input_file:
    header = input_file.readline().strip().split('\t')
    snps_index = header.index('SNPS')
    pval_index = header.index('P-VALUE')
    OR_or_beta_index = header.index('OR or BETA')
    risk_allele_index = header.index('STRONGEST SNP-RISK ALLELE')
    for line in input_file:
        line = line.strip().split('\t')
        trait = line[7]
        pval = line[pval_index]
        OR_or_beta = line[OR_or_beta_index]
        risk_allele = line[risk_allele_index]
        SNPs = line[snps_index].split(';')
        if len(SNPs) > 1:
            print(', '.join(SNPs))
            print(OR_or_beta)
            print(risk_allele)
            exit()
        for snp in SNPs:
            if 'rs' not in snp:
                continue
            if snp not in traits:
                traits[snp] = []
            traits[snp].append([trait, pval, OR_or_beta, risk_allele])
            SNPs_in_gwas.add(snp)
            snp_gwas_info

SNPs_tested_in_gwas = SNPs_tested & SNPs_in_gwas



eqtl_SNPs = {}
snp_seen = set([])
highest_pvalue = -1
print('start reading eQTL file '+args.eQTL_file)
print('start writing results to '+args.outfile)
with gzip.open(args.eQTL_file,'rt') as input_file, open(args.outfile,'wt') as out:
    out.write('SNP\tlowest_pval\ttrait\tparent_term\tmodified_parent_term\tgwas_pval\tgwas_OR_or_beta\tgwas_risk_allele\n')
    header = input_file.readline().strip('\n')
    split_head = header.split('\t')
    pval_index = split_head.index('PValue')
    snp_index = split_head.index('SNPName')
    x = 0
    y = 0
    for line in input_file:
        x += 1
        if x == 100000:
            print(x,'lines read,',y,'lines written')    
        line = line.strip('\n')
        snp = line.split('\t')[snp_index]
        pval = line.split('\t')[pval_index]
        highest_pvalue = pval
        if 'rs' in snp and ':' in snp:
            snp = snp.split(':')[2]
        if snp in snp_seen: # or snp not in SNPs_tested_in_gwas:
            continue
        snp_seen.add(snp)
        if snp not in traits:
            out.write(snp+'\t'+pval+'\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\n')
            continue
        traits_of_snp = traits[snp]
        for trait_gwas_info in traits_of_snp:
            trait = trait_gwas_info[0].lower()
            try:
                parent_term = trait_info[trait][0]
                modified_parent_term = trait_info[trait][1]
            except KeyError:
                print(trait,'not found')
                parent_term = 'UNKNOWN'
                modified_parent_term = 'UNKNOWN'
            out.write(snp+'\t'+pval+'\t'+trait_gwas_info[0]+'\t'+parent_term+'\t'+modified_parent_term+'\t'+trait_gwas_info[1]+'\t'+trait_gwas_info[2]+'\t'+trait_gwas_info[3]+'\n')
            y += 1
    print('following list of SNPs not in eqtl file, giving the highest p-value ('+str(highest_pvalue)+')')
    for snp in SNPs_tested_in_gwas - snp_seen:
        if snp not in traits:
            out.write(snp+'\t'+pval+'\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\tNotInGWAS\n')
            continue
        traits_of_snp = traits[snp]
        for trait_gwas_info in traits_of_snp:
            trait = trait_gwas_info[0].lower()
            try:
                parent_term = trait_info[trait][0]
                modified_parent_term = trait_info[trait][1]
            except KeyError:
                print(trait,'not found')
                parent_term = 'UNKNOWN'
                modified_parent_term = 'UNKNOWN'
            out.write(snp+'\t'+pval+'\t'+trait_gwas_info[0]+'\t'+parent_term+'\t'+modified_parent_term+'\t'+trait_gwas_info[1]+'\t'+trait_gwas_info[2]+'\t'+trait_gwas_info[3]+'\n')
        
print('written to '+args.outfile)
print('Done!')



