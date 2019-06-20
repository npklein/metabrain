import argparse

parser = argparse.ArgumentParser(description='Merge psychEncode with eQTLs from Sieberts et al.')
parser.add_argument('Sieberts_eqtls', help='eQTLs from Sieberts et al. study')
parser.add_argument('psychEncode_eqtls', help='eQTLs from psychEncode study')

args = parser.parse_args()


def parse_sieberts():
    n_in_Sieberts = 0
    sieberts_eQTLs = {}
    print('start reading Sieberts data')
    top_snp_per_gene = {}
    lowest_pval_gene = {}
    with open(args.Sieberts_eqtls) as input_file:
        for index, line in enumerate(input_file):
            if index == 10000:
                print(index)
            line = line.strip().split('\t')
            chr = line[0]
            pos = line[2]
            qtl_info = line[3].split(',')
            snp_name = qtl_info[2]
            gene = qtl_info[4]
            FDR = qtl_info[8]

        
            if gene not in top_snp_per_gene:
                top_snp_per_gene[gene] = chr+'_'+pos+'_'+gene
                lowest_pval_gene[gene] = float(FDR)
            elif float(FDR) < lowest_pval_gene[gene]:
                top_snp_per_gene[gene] = chr+'_'+pos+'_'+gene
                lowest_pval_gene[gene] = float(FDR)

            expr_increasing_allele = qtl_info[13]
            beta = qtl_info[9]
            allele1 = qtl_info[10]
            allele2 = qtl_info[11]
            if allele1 == expr_increasing_allele:
                if float(beta) > 0:
                    assessed_allele = allele1
                else:
                    assessed_allele = allele2
            elif allele2 == expr_increasing_allele:
                if float(beta) > 0:
                    assessed_allele  = allele2
                else:
                    assessed_allele = allele1
            else:
                raise RuntimeError('assessed alllele not in genotype')
            if chr+'_'+pos+'_'+gene in sieberts_eQTLs:
#                print(chr+'_'+pos+'_'+gene +' is a tri-allelic SNP, skip second line')
                continue
    #            raise RuntimeError('position_gene combination seen multiple times')
            sieberts_eQTLs[chr+'_'+pos+'_'+gene] = [snp_name, FDR, expr_increasing_allele, beta, assessed_allele, allele1+'/'+allele2]
            n_in_Sieberts += 1
    return(sieberts_eQTLs, n_in_Sieberts, top_snp_per_gene)

sieberts_eQTLs, n_in_Sieberts, top_snp_per_gene = parse_sieberts()


psychEncode_not_in_sieberts = 0
psychEncode_in_sieberts = 0
gene_seen = {}
print('done, start reading psychEncode and start writing')
with open(args.psychEncode_eqtls, 'rt') as input_file, open('2019-06-12.psychEncode-Sieberts-replication-FDR05.txt','w') as out:
    header = input_file.readline().strip().split('\t')
    out.write('chr\tpos\tSNP\tgene\tpsychEncode_FDR\tpsychEncode_zScore\tpsychEncode_zScore_swapped\tpsychEncode_alleleAssessed\tSieberts_FDR\tSieberts_beta\tSieberts_alleleAssessed\tpsychEncode_topSNP\tSieberts_topSNP\tsame_direction\n')
    for line in input_file:
        if len(line.strip()) == 0:
            continue
        line = line.strip().split('\t')
        FDR = line[-2]
        SNP = line[1]
        chr = line[2]
        pos = line[3]
        gene = line[4]
        if gene not in gene_seen:
            top_snp = 'True'
        else:
            top_snp = 'False'
        gene_seen[gene] = float(FDR)
        snpType = line[8]
        allele_assessed  = line[9]
        gene = line[4].split('.')[0]
        zscore = line[10]
        if chr+'_'+pos+'_'+gene in sieberts_eQTLs:
            if top_snp_per_gene[gene] == chr+'_'+pos+'_'+gene:
                top_snp_sieberts = 'True'
            else:
                top_snp_sieberts = 'False'
            sieberts = sieberts_eQTLs[chr+'_'+pos+'_'+gene]
            if SNP.startswith('rs') and SNP != sieberts[0]:
                print('rs IDs not the same, possibly merged?', SNP, sieberts[0])
            if snpType != sieberts[5] and sieberts[5].split('/')[1]+'/'+sieberts[5].split('/')[0] != snpType:
                print('Genotype not the same, skip: '+chr+'_'+pos+'_'+gene)
            
            if sieberts[4] != allele_assessed:
                zscore_swapped = -1*float(zscore)
            else:
                zscore_swapped = zscore
            same_direction = abs(float(zscore_swapped) + float(sieberts[3])) == abs(float(zscore_swapped)) + abs(float(sieberts[3]))
            out.write(chr+'\t'+pos+'\t'+SNP+'\t'+gene+'\t'+FDR+'\t'+zscore+'\t'+str(zscore_swapped)+'\t'+allele_assessed+'\t'+sieberts[1]+'\t'+sieberts[3]+'\t'+sieberts[4]+'\t'+top_snp+'\t'+top_snp_sieberts+'\t'+str(same_direction)+'\n')
            psychEncode_in_sieberts += 1
        else:
 #           print(chr+'_'+pos+'_'+gene, 'not in Sieberts')
            psychEncode_not_in_sieberts += 1

print('SNP-gene combis overlapping:',psychEncode_in_sieberts)
print('psychEncode not in Sieberts:',psychEncode_not_in_sieberts)
print('sieberts not in psychEncode:',n_in_Sieberts-psychEncode_in_sieberts)
