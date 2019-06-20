import argparse
import gzip
parser = argparse.ArgumentParser(description='Merge MetaBrain with eQTLs from Sieberts et al.')
parser.add_argument('Sieberts_eqtls', help='eQTLs from Sieberts et al. study')
parser.add_argument('name', help='name that is given in output columns')
parser.add_argument('MetaBrain_eqtls', help='eQTLs from MetaBrain study')
parser.add_argument('output_file', help='output file name')

args = parser.parse_args()

n_in_Sieberts = 0
sieberts_eQTLs = {}
print('start reading '+args.name+' data')
top_snp_per_gene = {}
lowest_pval_gene = {}
with open(args.Sieberts_eqtls) as input_file:
    for index, line in enumerate(input_file):
        if index == 10000:
            print(args.Sieberts_eqtls,'-',index)
        line = line.strip().split('\t')
        chr = line[0].replace('chr','')
        pos = line[2]
        qtl_info = line[3].split(',')
        snp_name = qtl_info[2]
        if len(qtl_info) == 17:
            index_modified = 1
        elif len(qtl_info) == 18:
            index_modified = 0
        else:
            print(args.Sieberts_eqtls)
            print(qtl_info)
            print(len(qtl_info))
            raise RuntimeError('Unexpected length')

        #chromosome,snpLocation,snpid,         gene,geneSymbol,statistic,pvalue,FDR,beta,A1,A2,A2freq,expressionIncreasingAllele,strand,geneBiotype,geneStartPosition,geneEndPosition
        #chromosome,snpLocation,snpid,snpLocId,gene,geneSymbol,statistic,pvalue,FDR,beta,A1,A2,A2freq,expressionIncreasingAllele,strand,geneBiotype,geneStartPosition,geneEndPosition
        gene = qtl_info[4-index_modified]
        FDR = qtl_info[8-index_modified]
        expr_increasing_allele = qtl_info[13-index_modified]
        beta = qtl_info[9-index_modified]
        allele1 = qtl_info[10-index_modified]
        allele2 = qtl_info[11-index_modified]

        if gene not in top_snp_per_gene:
            top_snp_per_gene[gene] = chr+'_'+pos+'_'+gene
            lowest_pval_gene[gene] = float(FDR)
        elif float(FDR) < lowest_pval_gene[gene]:
            top_snp_per_gene[gene] = chr+'_'+pos+'_'+gene
            lowest_pval_gene[gene] = float(FDR)
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
            print(allele1, allele2, expr_increasing_allele)
            raise RuntimeError('assessed alllele not in genotype')
        if chr+'_'+pos+'_'+gene in sieberts_eQTLs:
#            print(chr+'_'+pos+'_'+gene +' is a tri-allelic SNP, skip second line')
            continue
        sieberts_eQTLs[chr+'_'+pos+'_'+gene] = [snp_name, FDR, expr_increasing_allele, beta, assessed_allele, allele1+'/'+allele2]
        n_in_Sieberts += 1

metabrain_not_in_sieberts = 0
metabrain_in_sieberts = 0
gene_seen = {}
print(args.Sieberts_eqtls,'- done, start reading metabrain and start writing')
with gzip.open(args.MetaBrain_eqtls, 'rt') as input_file, open(args.output_file,'w') as out:
    header = input_file.readline().strip().split('\t')
    out.write('chr\tpos\tSNP\tgene\tmetaBrain_FDR\tmetaBrain_zScore\tmetaBrain_zScore_swapped\tmetaBrain_alleleAssessed')
    out.write('\t'+args.name+'_FDR\t'+args.name+'_beta\t'+args.name+'_beta_swapped\t'+args.name+'_alleleAssessed\tmetaBrain_topSNP\tSieberts_topSNP\tsame_direction\n')
    for index, line in enumerate(input_file):
        if index == 10000:
            print(args.Sieberts_eqtls,'-',index)
        line = line.strip().split('\t')
        FDR = line[-1]
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
#            if SNP.startswith('rs') and SNP != sieberts[0]:
#                print('rs IDs not the same, possibly merged?', SNP, sieberts[0])
            if snpType != sieberts[5] and sieberts[5].split('/')[1]+'/'+sieberts[5].split('/')[0] != snpType:
#                print('Genotype not the same, skip: '+chr+'_'+pos+'_'+gene)
                continue
            if sieberts[4] != allele_assessed:
                zscore_swapped = -1*float(zscore)
                beta_swapped = -1 * float(sieberts[3])
            else:
                zscore_swapped = zscore
                beta_swapped = sieberts[3]
            same_direction = abs(float(zscore_swapped) + float(sieberts[3])) == abs(float(zscore_swapped)) + abs(float(sieberts[3]))
            out.write(chr+'\t'+pos+'\t'+SNP+'\t'+gene+'\t'+FDR+'\t'+zscore+'\t'+str(zscore_swapped)+'\t'+allele_assessed+'\t'+sieberts[1]+'\t')
            out.write(sieberts[3]+'\t'+str(beta_swapped)+'\t'+sieberts[4]+'\t'+top_snp+'\t'+top_snp_sieberts+'\t'+str(same_direction)+'\n')
            metabrain_in_sieberts += 1
        else:
#            print(chr+'_'+pos+'_'+gene, 'not in Sieberts')
            metabrain_not_in_sieberts += 1

with open(args.output_file.replace('.txt','.overlapInfo.txt'),'w') as out:
    out.write('SNP-gene combis overlapping: '+str(metabrain_in_sieberts)+'\n')
    out.write('Metabrain not in Sieberts: '+str(metabrain_not_in_sieberts)+'\n')
    out.write('sieberts not in metabrain: '+str(n_in_Sieberts-metabrain_in_sieberts)+'\n')
