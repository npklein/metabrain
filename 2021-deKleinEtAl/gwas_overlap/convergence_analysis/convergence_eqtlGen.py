import glob
import gzip

def main():
    # TODO: change this to command line input
    for f in glob.glob('/groups/umcg-biogen/tmp03/umcg-ndeklein/gwas_overlap/biogen-cisandtrans-4neurodiseases/eqtlgen_disease_overlap/*trans*.txt.gz'):
        # filename has to be in format *-<trait>.txt
        trait = f.split('-')[-1].split('.txt')[0]
        cis_f = f.replace('_trans-','_cis-')
        out_file = '/groups/umcg-biogen/tmp03/umcg-ndeklein/gwas_overlap/biogen-cisandtrans-4neurodiseases/eqtlGen_convergent_genes/convergent_genes_'+trait+'.txt'
        get_convergent_genes(cis_f, f, out_file)

    for f in glob.glob('/groups/umcg-biogen/tmp03/umcg-ndeklein/gwas_overlap/biogen-cisandtrans-4neurodiseases/10pcs-trans/*txt.gz'):
        if '-' not in f.split('/')[-1]:
            continue
        trait = f.split('-')[-1].split('.txt')[0]
        cis_f = f.replace('-trans/','-cis/').replace('FDR0.05-','FDR0.05-ProbeLevel-')
        out_file = '/groups/umcg-biogen/tmp03/umcg-ndeklein/gwas_overlap/biogen-cisandtrans-4neurodiseases/MetaBrain_convergent_genes/convergent_genes_'+trait+'.txt'
        get_convergent_genes(cis_f, f, out_file)


def parse_eqtl_file(f):
    '''Get genes + their SNPs from eQTL file'''
    gene_snps = {}
    snp_location_per_chr = {}
    with gzip.open(f) as input_file:
        snp_include = set([])
        header = input_file.readline()
        for line in input_file:
            line = line.decode('utf-8').split('\t')
            gene = line[4]
            chr = line[2]

            if chr not in snp_location_per_chr:
                snp_location_per_chr[chr] = []
            pos = int(line[3])

            if gene not in gene_snps:
                gene_snps[gene] = set([])
            snp = line[1]
            
            # check that there is not other SNP from same loci (5Mb window)
            if snp not in snp_include:
                for other_pos in snp_location_per_chr[chr]:
                    if abs(other_pos-pos) < 5000000:
                        snp_location_per_chr[chr].append(pos)
                        continue
            snp_include.add(snp)
            snp_location_per_chr[chr].append(pos)
            gene_snps[gene].add(snp)
    return(gene_snps)

def get_convergent_genes(cis_file, trans_file, out_file):
    '''From a cis and trans file get convergent genes (genes with multiple eQTLs'''
    with open(out_file,'w') as out:
        out.write('gene\tcis_snps\ttrans_snps\n')

        cis_gene_snps = parse_eqtl_file(cis_file)
        trans_gene_snps = parse_eqtl_file(trans_file)
        # Loop over all unique genes from both the cis and trans eQTL file (keys() gets the keys, convert both to set, union is those that are present in both)
        for gene in set(cis_gene_snps.keys()).union(set(trans_gene_snps.keys())):
            n_snp_cis = 0
            cis_snps = ''
            if gene in cis_gene_snps:
                cis_snps = cis_gene_snps[gene]
                n_snp_cis = len(cis_snps)
            n_snp_trans = 0
            trans_snps = ''
            if gene in trans_gene_snps:
                trans_snps = trans_gene_snps[gene]
                n_snp_trans = len(trans_snps)


            # Check if cis + trans combined more than 1 snp affects some gene
            if n_snp_cis + n_snp_trans > 1:
                   out.write(gene+'\t'+','.join(cis_snps)+'\t'+','.join(trans_snps)+'\n')



if __name__ == "__main__":
    main()

