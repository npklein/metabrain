import gzip
import argparse
import os

parser = argparse.ArgumentParser(description='Remove genes from expression file that have exact same expression levels (likely duplicates), are only found on scaffolds, or have duplicate names.')
parser.add_argument('-g','--gtf', help='GTF file')
parser.add_argument('-e','--expression_file', help='expression file from which to filter genes')
parser.add_argument('-o','--outfile', help='Outfile to write')

args = parser.parse_args()


outdir = os.path.dirname(args.outfile)
print('write output to '+outdir)
if not os.path.exists(outdir):
    print('makedir '+outdir)
    os.makedirs(outdir)

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)

def parse_gtf():
    gene_chr = {}
    with openfile(args.gtf) as input_file:
        for line in input_file:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if line[2] != 'gene':
                continue
            gene_id = line[8].split('gene_id "')[1].split('"')[0]
            chr = line[0]
            if gene_id not in gene_chr:
                gene_chr[gene_id] = set()
            gene_chr[gene_id].add(chr)
    return(gene_chr)
gene_chr = parse_gtf()

accepted_chr = set(['1','10','11','12','13','14','15','16','17','18','19','2','20','21','22','3','4','5','6','7','8','9','MT','X','Y'])



expr_set = {}
genes_to_filter = set([])
print('start first read of file to find duplicate expression')
gene_on_scaffold = 0
genes_to_keep = set([])
no_expression = 0
gene_seen = set([])
average_expression_10 = 0
with openfile(args.expression_file,'rt') as input_file:
    header = input_file.readline()
    for line in input_file:
        spl_line = line.split('\t')
        gene = spl_line[0]

        if gene.startswith('LRG'):
            gene_on_scaffold += 1
            continue
        all_on_scaffolds = True
        for chr in gene_chr[gene]:
            if chr in accepted_chr:
                all_on_scaffolds = False
        if all_on_scaffolds:
            gene_on_scaffold += 1
            continue

        if gene in gene_seen:
            continue
        gene_seen.add(gene)

        all_expression = '_'.join(spl_line[1:])
        if all_expression not in expr_set:
            expr_set[all_expression] = set(gene)
        else:
            genes_to_filter.add(gene)

        has_expr = False
#        average_expression = 0
        for element in spl_line[1:]:
            if float(element) > 0:
                has_expr = True
#                average_expression += float(element)
                break
#        average_expression = average_expression / len(spl_line[1:])
        if not has_expr:
            no_expression += 1
            continue
#        if average_expression < 10:
#            average_expression_10 += 1
#            continue

        genes_to_keep.add(gene)

gene_seen = set([])
print('Done. Start second read of file to write filtered expression file')
duplicate_genes = 0
genes_written = 0
with openfile(args.expression_file,'rt') as input_file, openfile(args.outfile,'wt') as out:
    out.write(header)
    for line in input_file:
        gene = line.split('\t')[0]
        if gene not in genes_to_keep:
            continue

        if gene in gene_seen:
            duplicate_genes += 1
            continue
        gene_seen.add(gene)

        if gene not in genes_to_filter:
            out.write(line)
            genes_written += 1

with open(outdir+'/number_of_genes_filtered.txt','w') as out:
    out.write('Filtered genes for these reasons:\n')
    out.write('Gene on scaffold: '+str(gene_on_scaffold)+'\n')
    out.write('Genes with duplicate expression: '+str(len(genes_to_filter))+'\n')
    out.write('Duplicate gene IDs: '+str(duplicate_genes)+'\n')
    out.write('All samples 0 reads: '+str(no_expression)+'\n')
#    out.write('Average expression < 10: '+str(average_expression_10)+'\n')
    out.write('-'*20+'\n')
    out.write('Genes kept: '+str(genes_written)+'\n')

