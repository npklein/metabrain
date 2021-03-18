# Merge FeatureCount count files (counts per sample) to matrix with rows genes, columns samples
# NOTE: Needs Python version 3.5+
import glob
import argparse
import gzip
import sys
parser = argparse.ArgumentParser(description='Calculate the fraction each transcript of a gene is expressed')
parser.add_argument('feature_count_matrix', help='Feature count matrix with as rows samples and columns genes')
parser.add_argument('gtf', help='GTF file')
parser.add_argument('out_file', help='outfile')

args = parser.parse_args()
feature_gene_id = {}
feature_per_gene = {}
with open(args.gtf) as input_file:
    for line in input_file:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        feature = line[2]
        if feature != 'transcript':
            continue
        chr = line[0]
        start = line[3]
        end = line[4]
        info = line[8]
        if 'transcript_id "' not in info:
            print(line)
            raise RuntimeError('transcript_id not in info')
        id = info.split('transcript_id "')[1].split('"')[0]
        gene_id = info.split('gene_id "')[1].split('"')[0]
        feature_gene_id[id] = gene_id
        if gene_id not in feature_per_gene:
            feature_per_gene[gene_id] = 0
        feature_per_gene[gene_id] += 1
            
transcript_indices = {}
with open(args.feature_count_matrix) as input_file:
    header = input_file.readline().strip().split('\t')
    for index, element in enumerate(header[1:]):
        transcript_indices[index] = element
    counts_per_gene = {}
    print('first loop through file')
    for line in input_file:
        line = line.strip().split('\t')
        for index, element in enumerate(line[1:]):
            gene_id = feature_gene_id[transcript_indices[index]]
            if gene_id not in counts_per_gene:
                counts_per_gene[gene_id] = 0
            counts_per_gene[gene_id] += float(element)

   
with open(args.feature_count_matrix) as input_file, open(args.feature_count_matrix.replace('.txt','.RATIOS.txt','w')) as out:
    header = input_file.readline()
    out.write(header)
    print('second loop through file')
    for line in input_file:
        line = line.strip().split('\t')
        out.write(line[0])
        for index, element in enumerate(line[1:]):
            gene_id = feature_gene_id[transcript_indices[index]]
            new_count = float(element) / counts_per_gene[gene_id]
            out.write('\t'+str(new_count))
        out.write('\n')
