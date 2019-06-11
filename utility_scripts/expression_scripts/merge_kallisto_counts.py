import sys
import glob
import argparse
from multiprocessing import Pool


parser = argparse.ArgumentParser(description='Merge multiple kallisto count files into a matrix.')
parser.add_argument('kallisto_base_path', help='base path from where to search for kallisto abundance.tsv files')
parser.add_argument('transcript_to_gene_ID', help='File containing mapping of transcript ID to gene ID')
parser.add_argument('outfile', help='output file name')
parser.add_argument('--threads', help='Number of threads', default=1)


args = parser.parse_args()

transcript_to_gene = {}
with open(args.transcript_to_gene_ID) as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        transcript_to_gene[line[0]] = line[1]

def main():
    kallisto_files = glob.glob(args.kallisto_base_path+'/**/abundance.tsv', recursive=True)

    # dict to sum transcript coutns per gene per sample
    estimated_counts_per_gene = {}
    # loop over all kallisto files to get list of sample names
    for f in kallisto_files:
        sample = f.split('/')[-2]
        estimated_counts_per_gene[sample] = {}

    # get a set of genes (all files have the same list so only need to read 1 file)
    set_of_genes = set([])
    with open(kallisto_files[0]) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.split('\t')
            set_of_genes.add(transcript_to_gene[line[0]])

    # parallel process kallisto files. Per sample 1 thread
    with Pool(int(args.threads)) as pool:
        kallisto_data = pool.map(parse_kallisto_files, kallisto_files)
    
    print('Done reading kallisto file, start writing output matrix')
    sys.stdout.flush()

    with open(args.outfile,'w') as out:
        for result in kallisto_data:
            # result[0] is the sample name
            out.write('\t'+result[0])
        out.write('\n')
        for gene in set_of_genes:
            out.write(gene)
            for result in kallisto_data:
                # result[1] is the kallisto data for that specific sample
                out.write('\t'+str(result[1][gene]))
            out.write('\n')


def parse_kallisto_files(kallisto_abundance_file):
    # loop over all the kallisto abundance files. Take the estimated counts and 
    # sum counts from all transcript of certain gene to that gene
    sample_name = kallisto_abundance_file.split('/')[-2]
    print(kallisto_abundance_file)
    sys.stdout.flush()
    estimated_counts_per_gene = {}
    with open(kallisto_abundance_file) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            est_counts = float(line[3])
            transcript = line[0]
            gene = transcript_to_gene[transcript]
            if gene in estimated_counts_per_gene:
                estimated_counts_per_gene[gene] += est_counts
            else:
                estimated_counts_per_gene[gene] = est_counts
    return(sample_name, estimated_counts_per_gene)
            
if __name__ == '__main__':
    main()
