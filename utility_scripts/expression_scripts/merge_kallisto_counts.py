import sys
import glob
import argparse
from multiprocessing import Pool
import gzip

parser = argparse.ArgumentParser(description='Merge multiple kallisto count files into a matrix.')
parser.add_argument('kallisto_base_path', help='base path from where to search for kallisto abundance.tsv files')
parser.add_argument('transcript_to_gene_ID', help='File containing mapping of transcript ID to gene ID')
parser.add_argument('outfile_geneCounts', help='output file name gene counts')
parser.add_argument('outfile_transcriptTPMs', help='output file name transcript TPMs')
parser.add_argument('--threads', help='Number of threads', default=1)


args = parser.parse_args()

transcript_to_gene = {}
with open(args.transcript_to_gene_ID) as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        transcript_to_gene[line[0]] = line[1]

def openfile(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def main():
    kallisto_files = glob.glob(args.kallisto_base_path+'/**/abundance.tsv', recursive=True)

    # dict to sum transcript coutns per gene per sample
    estimated_counts_per_gene = {}
    print('loop over all kallisto files to get list of sample names')
    sys.stdout.flush()
    for f in kallisto_files:
        sample = f.split('/')[-2]
        estimated_counts_per_gene[sample] = {}
    print('done')
    sys.stdout.flush()
    # get a set of genes (all files have the same list so only need to read 1 file)
    set_of_genes = set([])
    set_of_transcripts = set([])
    with open(kallisto_files[0]) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.split('\t')
            set_of_genes.add(transcript_to_gene[line[0]])
            set_of_transcripts.add(line[0])
    print('start map')
    sys.stdout.flush()
    # parallel process kallisto files. Per sample 1 thread
    print('starting',args.threads,'processes')
    with Pool(int(args.threads)) as pool:
        kallisto_data = pool.map(parse_kallisto_files, kallisto_files)
    
    print('Done reading kallisto file, start writing output matrix')
    sys.stdout.flush()

    with openfile(args.outfile_geneCounts,'wt') as out_geneCounts:
        for result in kallisto_data:
            # result[0] is the sample name
            out_geneCounts.write('\t'+result[0])
        out_geneCounts.write('\n')
        for gene in set_of_genes:
            out_geneCounts.write(gene)
            for result in kallisto_data:
                # result[1] is the kallisto genecount data for that specific sample
                out_geneCounts.write('\t'+str(result[1][gene]))
            out_geneCounts.write('\n')

    with openfile(args.outfile_transcriptTPMs,'wt') as out_TPM:
        for result in kallisto_data:
            # result[0] is the sample name
            out_TPM.write('\t'+result[0])
        out_TPM.write('\n')
        for transcript in set_of_transcripts:
            out_TPM.write(transcript)
            for result in kallisto_data:
                # result[2] is the kallisto transcript TPM data for that specific sample
                out_TPM.write('\t'+str(result[2][transcript]))
            out_TPM.write('\n')



def parse_kallisto_files(kallisto_abundance_file):
    # loop over all the kallisto abundance files. Take the estimated counts and 
    # sum counts from all transcript of certain gene to that gene
    sample_name = kallisto_abundance_file.split('/')[-2]
    estimated_counts_per_gene = {}
    tpm_per_transcript = {}
    with open(kallisto_abundance_file) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            est_counts = float(line[3])
            tpm = line[4]
            transcript = line[0]
            gene = transcript_to_gene[transcript]
            if gene in estimated_counts_per_gene:
                estimated_counts_per_gene[gene] += est_counts
            else:
                estimated_counts_per_gene[gene] = est_counts
            tpm_per_transcript[transcript] = tpm
    return(sample_name, estimated_counts_per_gene, tpm_per_transcript)
            
if __name__ == '__main__':
    main()
