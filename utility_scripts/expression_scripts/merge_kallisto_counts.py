import os
import sys
import glob
import argparse
from multiprocessing import Pool
import gzip
from pathlib import Path
import sys

parser = argparse.ArgumentParser(description='Merge multiple kallisto count files into a matrix.')
parser.add_argument('kallisto_base_path', help='base path from where to search for kallisto abundance.tsv files')
parser.add_argument('gtf', help='GTF file containing mapping of transcript ID to gene ID')
parser.add_argument('outfile_geneCounts', help='output file name gene counts')
parser.add_argument('outfile_transcriptTPMs', help='output file name transcript TPMs')
parser.add_argument('outfile_transcriptCounts', help='output file name transcript counts')
parser.add_argument('--threads', help='Number of threads', default=1)


args = parser.parse_args()

geneCounts_dir = os.path.dirname(args.outfile_geneCounts)
if len(geneCounts_dir) > 0:
    Path(geneCounts_dir).mkdir(parents=True, exist_ok=True)

transcriptTPM_dir = os.path.dirname(args.outfile_transcriptTPMs)
if len(transcriptTPM_dir) > 0:
    Path(transcriptTPM_dir).mkdir(parents=True, exist_ok=True)

transcriptCounts_dir = os.path.dirname(args.outfile_transcriptCounts)
if len(transcriptCounts_dir) > 0:
    Path(transcriptCounts_dir).mkdir(parents=True, exist_ok=True)


if not os.path.isdir(args.kallisto_base_path):
    raise RuntimeError(args.kallisto_base_path+' does not exist')


transcript_to_gene = {}
print('read gtf file', flush=True)

with open(args.gtf) as input_file:
    for line in input_file:
        if 'transcript_id' in line:
            info = line.split('\t')[8]
            transcript_id = info.split('transcript_id "')[1].split('"')[0].split('.')[0]
            gene_id = info.split('gene_id "')[1].split('"')[0]
            transcript_to_gene[transcript_id] = gene_id
print('done', flush=True)

def openfile(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def main():
    print('start search for abundance.tsv in '+ args.kallisto_base_path, flush=True)
    kallisto_files = glob.glob(args.kallisto_base_path+'/**/abundance.tsv', recursive=True)

    print('found',len(kallisto_files),' kallisto files', flush=True)
    
    # dict to sum transcript coutns per gene per sample
    estimated_counts_per_gene = {}
    print('loop over all kallisto files to get list of sample names', flush=True)
    sys.stdout.flush()
    for f in kallisto_files:
        sample = f.split('/')[-2]
        estimated_counts_per_gene[sample] = {}
    print('done', flush=True)
    sys.stdout.flush()
    # get a set of genes (all files have the same list so only need to read 1 file)
    set_of_genes = set([])
    set_of_transcripts = set([])
    with open(kallisto_files[0]) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.split('\t')
            set_of_genes.add(transcript_to_gene[line[0].split('.')[0]])
            set_of_transcripts.add(line[0].split('.')[0])
    print('start map', flush=True)
    sys.stdout.flush()
    # parallel process kallisto files. Per sample 1 thread
    print('starting',args.threads,'processes', flush=True)
    with Pool(int(args.threads)) as pool:
        kallisto_data = pool.map(parse_kallisto_files, kallisto_files)
    
    print('Done reading kallisto file, start writing output matrix', flush=True)
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

    with openfile(args.outfile_transcriptCounts,'wt') as out_counts:
        for result in kallisto_data:
            # result[0] is the sample name
            out_counts.write('\t'+result[0])
        out_counts.write('\n')
        for transcript in set_of_transcripts:
            out_counts.write(transcript)
            for result in kallisto_data:
                # result[3] is the kallisto transcript count data for that specific sample
                out_counts.write('\t'+str(result[3][transcript]))
            out_counts.write('\n')

def parse_kallisto_files(kallisto_abundance_file):
    # loop over all the kallisto abundance files. Take the estimated counts and 
    # sum counts from all transcript of certain gene to that gene
    sample_name = kallisto_abundance_file.split('/')[-2]
    estimated_counts_per_gene = {}
    tpm_per_transcript = {}
    count_per_transcript = {}
    with open(kallisto_abundance_file) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            try:
                est_counts = float(line[3])
                tpm = line[4]
                transcript = line[0].split('.')[0]
            except IndexError:
                print('ERROR!! Line does not have enough value for: '+kallisto_abundance_file, flush=True)
                print(line, flush=True)
                raise
            gene = transcript_to_gene[transcript]
            if gene in estimated_counts_per_gene:
                estimated_counts_per_gene[gene] += est_counts
            else:
                estimated_counts_per_gene[gene] = est_counts
            tpm_per_transcript[transcript] = tpm
            count_per_transcript[transcript] = est_counts
    return(sample_name, estimated_counts_per_gene, tpm_per_transcript, count_per_transcript)
            
if __name__ == '__main__':
    main()
