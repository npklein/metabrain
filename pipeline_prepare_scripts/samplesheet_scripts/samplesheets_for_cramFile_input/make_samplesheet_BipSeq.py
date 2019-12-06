
# loop through all subdirs to find all fastq files
import os.path
import glob
import os
import re
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for BipSeq.')
parser.add_argument('samplesheet', help='BipSeq samplesheet from synapse')
parser.add_argument('fastq_dir', help='path to fastq file dir')
parser.add_argument('--outdir',help='Directory where output is written',default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

individual_per_sample = {}
samples_per_batch = {}
batch_count = {}

# BipSeq is small enough to not need batches
with open(args.samplesheet) as input_file, open(args.outdir+'samplesheet_BipSeq_RNA.all_jobs.txt','w') as out:
    header = input_file.readline().split(',')
    out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
    for line in input_file:
        line = line.strip().split(',')
        individual = line[1]
        sample = line[4]
        R1 = args.fastq_dir+'/'+sample+'.R1.fastq.gz'
        R2 = args.fastq_dir+'/'+sample+'.R2.fastq.gz'
        if 'R2986' in R1:
            print(R1)
        if not os.path.exists(R1):
            continue            
        if not os.path.exists(R2):
            raise RuntimeError('BipSeq should all be paired end')

        out.write(sample+',BipSeq,'+individual+','+R1+','+R2+'\n')
