# loop through all subdirs to find all fastq files
import os.path
import glob
import os
import re
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for BrianGVEx.')
parser.add_argument('samplesheet', help='BrianGVEx samplesheet from synapse')
parser.add_argument('fastq_dir', help='path to fastq file dir')
parser.add_argument('--outdir',help='Directory where output is written',default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

individual_per_sample = {}
samples_per_batch = {}
batch_count = {}

# BrianGVEx is small enough to not need batches
batch_size = 25
batch_number = 0
out = None
with open(args.samplesheet) as input_file:
    header = input_file.readline().split(',')
    for index, line in enumerate(input_file):
        if len(line.split(',')[1].strip()) == 0:
            continue
        if index % 25 == 0:
            if out:
                out.close()
            out = open(args.outdir+'samplesheet_UCLA_ASD_RNA.batch'+str(batch_number)+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
            batch_number += 1
        line = line.strip().split(',')
        individual = line[1]
        sample = line[1]
        R1 = args.fastq_dir+'/'+sample+'.R1.fastq.gz'
        R2 = args.fastq_dir+'/'+sample+'.R2.fastq.gz'

        out.write(sample+',UCLA_ASD,'+individual+','+R1+','+R2+'\n')
out.close()
