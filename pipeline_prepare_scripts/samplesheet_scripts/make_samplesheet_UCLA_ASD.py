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

batch_size = 25
batch_number = 0
out = None
fastq_files_initial = set(glob.glob(args.fastq_dir+'/*'))
fastq_files = set([])
for f in fastq_files_initial:
    while '//' in f:
        fastq_files.replace('//','/')
    fastq_files.add(f)

with open(args.samplesheet) as input_file:
    header = input_file.readline().split(',')
    for index, line in enumerate(input_file):
        if 'UMB5308' in line:
            print(line)
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
        sample = line[2].replace('Sample_','').replace('_1st','').replace('_2nd','').replace('_3rd','')
        R1 = args.fastq_dir+'/'+sample+'.R1.fastq.gz'
        R2 = args.fastq_dir+'/'+sample+'.R2.fastq.gz'
        if not os.path.isfile(R1):
            raise RuntimeError(R1+' does not exist')
        if not os.path.isfile(R2):
            raise RuntimeError(R2+' does not exist')
        out.write(sample+',UCLA_ASD,'+individual+','+R1+','+R2+'\n')

        while '//' in R1:
            R1 = R1.replace('//','/')
        try:
            fastq_files.remove(R1)
        except KeyError:
            if '1st' in line[2] or '2nd' in line[2] or '3rd' in line[2]:
                pass
            else:
                raise
        while '//' in R2:
            R2 = R2.replace('//','/')
        try:
            fastq_files.remove(R2)
        except KeyError:
            if '1st' in line[2] or '2nd' in line[2] or '3rd' in line[2]:
                pass
            else:
                raise

out.close()
print('Missing samples:')
print(fastq_files)
