import os.path
import os
import glob
import re
import argparse

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for Brainseq.')
parser.add_argument('fastq_dir', help='directory containing fastq files')
parser.add_argument('phenotype_file', help='file with sample info')
parser.add_argument('outdir',help='Directory where output is written', default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

samples = []
with open('/groups/umcg-biogen/tmp03/input/rawdata/Brainseq/phenotype_data/phenotypeFile_LIBD_szControl.csv') as input_file:
    header = input_file.readline().split(',')
    for line in input_file:
        line = line.strip().split(',')
        rnum = line[0].strip('"')
        bnum = line[1].strip('"')
        samples.append(rnum+'<<>>'+bnum)


out = None
number_of_samples = 0
for sample in samples:
    sample_name = sample.split('<<>>')[0]
    individual_id = sample.split('<<>>')[1].replace('/','_')
    R1 = args.fastq_dir+sample_name+'.R1.fastq.gz'
    R2 = args.fastq_dir+sample_name+'.R2.fastq.gz'
    # make batches of size 25
    if number_of_samples % 25 == 0:
        if out:
            out.close()
        batch = 'batch'+str(int(number_of_samples/25))
        out = open(args.outdir+'/samplesheet_Brainseq_RNA.'+batch+'.txt','w')
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
           
    out.write(sample_name+',Braineac,'+individual_id+','+R1+','+R2+'\n')
    number_of_samples += 1

out.close()
