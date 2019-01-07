import sys
import os.path
import glob
import os
import glob
import re

# TODO: remove hard coded path
fastq_dir = '/groups/umcg-biogen/tmp03/input/rawdata/Brainseq/RNAseq_merged/'
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
    R1 = fastq_dir+sample_name+'.R1.fastq.gz'
    R2 = fastq_dir+sample_name+'.R2.fastq.gz'
    # make batches of size 25
    if number_of_samples % 25 == 0:
        if out:
            out.close()
        batch = 'batch'+str(int(number_of_samples/25))
        out = open('Public_RNA-seq_QC/samplesheets/samplesheet_Brainseq_RNA.'+batch+'.txt','w')
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
           
    out.write(sample_name+',Braineac,'+individual_id+','+R1+','+R2+'\n')
    number_of_samples += 1

out.close()
