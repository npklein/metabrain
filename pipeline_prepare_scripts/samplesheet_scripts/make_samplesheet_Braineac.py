# loop through all subdirs to find all fastq files
import sys
import os.path
import glob
import os
import glob
import re

fastq_dir = '/groups/umcg-biogen/tmp03/input/ucl-upload-biogen/data/fastq/'
samples = set([])
all_fastq_files = set([])
with open('data/original_fastq_files.txt') as input_file:
    for line in input_file:
        sample = '_'.join(line.split('/')[-1].split('_')[:2])
        if sample == 'A653-1341':
            sample = 'A653_1341'
        samples.add(sample)
        all_fastq_files.add(fastq_dir+sample+'.R1.fastq.gz')

samples_per_biobank = {}
number_of_samples_per_biobank = {}
rin = {}
with open('SampleInfoBiogen.csv') as input_file:
    header = input_file.readline().split(',')
    for line in input_file:
        line = line.strip().split(',')
        sample = line[0]
        biobank = line[3]
        individual_id = line[2]
        if biobank not in samples_per_biobank:
            samples_per_biobank[biobank] = []
            number_of_samples_per_biobank[biobank] = 0
        samples_per_biobank[biobank].append(sample+'<<<>>>'+individual_id)         
        rin[sample] = line[-1]

out = None
for biobank in samples_per_biobank:
    for sample in samples_per_biobank[biobank]:
        sample_name = sample.split('<<<>>>')[0]
        individual_id = sample.split('<<<>>>')[1].replace('/','_')
        if sample_name == 'A653-1341':
            sample_name = 'A653_1341'
        R1 = fastq_dir+sample_name+'.R1.fastq.gz'
        R2 = fastq_dir+sample_name+'.R2.fastq.gz'
        if sample_name not in samples:
            continue
        if number_of_samples_per_biobank[biobank] % 25 == 0:
            if out:
                out.close()
            batch = biobank+'_batch'+str(int(number_of_samples_per_biobank[biobank]/25))
            out = open('Public_RNA-seq_QC/samplesheets/samplesheet_Braineac_RNA.'+batch+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
           
        out.write(sample_name+',Braineac,'+individual_id+','+R1+','+R2+'\n')
        all_fastq_files.remove(R1)
        number_of_samples_per_biobank[biobank] += 1

if len(all_fastq_files) > 0:
    raise RuntimeError('Not all fastq files included')
out.close()
