import os.path
import os
import glob
import re
import argparse

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for Braineac.')
parser.add_argument('fastq_dir', help='directory containing fastq files')
parser.add_argument('original_fastq_files', help='file with list of downloaded fastq files')
parser.add_argument('sample_info_biogen',help='samplesheet provided by biogen')
parser.add_argument('outdir',help='Directory where output is written', default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

samples = set([])
all_fastq_files = set([])
with open(args.original_fastq_files) as input_file:
    for line in input_file:
        sample = '_'.join(line.split('/')[-1].split('_')[:2])
        # sample name is different between file and samplesheet, so adjust
        if sample == 'A653-1341':
            sample = 'A653_1341'
        samples.add(sample)
        all_fastq_files.add(args.fastq_dir+sample+'.R1.fastq.gz')

samples_per_biobank = {}
number_of_samples_per_biobank = {}
rin = {}
# SampleInfoBiogen is the samplesheet as provided to us, get information on individuals and fastq files from there
with open(args.sample_info_biogen) as input_file:
    header = input_file.readline().split(',')
    for line in input_file:
        line = line.strip().split(',')
        sample = line[0]
        biobank = line[3]
        individual_id = line[2]
        if biobank not in samples_per_biobank:
            samples_per_biobank[biobank] = []
            number_of_samples_per_biobank[biobank] = 0
        # <<<>>> is just something that can split on later
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
        # make batches of size 25
        if number_of_samples_per_biobank[biobank] % 25 == 0:
            if out:
                out.close()
            batch = biobank+'_batch'+str(int(number_of_samples_per_biobank[biobank]/25))
            out = open(args.outdir+'samplesheet_Braineac_RNA.'+batch+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
           
        out.write(sample_name+',Braineac,'+individual_id+','+R1+','+R2+'\n')
        all_fastq_files.remove(R1)
        number_of_samples_per_biobank[biobank] += 1

if len(all_fastq_files) > 0:
    raise RuntimeError('Not all fastq files included')
out.close()
