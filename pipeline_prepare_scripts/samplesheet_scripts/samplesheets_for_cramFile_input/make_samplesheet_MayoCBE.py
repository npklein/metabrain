
# loop through all subdirs to find all fastq files
import os.path
import re
import argparse 
import glob

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for MayoCBE.')
parser.add_argument('samplesheet', help='MayoCBE samplesheet from synapse')
parser.add_argument('cram_dir', help='path to cram file dir for input')
parser.add_argument('fastq_dir', help='path to fastq file dir for output')
parser.add_argument('samplesheet_dir',help='Directory where samplesheet is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')


args = parser.parse_args()

batch_number = 0
batch = None
out = None
samples = set([])
imdex = 0
with open(args.samplesheet) as input_file:
    header = input_file.readline().split(',')
    index = 0
    for line in input_file:
        line = line.strip().split(',')

        sample_id = line[0]
        individual_id = sample_id
        if sample_id in samples:
            raise RuntimeError('Sample already seen')
        if sample_id in samples:
            samples.add(sample_id)
        fastq1 = args.fastq_dir+'/'+individual_id+'_'+sample_id+'.r1.fastq.gz'
        fastq2 = fastq1.replace('.r1.','.r2.')

        cram = glob.glob(args.cram_dir+'/'+sample_id+'.*bam')
        if len(cram) == 0:
            continue
        bam = cram[0]
        print(cram)

        if index % 25 == 0:
            batch_number += 1
            batch = 'batch'+str(batch_number)
            if out:
                out.close()
            out = open(args.samplesheet_dir+'samplesheet_MayoCBE_RNA.'+batch+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')
    
        out.write(sample_id+',MayoCBE,'+individual_id+','+fastq1+','+fastq2+','+bam+'\n')

        index += 1
out.close()
