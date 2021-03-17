
# loop through all subdirs to find all fastq files
import os.path
import re
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for MayoTCX.')
parser.add_argument('samplesheet', help='MayoTCX samplesheet from synapse')
parser.add_argument('fastq_dir', help='path to fastq file dir for input files')
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
        if line[1] == 'NA':
            continue
        sample_id = line[0]
        individual_id = sample_id
        if sample_id in samples:
            raise RuntimeError('Sample already seen')
        if sample_id in samples:
            samples.add(sample_id)
        fastq1 = args.fastq_dir+'/'+sample_id+'.r1.fastq.gz'
        fastq2 = fastq1.replace('.r1.','.r2.')

        # Have to point to the bamfile
        if not os.path.exists(fastq1):
            raise RuntimeError(fastq1+' does not exist')
        if not os.path.exists(fastq2):
            raise RuntimeError(fastq2+' does not exist')

        if index % 25 == 0:
            batch_number += 1
            batch = 'batch'+str(batch_number)
            if out:
                out.close()
            out = open(args.samplesheet_dir+'samplesheet_MayoTCX_RNA.'+batch+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
    
        out.write(sample_id+',MayoTCX,'+individual_id+','+fastq1+','+fastq2+'\n')

        index += 1
out.close()
