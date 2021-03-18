
# loop through all subdirs to find all fastq files
import os.path
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for BipSeq.')
parser.add_argument('samplesheet', help='BipSeq samplesheet from synapse')
parser.add_argument('cram_dir', help='path to cram file dir for input')
parser.add_argument('fastq_dir', help='path to fastq file dir for output')
parser.add_argument('samplesheet_dir',help='Directory where samplesheet is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

individual_per_sample = {}
samples_per_batch = {}
batch_count = {}

batch = 0
out = None
# BipSeq is small enough to not need batches
with open(args.samplesheet) as input_file:
    header = input_file.readline().split(',')
    for index, line in enumerate(input_file):
        if index % 25 == 0:
            batch += 1 
            if out:
                out.close()
            batchname = 'batch'+str(batch)
            out = open(args.samplesheet_dir+'/samplesheet_BipSeq_RNA.'+batchname+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')

        line = line.strip().split(',')
        individual = line[1]
        sample = line[4]

        R1 = args.fastq_dir+'/'+sample+'.R1.fastq.gz'
        R2 = args.fastq_dir+'/'+sample+'.R2.fastq.gz'
        cram = args.cram_dir+'/'+individual+'_'+sample+'.cram'
        if not os.path.exists(cram):
            print(cram,'does not exist')

        out.write(sample+',BipSeq,'+individual+','+R1+','+R2+','+cram+'\n')
out.close()
