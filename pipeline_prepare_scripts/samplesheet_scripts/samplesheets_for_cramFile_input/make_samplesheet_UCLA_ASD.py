# loop through all subdirs to find all fastq files
import os.path
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for BrianGVEx.')
parser.add_argument('samplesheet', help='UCLA_ASD samplesheet from synapse')
parser.add_argument('cram_dir', help='path to cram file dir for input')
parser.add_argument('fastq_dir', help='path to fastq file dir for output')
parser.add_argument('samplesheet_dir',help='Directory where samplesheet is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')


args = parser.parse_args()

individual_per_sample = {}
samples_per_batch = {}
batch_count = {}

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
            out = open(args.samplesheet_dir+'/samplesheet_UCLA_ASD_RNA.batch'+str(batch_number)+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')
            batch_number += 1
        line = line.strip().split(',')
        individual = line[1]
        sample = line[2].replace('Sample_','').replace('_1st','').replace('_2nd','').replace('_3rd','')
        R1 = args.fastq_dir+'/'+sample+'.R1.fastq.gz'
        R2 = args.fastq_dir+'/'+sample+'.R2.fastq.gz'
        cram = args.cram_dir+'/'+individual+'_'+sample+'.cram'

        if not os.path.isfile(cram):
            raise RuntimeError(cram+' does not exist')
        out.write(sample+',UCLA_ASD,'+individual+','+R1+','+R2+','+cram+'\n')

out.close()
