
# loop through all subdirs to find all fastq files
import os.path
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for CMC_HBCC.')
parser.add_argument('samplesheet', help='CMC samplesheet from synapse')
parser.add_argument('cram_dir', help='path to cram file dir for input')
parser.add_argument('fastq_dir', help='path to fastq file dir for output')
parser.add_argument('samplesheet_dir',help='Directory where samplesheet is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')
args = parser.parse_args()

individual_per_sample = {}
samples_per_batch = {}
batch_count = 0
batch_number = 0
batch = 'batch0'
with open(args.samplesheet) as input_file:
    header = input_file.readline().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        sample_id = line[0]
        individual_id = line[1]

        if sample_id in individual_per_sample:
            raise RuntimeError('Sample ID in multiple times')
        individual_per_sample[sample_id] = individual_id

        fastq_path = args.fastq_dir+'/'+sample_id+'_R1.fq.gz'

        batch_count += 1
        if batch_count % 25 == 0:
            batch_number += 1
            batch = 'batch'+str(batch_number)

        if batch not in samples_per_batch:
            samples_per_batch[batch] = []

        cram = args.cram_dir+'/individualID.'+individual_id+'_specimenID.'+sample_id+'.cram'
        if not os.path.exists(cram):
#            print(cram+' does not exist')
            raise RuntimeError(cram+' does not exist')

        samples_per_batch[batch].append([fastq_path, sample_id, cram])

for batch in samples_per_batch:
    with open(args.samplesheet_dir+'/samplesheet_CMC_HBCC_RNA.'+batch+'.txt','w') as out:
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')
        for data in samples_per_batch[batch]:
            out.write('specimenID.'+data[1]+',CMC_HBCC,'+'individualID.'+individual_per_sample[data[1]]+','+data[0]+','+data[0].replace('R1','R2')+',')
            out.write(data[2]+'\n')
