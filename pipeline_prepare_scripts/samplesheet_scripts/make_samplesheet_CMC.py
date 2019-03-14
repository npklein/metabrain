
# loop through all subdirs to find all fastq files
import os.path
import glob
import os
import re
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for CMC.')
parser.add_argument('samplesheet', help='CMC samplesheet from synapse')
parser.add_argument('fastq_dir', help='path to fastq file dir')
parser.add_argument('outdir',help='Directory where output is written',default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

individual_per_sample = {}
samples_per_batch = {}
batch_count = {}

with open(args.samplesheet) as input_file:
    header = input_file.readline().split('\t')
    for line in input_file:
        line = line.strip().replace('"','').split('\t')
        aligned_bam = line[4]
        if not aligned_bam.endswith('.accepted_hits.sort.coord.bam'):
            continue
        sample_id = line[-1]
        individual_id = line[-2]
        if sample_id in individual_per_sample:
            raise RuntimeError('Sample ID in multiple times')
        individual_per_sample[sample_id] = individual_id

        # The reads are downloaded as BAM files, one with aligned reads, one with the unaligned reads. These will be
        # merged together later, so need to have both of the file names
        unaligned_bam = aligned_bam.replace('BamAlignedReadData','BamUnmappedReadData').replace('accepted_hits.sort.coord.bam','unmapped.bam')

        fastq_path = args.fastq_dir+'/'+sample_id+'_R1.fq.gz'
        batch = re.search('(.*?)_\d+', sample_id).group(1)
        if batch not in batch_count:
            batch_count[batch] = 0
        batch_count[batch] += 1
        # Make max 5 batches, where the first 4 have 25 samples and the 5th has al other sampls
        if batch_count[batch] > 25 and batch_count[batch] < 50:
            batch = batch+'_batch2'
        elif batch_count[batch] >= 50 and batch_count[batch] < 75:
            batch = batch+'_batch3'
        elif batch_count[batch] >= 75 and batch_count[batch] < 100:
            batch = batch+'_batch4'
        elif batch_count[batch] >= 100:
            batch = batch+'_batch5'
        print(batch)
        if batch not in samples_per_batch:
            samples_per_batch[batch] = []
        samples_per_batch[batch].append([fastq_path, sample_id, aligned_bam, unaligned_bam])

for batch in samples_per_batch:
    print(batch)
    with open(args.outdir+'samplesheet_CMC_RNA.'+batch+'.txt','w') as out:
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBam,unalignedBam\n')
        for data in samples_per_batch[batch]:
            print(data[1], individual_per_sample[data[1]])
            out.write('specimenID.'+data[1]+',CMC,'+'individualID.'+individual_per_sample[data[1]]+','+data[0]+','+data[0].replace('R1','R2')+',')
            out.write('/groups/umcg-biogen/tmp03/input/CMC/RNAseq/DorsolateralPrefrontalCortex/Raw/BamAlignedReadData/'+data[2])
            out.write(',/groups/umcg-biogen/tmp03/input/CMC/RNAseq/DorsolateralPrefrontalCortex/Raw/BamUnmappedReadData/'+data[3]+'\n')
