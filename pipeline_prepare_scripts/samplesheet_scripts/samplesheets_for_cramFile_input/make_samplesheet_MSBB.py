
# loop through all subdirs to find all fastq files
import os.path
import re
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for MSBB.')
parser.add_argument('samplesheet', help='MSBB samplesheet from synapse')
parser.add_argument('bam_dir', help='path to bam file dir for input')
parser.add_argument('fastq_dir', help='path to fastq file dir for output')
parser.add_argument('samplesheet_dir',help='Directory where samplesheet is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

batch_number = 0
batch = None
out = None
samples = set([])
with open(args.samplesheet) as input_file:
    header = input_file.readline().split(',')
    for index, line in enumerate(input_file):
        line = line.strip().split(',')
        sample_id = line[1].replace('"','')
        individual_id = line[6].replace('"','')
        line[2] = line[2].replace('"','')
        if not line[2].endswith('.bam'):
            continue
        if sample_id in samples:
            raise RuntimeError('Sample ID '+sample_id+' in multiple times')
        samples.add(sample_id)

        fastq1 = args.fastq_dir+'/'+sample_id+'.fq.gz'
        # MSBB is single end
        fastq2 = ''

        # Have to point to the bamfile
        bam = args.bam_dir+'/'+line[2].replace('"','').replace('sort.coord.bam','sort.coordAligned.out.bam')
        if not os.path.exists(bam):
            print(sample_id+' '+bam+' does not exist')
            continue
#            raise RuntimeError(bam+' does not exist')

        if index % 25 == 0:
            batch_number += 1
            batch = 'batch'+str(batch_number)
            if out:
                out.close()
            out = open(args.samplesheet_dir+'samplesheet_MSBB_RNA.'+batch+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')
    
        out.write(sample_id+',MSBB,'+individual_id+','+fastq1+','+fastq2+','+bam+'\n')

        if sample_id.endswith('_resequenced'):
            with open(args.samplesheet_dir+'samplesheet_MSBB_RNA.mergedSamples.txt','a') as out_reseq:
                out_reseq.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')
                bam = bam.replace('_resequenced','').replace('.bam','.mergedWithResequenced.bam')
                if not os.path.exists(bam):
                    pass
#                    raise RuntimeError(sample_id+' '+bam+' does not exist')
                out_reseq.write(sample_id.replace('_resequenced','_mergedWithResequenced')+',MSBB,'+individual_id+','+fastq1+','+fastq2+','+bam+'\n')

out.close()
