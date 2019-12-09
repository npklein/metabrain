
# loop through all subdirs to find all fastq files
import os.path
import re
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for MSBB.')
parser.add_argument('samplesheet', help='MSBB samplesheet from synapse')
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
        sample_id = line[1].replace('"','')
        individual_id = line[6].replace('"','')
        line[2] = line[2].replace('"','')
        if line[2].endswith('.bam'):
            continue
        if sample_id in samples:
            raise RuntimeError('Sample already seen')
        if sample_id in samples:
            samples.add(sample_id)

        fastq1 = args.fastq_dir+'/'+sample_id+'.accepted_hits.sort.coord.combined.fastq.gz'
        # MSBB is single end
        fastq2 = ''

        # Have to point to the bamfile
        if not os.path.exists(fastq1):
            print(fastq1+' does not exist')
#            raise RuntimeError(fastq1+' does not exist')

        if index % 25 == 0:
            batch_number += 1
            batch = 'batch'+str(batch_number)
            if out:
                out.close()
            out = open(args.samplesheet_dir+'samplesheet_MSBB_RNA.'+batch+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
    
        out.write(sample_id+',MSBB,'+individual_id+','+fastq1+','+fastq2+'\n')

        if sample_id.endswith('_resequenced'):
            with open(args.samplesheet_dir+'samplesheet_MSBB_RNA.mergedSamples.txt','a') as out_reseq:
                out_reseq.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
                fastq1 = fastq1.replace('_resequenced','').replace('.fastq.gz','.mergedWithResequenced.fastq.gz')
                if not os.path.exists(fastq1):
                    print(fastq1+' does not exist')
#                    raise RuntimeError(fastq1+' does not exist')
                out_reseq.write(sample_id.replace('_resequenced','_mergedWithResequenced')+',MSBB,'+individual_id+','+fastq1+','+fastq2+'\n')
        index += 1
out.close()
