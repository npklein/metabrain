# Parse GTEx cram files to make molgenis-compute samplesheet

import os.path
import argparse

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for GTEx.')
parser.add_argument('samplesheet', help='Enter delimited file with List of GTEx cramfiles')
parser.add_argument('fastq_dir', help='path to fastq file dir for output')
parser.add_argument('samplesheet_dir',help='Directory where samplesheet is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()


# This simply uses all the cram files, not a samplesheet. 
# Better to rewrite to use a samplesheet. However, such samplesheet
# was not available at time of writing this script

batch_size = 25
batch = 0
out = None
with open(args.samplesheet) as input_file:
    for index,line in enumerate(input_file):
        if index % batch_size == 0:
            batch += 1
            batchname = 'batch'+str(batch)
            print(batchname)
            if out:
                out.close()
            out = open('Public_RNA-seq_QC/samplesheets/samplesheet_GTEx_RNA.'+batchname+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')

        sampleName = line.strip().split('/')[-1].split('.cram')[0].split('_')[0]
        cram = line.strip()
        if not os.path.exists(cram):
            cram = cram.replace('/patch','/no_patch')
            if not os.path.exists(cram):
                raise RuntimeError(cram+' does not exist')
        fastq1 = args.fastq_dir+'/'+sampleName+'_1.fastq.gz'
        fastq2 = args.fastq_dir+'/'+sampleName+'_2.fastq.gz'
        out.write(sampleName+',GTEx,'+sampleName+','+fastq1+','+fastq2+','+cram+'\n')
