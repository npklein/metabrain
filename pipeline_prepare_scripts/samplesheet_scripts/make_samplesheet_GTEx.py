# Parse GTEx fastq files to make molgenis-compute samplesheet
import sys
import os.path
import glob

# This simply uses all the fq files, not a smaplesheet. Better to rewrite to use a samplesheet
print('Not implemeted yet')
exit()

# sys.argv[1] should be directory with the batches that contain fastq files
if len(sys.argv) != 3:
    raise RuntimeError('Not correct number of arguments')

for d in glob.glob(sys.argv[1]+'/'+sys.argv[2]):
    print(d)
    batch = d.split('/')[-1]
    print(batch)
    with open('Public_RNA-seq_QC/samplesheets/samplesheet_GTEx_RNA.'+batch+'.txt', 'w') as out:
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
        for f in glob.glob(d+'/*_1.fastq.gz'):
            sampleName = f.split('/')[-1].split('_1.fastq.gz')[0]
            fastq1 = f
            fastq2 = f.replace('_1.fastq.gz','_2.fastq.gz')
            out.write(sampleName+',GTEx,'+sampleName+','+fastq1+','+fastq2+'\n')
