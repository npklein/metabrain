# Parse TargetALS samplesheet to make molgenis-compute samplesheet
import os.path
import argparse

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for TargetALS.')
parser.add_argument('samplesheet', help='TargetALS samplesheet')
parser.add_argument('cram_dir', help='path to cram file dir for input')
parser.add_argument('fastq_dir', help='path to fastq file dir for output')
parser.add_argument('samplesheet_dir',help='Directory where samplesheet is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

cram_seen = set([])
info_per_quote = {}
with open(args.samplesheet) as input_file:
    header = input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        quote = line[0]
        external_sample_id = line[1]
        external_subject_id = line[2]


        if quote not in info_per_quote:
            info_per_quote[quote] = []
        external_subject_id = external_subject_id.replace(' ','_')
        external_sample_id = external_sample_id.replace('-HRA-','_HRA_')

#/groups/umcg-biogen/tmp04/input/rawdata/TargetALS/pipelines/patch_chromosomes/results/cramFiles/NEUKX037NKA_CGND_HRA_01218.cram
#/groups/umcg-biogen/tmp04/input/rawdata/TargetALS/pipelines/patch_chromosomes/results/cramFiles/NEUKX037NKA_CGND-HRA-01218.cram

        fastq1 = args.fastq_dir+'/'+external_subject_id+'_'+external_sample_id+'.R1.fastq.gz'
        fastq2 = args.fastq_dir+'/'+external_subject_id+'_'+external_sample_id+'.R2.fastq.gz'
        cram1 = args.cram_dir+'/'+external_subject_id+'_'+external_sample_id+'.cram'

        for i in range(0,10):
            cram1 = cram1.replace('JHU'+str(i),'JHU_'+str(i))
        # below is probably the worst code ever written. If you have to debug it, I'm sorry.
        # cram file names were sometimes wrong due to a mistake in the samplesheet. 
        # the mistake in the samplesheet has been fixed, but now the cram files
        # are not matching with the samplesheet. Have to check which is the correct one
        if not os.path.exists(cram1):
            cram3 = cram1.replace('_HRA_','-HRA-')
            if not os.path.exists(cram3):
                try:
                    cram2 = args.cram_dir+'/'+line[3]+'_'+external_sample_id+'.cram'
                except:
                    print(cram1)
                    print(cram3)
                    raise RuntimeError(cram1+' does not exist')
                if not os.path.exists(cram2):
                    cram4 = cram2.replace('_HRA_','-HRA-')
                    if not os.path.exists(cram4):
                        print(cram1)
                        print(cram2)
                        print(cram3)
                        print(cram4)
                        print('all not found')
                        print('-'*20)
#                        raise RuntimeError(cram2+' does not exist')
                    cram = cram4
                cram = cram2
            else:
                cram = cram3
        else:
            cram = cram1
#            raise RuntimeError(cram+' does not exist')
        if cram in cram_seen:
#            print(cram+' already added, continue')
            continue

        cram_seen.add(cram)
        info_per_quote[quote].append(external_sample_id+',TargetALS,'+external_subject_id+','+fastq1+','+fastq2+','+cram+'\n')

for quote in info_per_quote:
    with open(args.samplesheet_dir+'/samplesheet_TargetALS_RNA.'+quote+'.txt', 'w') as out:
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')
        for info in info_per_quote[quote]:
            out.write(info)
