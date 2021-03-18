# make merge_fastq.sh which merges FastQ files together for Brainseq
# FastQ files and phenotype file are downloaded from Synapse

import glob

samples = []
# Get the sample names from the phenotype file
with open('phenotype_data/phenotypeFile_LIBD_szControl.csv') as input_file:
    header = input_file.readline().split(',')
    for line in input_file:
        line = line.strip().split(',')
        rnum = line[0]
        samples.append(rnum.strip('"'))
    
with open('merge_fastq.sh','w') as out:
    for sample in samples:
        fastq_files_R1 = glob.glob('RNAseq/'+sample+'_*_R1_*gz')
        out.write('cat '+' '.join(fastq_files_R1)+' > RNAseq_merged/'+sample+'.R1.fastq.gz\n')
        
        # this way the fastq files are in same order, which is necesarry because the paired ends need to be in same order
        fastq_files_R2 = [x.replace('R1','R2') for x in fastq_files_R1]
        out.write('cat '+' '.join(fastq_files_R2)+' > RNAseq_merged/'+sample+'.R2.fastq.gz\n')
        for f in fastq_files_R1:
            out.write('rm '+f+'\n')
        for f in fastq_files_R2:
            out.write('rm '+f+'\n')
        out.write('\n')
