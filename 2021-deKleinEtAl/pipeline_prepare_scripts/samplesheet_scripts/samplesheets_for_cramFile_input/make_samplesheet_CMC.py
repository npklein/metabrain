
# loop through all subdirs to find all fastq files
import os.path
import re
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for CMC.')
parser.add_argument('samplesheet', help='CMC samplesheet from synapse')
parser.add_argument('cram_dir', help='path to cram file dir for input')
parser.add_argument('fastq_dir', help='path to fastq file dir for output')
parser.add_argument('samplesheet_dir',help='Directory where samplesheet is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

individual_per_sample = {}
samples_per_batch = {}
study_batch_order = []
with open(args.samplesheet) as input_file:
    header = input_file.readline().split('\t')
    batch = {}
    index = {}
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

        fastq_path = args.fastq_dir+'/'+sample_id+'_R1.fastq.gz'

        # Have to point to the cramfile
        cram = args.cram_dir+'/individualID.'+individual_id+'_specimenID.'+sample_id+'.cram'
        if not os.path.exists(cram):
            cram = args.cram_dir.replace('no_patch','patch')+'/individualID.'+individual_id+'_specimenID.'+sample_id+'.cram'
            if not os.path.exists(cram):
                raise RuntimeError(cram+' does not exist')

        study_subset = re.search('(.*?)_\d+', sample_id).group(1)
        if study_subset not in index:
            index[study_subset] = 0
            batch[study_subset] = 0

        # Make batch of 25 samples
        if index[study_subset] % 25 == 0:
            batch[study_subset] += 1
        index[study_subset] += 1
        batch_name = study_subset+'_batch'+str(batch[study_subset])
        if batch_name not in samples_per_batch:
            samples_per_batch[batch_name] = []
            study_batch_order.append(batch_name)
        samples_per_batch[batch_name].append([fastq_path, sample_id,cram])

for batch in study_batch_order:
    print('batch: '+batch)
    with open(args.samplesheet_dir+'samplesheet_CMC_RNA.'+batch+'.txt','w') as out:
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n')
        for data in samples_per_batch[batch]:
#            print(data[1], individual_per_sample[data[1]])
            out.write(data[1]+',CMC,'+individual_per_sample[data[1]]+','+data[0]+','+data[0].replace('R1','R2')+','+data[2]+'\n')
