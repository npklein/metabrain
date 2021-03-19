import sys
import os.path
import glob
import os
import glob
import re

batch_size = 25

samples_per_individual = {}
study_individuals = {}
assay_and_library_type = {}
with open(sys.argv[1]) as input_file:
    header = input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        assay_type = line[0]
        individual = line[3]
        sample = line[14]
        library = line[7]
        bio_project = line[2]
        
        if assay_type == 'WXS':
            continue
        
        # samples per individual
        if individual not in samples_per_individual:
            samples_per_individual[individual] = set([])
        samples_per_individual[individual].add(sample)

        # assay and library type per sample
        assay_and_library_type[sample] = [assay_type, library]


        # individual per bio_project(study)
        if bio_project not in study_individuals:
            study_individuals[bio_project] = set([])
        study_individuals[bio_project].add(individual)


fastq_dir = sys.argv[2]
seen = set([])
# Divide the jobs up per bio project, easier to keep track later 
# Since some studies are very large, make batches per <batch_size> samples, but make sure that samples from same individual are in the same batch
for study in study_individuals:
    total_number_of_samples = 0
    current_number_of_samples = 0
    out = None
    individuals_in_current_batch = set([])
    batch_number = 0
    batch = ''
    prev_individual = None
    individuals = sorted(study_individuals[study])
    for individual in individuals:
        # check if number_of_samples is equal to batch size
        # or that if the total number of samples for next batch plus current number in batch > batch_size (so no batch gets larger than batch_size)
        # This to keep number of jobs low enough for aspera
        if total_number_of_samples % batch_size == 0 or current_number_of_samples + len(samples_per_individual[individual]) > batch_size:
            if out:
                out.close()

            batch_number = int(total_number_of_samples / batch_size)
            batch = 'batch'+str(batch_number)
            out = open('Public_RNA-seq_QC/samplesheets/samplesheet_NABEC_RNA.'+study+'_'+batch+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
            current_number_of_samples = 0

        # Because some of the samples from an indivudal might not have been added to this batch, do it now before closing the file
        # samples_per_individual <- all samples of individual
        # samples_included_per_individual <- samples includeded per individual in this batch
        # The difference is the samples that have not yet been included
        for sample in sorted(samples_per_individual[individual]):
            if sample in seen:
                raise RuntimeError(sample+' seen for 2nd time')
            seen.add(sample)
            R1 = ''
            R2 = ''
            if assay_and_library_type[sample][1] == 'PAIRED':
                R1 = fastq_dir+'/'+sample+'_1.fastq.gz'
                R2 = fastq_dir+'/'+sample+'_2.fastq.gz'
            elif assay_and_library_type[sample][1] == 'SINGLE':
                R1 = fastq_dir+'/'+sample+'.fastq.gz'

            out.write(sample+',NABEC,'+individual+','+R1+','+R2+'\n')
        total_number_of_samples += len(samples_per_individual[individual])
        current_number_of_samples += len(samples_per_individual[individual])

    out.close()

           

