import sys
import os.path
import glob
import os
import glob
import re

# batch size is set to 10 because aspera does not allow more than 10 connections at the same time 
# (or can not get more than 10 through the port)
batch_size = 10

# argument 1: location of the samplesheet
# argument 2: location where to put the fastq files
if len(sys.argv) != 3:
    raise RuntimeError('Not correct number of arguments')


def get_all_indices(list_to_index):
    '''Get the indexes of all the items in a list and put them in a dict with key: element, value: index
        
       list_to_index(list)    List to get index from all elements from
    '''
    list_indexes = {}
    i = 0
    for element in list_to_index:
        list_indexes[element] = i
        i += 1
    return list_indexes

fastq_dir = sys.argv[2]


seen = set([])
study_individuals = {}
samples_per_individual = {}
fq_files_per_sample = {}
with open(sys.argv[1]) as input_file:
    header = input_file.readline().strip().split('\t')
    # Because the samplesheet has many headers that might change, get the index of each of hte headers
    header_by_index = get_all_indices(header)
    for line in input_file:
        line = line.strip().split('\t')
        # individual can have multiple samples, need to know for merging bams later
        sample = line[header_by_index['run_accession']]
        individual = line[header_by_index['sample_accession']]
        # need the links to see if it is single or paired end sequencing
        fastq_aspera_links = line[header_by_index['fastq_aspera']].rstrip(';').split(';')
        study = line[0]
        if study not in study_individuals:
            study_individuals[study] = set([])
        study_individuals[study].add(individual)
        fq_files_per_sample[sample] = fastq_aspera_links
        # need to know which samples belong to which individual so they can be added to the sameb batch
        if individual not in samples_per_individual:
            samples_per_individual[individual] = set([])
        samples_per_individual[individual].add(sample)

        # samples should be unique
        if sample in seen:
            raise RuntimeError('Sample '+sample+' already seen')
        seen.add(sample)

seen = set([])
# Divide the jobs up per study, easier to keep track later (although some studies have only 1 sample)
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
        if total_number_of_samples % batch_size == 0 or current_number_of_samples + len(samples_per_individual[individual]) > 10:
            if out:
                out.close()
            
            batch_number = int(total_number_of_samples / batch_size)
            batch = 'batch'+str(batch_number)
            out = open('Public_RNA-seq_QC/samplesheets/samplesheet_ENA_RNA.'+study+'_'+batch+'.txt','w')
            out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,sortedBamFile\n')
            current_number_of_samples = 0

        # Because some of the samples from an indivudal might not have been added to this batch, do it now before closing the file
        # samples_per_individual <- all samples of individual
        # samples_included_per_individual <- samples includeded per individual in this batch
        # The difference is the samples that have not yet been included
        for sample in sorted(samples_per_individual[individual]):
            if sample in seen:
                raise RuntimeError(sample+' seen for 2nd time')
            seen.add(sample)
            fastq_files = fq_files_per_sample[sample]
            R1 = fastq_dir+'/'+study+'/'+fastq_files[0].split('/')[-1]
            R2 = ''
            if len(fastq_files) > 1:
                R2 = fastq_dir+'/'+study+'/'+fastq_files[1].split('/')[-1]
            out.write(individual+',ENA,'+sample+','+R1+','+R2+',${sortedBam}\n')
        total_number_of_samples += len(samples_per_individual[individual])
        current_number_of_samples += len(samples_per_individual[individual])

    out.close()
