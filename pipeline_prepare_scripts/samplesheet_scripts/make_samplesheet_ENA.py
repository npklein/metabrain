import sys
import os.path
import glob
import os
import glob
import re

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
samples_to_add = set([])
fq_files = {}
with open(sys.argv[1]) as input_file:
    header = input_file.readline().strip().split('\t')
    header_by_index = get_all_indices(header)
    for line in input_file:
        line = line.strip().split('\t')
        sample = line[header_by_index['run_accession']]
        individual = line[header_by_index['sample_accession']]
        # need the links to see if it is single or paired end sequencing
        fastq_aspera_links = line[header_by_index['fastq_aspera']].rstrip(';').split(';')
        # <<<>>> is just a separator
        samples_to_add.add(sample+'<<<>>>'+individual)
        fq_files[sample+'<<<>>>'+individual] = fastq_aspera_links

        if sample+'<<<>>>'+individual in seen:
            raise RuntimeError('Combination of sample+individual already seen')
        seen.add(sample+'<<<>>>'+individual)

out = None
number_of_samples = 0
for sample in samples_to_add:
    sample_name = sample.split('<<<>>>')[0]
    individual_id = sample.split('<<<>>>')[1]

    

    # make batches of size 25
    if number_of_samples % 25 == 0:
        if out:
            out.close()
            exit()
        batch = 'batch'+str(int(number_of_samples/25))
        R1 = fastq_dir+'/'+batch+'/'+fq_files[sample][0].split('/')[-1]
        R2 = ''
        if len(fq_files[sample]) > 1:
            R2 = fastq_dir+'/'+batch+'/'+fq_files[sample][1].split('/')[-1]
        out = open('Public_RNA-seq_QC/samplesheets/samplesheet_ENA_RNA.'+batch+'.txt','w')
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz\n')
           
    out.write(individual_id+',ENA,'+sample_name+','+R1+','+R2+'\n')
    number_of_samples += 1

out.close()
