# Read through all .txt files in ../brain_samples/ and filter ena.txt to only include those samples

import glob
import datetime

samples_to_keep = set([])

for f in glob.glob('../brain_samples/*txt'):
    with open(f) as input_file:
        header = input_file.readline().strip()
        if header != 'brainSamplesToInclude':
            raise RuntimeError('Different header than expected. Change or add header "brainSamplesToInclude"')
        for line in input_file:
            line = line.strip()
            samples_to_keep.add(line)

today = datetime.datetime.today().strftime('%Y%m%d')
with open('ena.txt') as input_file, open('samplesheet_ENA_'+today+'.txt','w') as out:
    out.write(input_file.readline())
    for line in input_file:
        split_line = line.strip().split('\t')
        run_accession = split_line[5]
        if run_accession in samples_to_keep:
            out.write(line)
