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
samples_added = set([])
x = 0
with open('ena.txt') as input_file, open('samplesheet_ENA_'+today+'.txt','w') as out:
    out.write(input_file.readline())
    for line in input_file:
        split_line = line.strip().split('\t')
        run_accession = split_line[5]
        if run_accession in samples_to_keep:
            out.write(line)
            samples_added.add(run_accession)
            x += 1

print('Added',str(x),'of',str(len(samples_to_keep)),'samples')
samples_missed = samples_to_keep-samples_added

if len(samples_missed) > 0:
    print('Missed samples:')
    print('\n'.join(samples_missed))
    with open('missing_samples_'+today+'.txt','w') as out:
        out.write('# These samples are in the sample list, but not in the ENA samplesheet. Possibly removed\n')
        for sample in samples_missed:
            out.write(sample+'\n')
