# merge psi (percent spliced-in)  tables together

import glob
events = ['A3SS','A5SS','MXE','RI','SE']


incl_per_sample_per_event = {'A3SS':{},'A5SS':{},'MXE':{},'RI':{},'SE':{}}
exons_per_event = {}
samples = set([])
genes = {}

#ID    GeneID    geneSymbol    chr    strand    longExonStart_0base    longExonEnd    shortES    shortEE    flankingES    flankingEE    ID    IJC_NEUNE940RZC_CGND_HRA_00399    SJC_NEUNE940RZC_CGND_HRA_00399    IJC__NEUVB305VU7_CGND-HRA-01231    SJC__NEUVB305VU7_CGND-HRA-01231    IncFormLen    SkipFormLen    IncLevel1    IncLevel2    IncLevelDifference
for event in events:
    if event not in exons_per_event:
        exons_per_event[event] = set([])
    for f in glob.iglob('*/'+event+'*txt'):
        print(f)
        sample1 = f.split('___')[0]
        sample2 = f.split('___')[1]
        samples.add(sample1)
        samples.add(sample2)
        with open(f) as input_file:
            header = input_file.readline()
            for line in input_file:
                line = line.strip().split('\t')
                gene = line[1]
                exon = line[3]+'_'+line[5]+'_'+line[6]
                genes[exon] = gene
                exons_per_event[event].add(exon)
                sample1_incl_level = line[-3]
                sample2_incl_level = line[-2]
                if exon not in incl_per_sample_per_event[event]:
                    incl_per_sample_per_event[event][exon] = {}
                incl_per_sample_per_event[event][exon][sample1] = sample1_incl_level
                incl_per_sample_per_event[event][exon][sample2] = sample2_incl_level

samples_sorted = sorted(samples)

for event in events:
    exons_sorted = sorted(exons_per_event[event])
    with open(event+'_merged.txt','w') as out:
        out.write('gene\texon')
        for sample in samples_sorted:
            out.write('\t'+sample)
        out.write('\n')

        for exon in exons_sorted:
            out.write(genes[exon]+'\t'+exon)
            for sample in samples_sorted:
                try:
                    out.write('\t'+incl_per_sample_per_event[event][exon][sample])
                except:
                    out.write('\tNA')
            out.write('\n')
