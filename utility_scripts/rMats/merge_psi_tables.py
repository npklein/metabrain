# merge psi (percent spliced-in)  tables together
# Round incl. level to 4 decimal behind ocmma
import glob
import argparse
import os

parser = argparse.ArgumentParser(description='Generate rMATs files.')
parser.add_argument('psi_base_dir', help='Base dir which contains all subdirectories from which files need to be merged')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('file_prefix', help='Prefix for the file names. Will create files in output_directory with <prefix>.<event>.txt')


args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)


events = ['A3SS','A5SS','MXE','RI','SE']


incl_per_sample_per_event = {'A3SS':{},'A5SS':{},'MXE':{},'RI':{},'SE':{}}
exons_per_event = {}
samples = set([])
genes = {}
samples = []
#ID    GeneID    geneSymbol    chr    strand    longExonStart_0base    longExonEnd    shortES    shortEE    flankingES    flankingEE    ID    IJC_NEUNE940RZC_CGND_HRA_00399    SJC_NEUNE940RZC_CGND_HRA_00399    IJC__NEUVB305VU7_CGND-HRA-01231    SJC__NEUVB305VU7_CGND-HRA-01231    IncFormLen    SkipFormLen    IncLevel1    IncLevel2    IncLevelDifference
print('read data')
for event in events:
    print(event)
    if event not in exons_per_event:
        exons_per_event[event] = set([])
    event_files = glob.glob(args.psi_base_dir+'/**/'+event+'*txt', recursive=True)
    n_event_files = len(event_files)
    print('parsing through',n_event_files,'files')
    for index, f in enumerate(event_files):
        if index % 100 == 0:
            print(str(index)+'/'+str(n_event_files))
        sample = f.split('/')[-2]
        samples.append(sample)
        with open(f) as input_file:
            header = input_file.readline()
            for line in input_file:
                line = line.strip().split('\t')
                gene = line[0]
                exon = line[4]+'_'+line[5]+'_'+line[6]+'_'+line[7]+'_'+line[8]
                genes[exon] = gene
                exons_per_event[event].add(exon)
                incl_level = line[-1]
                if incl_level == 'NA':
                    incl_level = ''
                else:
                    incl_level = str(round(float(incl_level),4))
                if exon not in incl_per_sample_per_event[event]:
                    incl_per_sample_per_event[event][exon] = {}
                incl_per_sample_per_event[event][exon][sample] = incl_level


for event in events:
    exons_sorted = sorted(exons_per_event[event])
    outfile = args.output_directory+'/'+args.file_prefix+'.'+event+'.txt'
    with open(outfile, 'w') as out:
        out.write('gene\texon')
        for sample in samples:
            out.write('\t'+sample)
        out.write('\n')

        for exon in exons_sorted:
            out.write(genes[exon]+'\t'+exon)
            for sample in samples:
                try:
                    out.write('\t'+incl_per_sample_per_event[event][exon][sample])
                except:
                    out.write('\t')
            out.write('\n')
    print('file written to',outfile)    
