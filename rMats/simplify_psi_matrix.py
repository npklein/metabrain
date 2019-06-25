import glob
import argparse
import gzip
import os
import multiprocessing

parser = argparse.ArgumentParser(description='Merge (gzipped) matrices by column')
parser.add_argument('psi_dir', help='Directory with psi files')
parser.add_argument('outdir', help='outdir to write to')
parser.add_argument('--multiprocess', help='if multiprocessing should be used', action='store_true')
    
args = parser.parse_args()



if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

def openfile(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def simplify_matrix(f):

    none_samples = set(['rowname','gene','splice_event','short','longExon','flanking','downstream','upstream','exon','riExon','downstread','1stExon','2ndExon'])
    samples = []
    identifiers = []
    psi_per_exon = {}
    print(f)
    with openfile(f,'rt') as input_file:
        header = input_file.readline().strip('\n').split('\t')
        for element in header:
            if element not in none_samples and element not in samples:
               samples.append(element)

        x = 0        
        for line in input_file:
            x += 1
            if x % 1000 == 0:
                print(f,'-',x,'lines processed')
            line = line.strip('\n').split('\t')
            if 'MXE' in f:
                identifier = line[3]+'_'+line[4]
            else:
                identifier = line[3]
            if identifier not in psi_per_exon:
                psi_per_exon[identifier] = {}
                identifiers.append(identifier)
            for index, element in enumerate(line):
                if header[index] in none_samples:
                    continue
                else:
                    if header[index] not in psi_per_exon[identifier]:
                        psi_per_exon[identifier][header[index]] = set([])
                    if element != '':
                       psi_per_exon[identifier][header[index]].add(element)

    print(f,' - done reading, start writing')
    print(f, ' - ',len(samples),'samples')
    print(f, ' - ',len(identifiers),'identifiers')
    with openfile(args.outdir+'/'+f.split('/')[-1].replace('.txt','.simplified.txt'),'wt') as out:
        out.write('exon')
        for sample in samples:
            out.write('\t'+sample)
        out.write('\n')
        for index, identifier in enumerate(identifiers):
            if index % 1000 == 0:
                print(f,'- written ',index)
            out.write(identifier)
            for sample in samples:
                if len(psi_per_exon[identifier][sample]) > 0:
                    psi = 0
                    for value in psi_per_exon[identifier][sample]:
                        psi += float(value)
                    psi /= len(psi_per_exon[identifier][sample])
                    out.write('\t'+str(psi))
                else:
                    out.write('\t')
            out.write('\n')
    return

input_files = glob.glob(args.psi_dir+'/*gz')
input_files_to_run = []
for f in input_files:
    outfile = args.outdir+'/'+f.split('/')[-1].replace('.txt','.simplified.txt')
    if not os.path.isfile(outfile):
        input_files_to_run.append(f)

if args.multiprocess:
    pool = multiprocessing.Pool(len(input_files_to_run))
    pool.map(simplify_matrix, input_files_to_run)
    pool.close()
else:
    for f in input_files_to_run:
        simplify_matrix(f)
