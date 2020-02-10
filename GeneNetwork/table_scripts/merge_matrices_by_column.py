import glob
import re
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser(description='Merge matrices by columns (no check is done on rownames, so rows have to be in same order!).')
parser.add_argument('-o','--outfile_prefix', help='Prefix that will be added to the output file name', required=True)
parser.add_argument('-r', '--rootdir', help='Root directory that contins the prediction files', required = True)
parser.add_argument('-t','--type', help='List of types, has to be same as the subdirectory in root. E.g. --type go_C go_F go_P',
                    nargs='+', required=True, type=list)

args = parser.parse_args()
print(args.accumulate(args.integers))

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in _nsre.split(s)]

def merge(type):
    print(type)
    print("search in "+args.rootdir)
    column_header = []
    gene_list = []
    gene_set = set([])
    gene_values = {}
    for f in sorted(glob.glob(args.rootdir+'/'+type+'/*txt'), key=natural_sort_key):
        with open(f) as input_file:
            header = input_file.readline().strip().split('\t')
            for element in header[1:]:
                column_header.append(element)
            for line in input_file:
                line = line.strip().split('\t')
                if line[0] not in gene_set:
                    gene_list.append(line[0])
                    gene_set.add(line[0])
                if line[0] not in gene_values:
                    gene_values[line[0]] = []
                gene_values[line[0]].extend(line[1:])

    with open(args.rootdir+'/'args.outfile_prefix+type+'_predictions.txt','w') as out:
        out.write('-\t'+'\t'.join(column_header)+'\n')
        for gene in gene_list:
            out.write(gene+'\t'+'\t'.join(gene_values[gene])+'\n')


p = Pool(len(args.type))
p.map(merge, args.type)
p.close()
p.join()
