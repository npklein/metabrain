import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description='Merge (gzipped) matrices by column')
parser.add_argument('output_file', help='File to write output to')
parser.add_argument('input_files', help='files to merge', nargs='+')

args = parser.parse_args()

def openfile(filename, mode='rt',gzip=True):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)


value_per_column_per_row = {}
row_set = set([])
columns = []
for f in ars.input_files:
    print(f)
    column_index = {}
    with openfile(f,'rt') as input_file:
        header = input_file.readline().rstrip('\n').split('\t')
        for index, element in enumerate(header[1:]):
            column_index[index] = element
            columns.append(element)
        for line in input_file:
            line = line.rstrip('\n').split('\t')
            rowname = line[0]
            row_set.add(rowname)
            if rowname not in value_per_column_per_row:
                value_per_column_per_row[rowname] = {}
            for index, element in enumerate(line[1:]):
                value_per_column_per_row[rowname][column_index[index]] = element

sorted_rownames = sorted(row_set)
print('done reading, start writing')
with openfile(outfile,'wt') as out:
    out.write('-')
    for sample in columns:
        out.write('\t'+sample)
    out.write('\n')
    for rowname in sorted_rownames:
        out.write(rowname)
        for sample in columns:
            out.write('\t'+value_per_column_per_row[rowname][sample])
        out.write('\n')


