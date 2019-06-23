import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description='Merge (gzipped) matrices by column')
parser.add_argument('outfile', help='File to write output to')
parser.add_argument('input_files', help='files to merge', nargs='+')
parser.add_argument('--rowname_column', help='col number of column to merge on', default = 0)

args = parser.parse_args()

rowname_column = int(args.rowname_column)
def openfile(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)


value_per_column_per_row = {}
row_set = set([])
columns = []
cols_before_merge = {}
for f in args.input_files:
    print(f)
    column_index = {}
    with openfile(f,'rt') as input_file:
        header = input_file.readline().rstrip('\n').split('\t')
        print('merging on',header[rowname_column])
        for index, element in enumerate(header[rowname_column+1:]):
            column_index[index] = element
            columns.append(element)
        for line in input_file:
            line = line.rstrip('\n').split('\t')
            rowname = line[rowname_column]
            if rowname_column != 0:
                for i in range(0,rowname_column):
                    rowname = line[i]+'_'+rowname
            row_set.add(rowname)
            if rowname not in value_per_column_per_row:
                value_per_column_per_row[rowname] = {}
            for index, element in enumerate(line[rowname_column+1:]):
                value_per_column_per_row[rowname][column_index[index]] = element
sorted_rownames = sorted(row_set)
print('done reading, start writing')
with openfile(args.outfile,'wt') as out:
    out.write('-')
    for colname in columns:
        out.write('\t'+colname)
    out.write('\n')
    for rowname in sorted_rownames:
        out.write(rowname)
        for colname in columns:
            if colname in value_per_column_per_row[rowname]:
                out.write('\t'+value_per_column_per_row[rowname][colname])
            else:
                out.write('')
        out.write('\n')


