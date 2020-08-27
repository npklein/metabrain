
import argparse
import os
import gzip

parser = argparse.ArgumentParser(description='Check if genes in first column have a version number and if so, remove')
parser.add_argument('matrix', 
                    help='Matrix to remove versions of genes (inplace)')

args = parser.parse_args()

if not os.path.exists(args.matrix):
    if not os.path.exists(args.matrix+'.gz'):
        raise RuntimeError(args.matrix+' and '+args.matrix+'.gz do not exist')
    args.matrix = args.matrix+'.gz'

if args.matrix.endswith('.gz'):
    f = gzip.open(args.matrix,'rt')
else:
    f = open(args.matrix)

header = f.readline()
line1 = f.readline().split('\t')
if '.' in line1[0]:
   print('genes includes version number, remove these')
else:
    print('genes do not have version number, nothing to change')
    exit(0)

with open(args.matrix+'.tmp','w') as out:
    out.write(header)
    line1[0] = line1[0].split('.')[0]
    out.write('\t'.join(line1))
    for line in f:
        line = line.split('\t')
        line[0] = line[0].split('.')[0]
        out.write('\t'.join(line))
f.close()
os.rename(args.matrix+'.tmp', args.matrix)
