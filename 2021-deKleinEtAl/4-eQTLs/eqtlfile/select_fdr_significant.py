import gzip
import glob
import sys

f=sys.argv[1]
outfile=sys.argv[2]

with gzip.open(f,'rt') as input_file, gzip.open(outfile,'wt') as out:
        x = 0
        out.write(input_file.readline())
        for line in input_file:
            x += 1
            if x % 1000000 == 0:
                print(x)
            fdr = float(line.strip().split('\t')[-1])
            if fdr < 0.05:
                out.write(line)

