

import gzip
import sys
import os

fh = gzip.open(sys.argv[1],'rt')
fho = gzip.open(sys.argv[2],'wt')


print("in: "+sys.argv[1])
print("out: "+sys.argv[2])

for line in fh:
	line = line.replace('chr','')
	fho.write(line)

fh.close()
fho.close()

os.rename(sys.argv[2],sys.argv[1])
