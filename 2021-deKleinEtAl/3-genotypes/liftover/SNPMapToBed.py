import gzip
import sys

## converts SNPMappings.txt.gz to bed

infile = sys.argv[1]  # SNPMappings.txt.gz file
outfile = sys.argv[2] # bed out

fh = gzip.open(infile,'rt')
fho = open(outfile,'w')

for line in fh:
	# dostuff
	# 1_10497_C_T_b37
	line = line.rstrip()
	elems = line.split("\t")
	fho.write("chr"+elems[0]+"\t"+elems[1]+"\t"+str(int(elems[1])+1)+"\t"+elems[2]+"\n")


fh.close()
fho.close()
