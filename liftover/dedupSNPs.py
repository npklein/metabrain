import gzip
import sys

## replaces duplicate IDs in SNPs.txt.gz with id_number

infile = sys.argv[1]  # SNPs.txt.gz with ids
outfile = sys.argv[2] # Output filename prefix

fh = gzip.open(infile,'rt')
fho = gzip.open(outfile+".txt.gz",'wt')

set = set()

for line in fh:
	# dostuff
	# 1_10497_C_T_b37
	outid = line.rtrim()
	if outid not in set:
		fho.write(outid+"\n")
		set.add(outid)
	else:
		add = 2
		outid2 = outid+"_"+str(add)
		while outid2 in set:
			outid2 = outid+"_"+str(add)
			add = add + 1
		
		fho.write(outid2+"\n")
		set.add(outid2)

fh.close()
fho.close()
