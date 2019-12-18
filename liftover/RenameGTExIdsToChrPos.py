import gzip
import sys

## converts GTEx 1_10497_C_T_b37 ids to chr:pos

infile = sys.argv[1]  # SNPs.txt.gz with GTEx ids
outfile = sys.argv[2] # Output filename prefix

fh = gzip.open(infile,'rt')
fho = gzip.open(outfile+".txt.gz",'wt')
fho2 = gzip.open(outfile+"-mappings.txt.gz", 'wt')

set = set()

for line in fh:
	# dostuff
	# 1_10497_C_T_b37
	elems = line.split("_")
	outid = elems[0]+":"+elems[1]
	if outid not in set:
		output = elems[0]+":"+elems[1]
		fho.write(outid+"\n")
		fho2.write(elems[0]+"\t"+elems[1]+"\t"+outid+"\n")
		set.add(outid)
	else :
		add = 2
		outid2 = outid+"_"+str(add)
		while outid2 in set:
			outid2 = outid+"_"+str(add)
			add = add + 1
		
		fho.write(outid2+"\n")
		fho2.write(elems[0]+"\t"+elems[1]+"\t"+outid2+"\n")
		set.add(outid2)


fh.close()
fho.close()
fho2.close()
