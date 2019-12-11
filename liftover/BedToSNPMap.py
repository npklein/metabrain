import gzip
import sys

## writes out lifted positions for original list of SNPs

infile = sys.argv[1]  # Lifted BED file
infile2 = sys.argv[2] # SNPs.txt.gz to be updated
outfile = sys.argv[3] # outfile prefix

fh = open(infile,'r')

liftedmap = {}

# parse lifted bed file
for line in fh:
	line = line.rstrip()
	elems = line.split("\t")
	rs = elems[3]
	chr = elems[0]
	chr = chr.replace("chr","")
	pos = elems[1]
	liftedmap[rs] = chr+"\t"+pos
fh.close()

# parse SNPs.txt.gz
fh2 = gzip.open(infile2, 'rt')
fho1 = gzip.open(outfile+".txt.gz", 'wt')
fho2 = gzip.open(outfile+"-mappings.txt.gz", 'wt')

for line in fh2:
	line = line.rstrip()
	newpos = liftedmap.get(line)
	if newpos is not None:
		fho1.write(line+"\n")
		fho2.write(newpos+"\t"+line+"\n")
	else:
		# map to 0	
		fho1.write(line+"_unlifted\n")
		fho2.write("0\t0\t"+line+"_unlifted\n")

fh2.close()
fho1.close()
fho2.close()
