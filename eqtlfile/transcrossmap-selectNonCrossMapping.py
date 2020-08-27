import gzip
import sys


if len(sys.argv) < 4:
	print("Usage: eQTLsFDR0.05.txt.gz eQTLCrossMapFile.txt.gz percentagecrossmap outfile.txt.gz")
else:
	infile=sys.argv[1]
	crossmapfile=sys.argv[2]
	crossmapthreshold=float(sys.argv[3])
	outfile=sys.argv[4]
	
	selectedset = set()
	fh = gzip.open(crossmapfile, 'rt')
	fh.readline()
	for line in fh:
		elems = line.strip().split()
		snp = elems[0]
		gene = elems[1]
		perc = float(elems[4])
		if perc < crossmapthreshold:
			eqtl = snp+"-"+gene
			selectedset.add(eqtl)
	fh.close()

	fh = gzip.open(infile, 'rt')
	fho = gzip.open(outfile, 'wt')
	fho.write(fh.readline())
	for line in fh:
		elems = line.split("\t")
		snp = elems[1]
		gene = elems[4]
		eqtl = snp+"-"+gene
		if eqtl in selectedset:
			fho.write(line)
	fh.close()
	fho.close()
