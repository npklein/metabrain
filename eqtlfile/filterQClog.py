import sys
import gzip

if len(sys.argv) < 5:
	print("Usage: chr sta sto out")
else:

	chr = int(sys.argv[1])
	sta = int(sys.argv[2])
	sto = int(sys.argv[3])
	outp = sys.argv[4]

	fh = gzip.open('./SNPQCLog.txt.gz','rt')
	fho = gzip.open(outp,'wt')

	ctr = 0
	written = 0
	for line in fh:
		line = line.strip()
		if ctr < 2:
			fho.write(line+"\n")
		else:
			elems = line.split()
			snpelems = elems[0].split(":")
			chrs = int(snpelems[0])
			pos = int(snpelems[1])
			if chrs == chr and pos > sta and pos < sto:
				fho.write(line+"\n")
				written = written + 1
		ctr = ctr + 1
		if ctr % 10000 == 0:
			print(str(ctr)+" lines read, "+ str(written) +" lines written")
	fh.close()
	fho.close()
