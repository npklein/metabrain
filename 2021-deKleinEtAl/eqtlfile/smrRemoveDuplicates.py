

import sys
import gzip

input = sys.argv[1]
output = sys.argv[2]

fh = gzip.open(input, 'rt')
fho = gzip.open(output,'wt')
fho.write(fh.readline())

ids = {}
eqtls = set()

print("Parsing: "+input)
print("Output: "+output)
written = 0
ctr = 0
for line in fh:
	elems = line.strip().split("\t")
	rsid = elems[0]
	a1 = elems[3]
	a2 = elems[4]
	alleles = ids.get(rsid)
	include = 0
	if alleles is None:
		alset = set()
		alset.add(a1)
		alset.add(a2)
		ids[rsid] = [alset,a1]
		include = 1
	else:
		alset = alleles[0]
		refa1 = alleles[1]
		if a1 in alset and a2 in alset and a1 == refa1:
			include = 1
	eqtl = rsid + "_"+ elems[6]
	if include == 1 and eqtl not in eqtls:
		fho.write("\t".join(elems)+"\n")	
		written = written + 1
		eqtls.add(eqtl)
	ctr = ctr + 1
	if ctr % 1000000 == 0:
		print("{} lines read, {} lines written".format(ctr,written))
fho.close()
fh.close()
print("Done. {} lines read, {} lines written".format(ctr,written))
