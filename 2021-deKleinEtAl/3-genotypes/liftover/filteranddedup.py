import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
fileout = sys.argv[3]


fh = open(file1, "r")

snps = {}

for line in fh:
	elems = line.split()
	coord = elems[0] + "_"+elems[3]
	snps[elems[1]] = coord
#1       rs12407459      0       39462571        G       A

fh.close()

fh = open(file2,"r")
fho = open(fileout,'w')

for line in fh:
	elems = line.split()
	snp = elems[1]
	if snp in snps:
		coord = elems[0] + "_"+elems[3]
		if coord == snps[snp]:
			fho.write(snp+"\n")

fh.close()	
fho.close()
