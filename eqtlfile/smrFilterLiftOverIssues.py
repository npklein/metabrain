import gzip
import sys

input = sys.argv[1]
output = sys.argv[2]
chr = sys.argv[3]

print(input)
print(output)
print(chr)

fh = gzip.open(input, 'rt')
fho = gzip.open(output,'wt')
fho.write(fh.readline())
for line in fh:
	elems = line.split("\t")
	snp  = elems[0]
	snpelems = snp.split(":")
	if snpelems[0] == chr:
		fho.write(line)
fh.close()
fho.close()
