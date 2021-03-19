import sys
import gzip


if len(sys.argv) < 4:
	ctr = 0
	for arg in sys.argv:
		print("{}\t{}".format(ctr,arg))
		ctr = ctr + 1
	print("Usage: input.txt.gz output.txt.gz filter.txt")
	exit()



input=sys.argv[1]
output=sys.argv[2]
filter=sys.argv[3]

print("Input "+input)
print("Output "+output)
print("Filter "+filter)

filterset=set()
fh = open(filter,'r')
for line in fh:
	elems = line.strip().split("\t")
	if len(elems) > 1:
		filterset.add(elems[0]+"-"+elems[1])
fh.close()


fh = gzip.open(input,'rt')
fho = gzip.open(output,'wt')
fho.write(fh.readline())

for line in fh:
	elems = line.strip().split("\t")
	eqtl = elems[1]+"-"+elems[4]
	if eqtl in filterset:
		fho.write(line)
fho.close()
fh.close()
