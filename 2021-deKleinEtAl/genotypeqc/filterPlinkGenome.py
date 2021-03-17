import sys

input = sys.argv[1]
output = sys.argv[2]


print("in: "+input)
print("out: "+output)

fh = open(input,'r')
fho = open(output, 'w')
lnctr = 0
written = 0
for line in fh:
	line = line.strip()
	while "  " in line:
		line = line.replace("  "," ")
	elems = line.split()
	if lnctr == 0:
		fho.write(line+"\n")
	else:
		pihat = float(elems[9])
		if pihat >= 0.125:
			fho.write(line+"\n")
			written = written + 1
	lnctr = lnctr + 1
	if lnctr % 1000 == 0:
		print(str(lnctr)+ " lines parsed, "+str(written)+" lines written")
fh.close()
fho.close()
