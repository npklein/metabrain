import gzip
import sys
import os

inputfile = sys.argv[1]
tmpfile = inputfile+"_tmp.txt.gz"

print("Stripping non-rs IDs from: "+inputfile)

fh = gzip.open(inputfile, 'rt')
fho = gzip.open(tmpfile, 'wt')

rsvisit = set()

norsctr = 0

for line in fh:
	elems = line.split(":")
	if len(elems) < 3:
		print("Error parsing line:")
		print(line)
		sys.exit()
	else:
		rs = elems[2]

		ctr = 2
		if rs == 'nors':
			norsctr = norsctr + 1
			fho.write("nors_"+str(norsctr)+"\n")
		elif rs in rsvisit:
			tmprs = rs+"_"+str(ctr)
			while tmprs in rsvisit:
				tmprs = rs+"_"+str(ctr)
				ctr = ctr + 1
			fho.write(tmprs+"\n")
			rsvisit.add(tmprs)
		else:
			fho.write(rs+"\n")
			rsvisit.add(rs)


fh.close()
fho.close()


os.rename(tmpfile,inputfile)

