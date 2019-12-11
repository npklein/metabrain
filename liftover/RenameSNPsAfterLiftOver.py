import gzip
import sys

## rewrites lifted RSIds and SNPMap using a reference set
infile = sys.argv[1]   # SNPMappings.txt.gz after liftover
infile2 = sys.argv[2]  # Reference SNPMappings.txt.gz (e.g. from dbSNP)
infile3 = sys.argv[3]  # SNPs.txt.gz after liftover
outfile = sys.argv[4]  # output file prefix



chrpostoread = set()
oldrsidmap = {}

# read lifted SNPMappings file
fh1 = gzip.open(infile, 'rt')
lncounter = 0
for line in fh1:
	line = line.rstrip()
	elems = line.split("\t")
	rs = elems[2]
	try:	
		chr = int(elems[0])
		chrid = elems[0]+"\t"+elems[1]
		chrpostoread.add(chrid)		
		oldrsidmap[rs] = chrid	
	except ValueError:
		oldrsidmap[rs] = "0\t0"
	lncounter = lncounter + 1
	if (lncounter % 1000000 == 0):
		print(str(lncounter)+ " lines parsed")
				
fh1.close()
print("done parsing: "+infile)

# read new rsIDS from reference
chrpostorsmap = {}
rsidset = set()
fh2 = gzip.open(infile2, 'rt')
lncounter = 0
for line in fh2:
	line = line.rstrip()
	elems = line.split("\t")
	rs = elems[2]
	chrid = elems[0]+"\t"+elems[1]
	if chrid in chrpostoread:
		if rs not in rsidset:
			chrpostorsmap[chrid] = rs
			rsidset.add(rs)
	lncounter = lncounter + 1
        if (lncounter % 1000000 == 0):
                print(str(lncounter)+ " lines parsed, "+str(len(rsidset)) +" loaded out of "+str(len(chrpostoread)))

fh2.close()

print("Done parsing: "+infile2)


# process SNPs.txt.gz
fh3 = gzip.open(infile3, 'rt')
fho1 = gzip.open(outfile+".txt.gz",'wt')
fho2 = gzip.open(outfile+"-mappings.txt.gz", 'wt')
for line in fh3:
	line = line.rstrip()
	chrid = oldrsidmap.get(line)
	chridelems = chrid.split("\t")
	if chridelems[0] is not "0":
		newrs = chrpostorsmap.get(chrid)
		if newrs is not None:
			fho1.write(newrs+"\n")
			fho2.write(chrid+"\t"+newrs+"\n")
		else:
			# position not found in reference; think of something else	
			fho1.write(chridelems[0]+":"+chridelems[1]+"\n")
			fho2.write("0\t0\t"+chridelems[0]+":"+chridelems[1]+"\n")
	else:
		# rs id not lifted (chr == 0); write something else.
		fho1.write(line+"\n")
		fho2.write(chrid+"\t"+line+"\n")
		

fh3.close()
fho1.close()
fho2.close()

		
