import gzip
import sys

## converts chr:pos to rsids

infile = sys.argv[1]   # SNPs.txt.gz containing chr:pos as ids 
infile2 = sys.argv[2]  # SNPMappings.txt.gz with appropriate rsids
outfile = sys.argv[3]  # output file prefix

print("reading "+infile)
fh = gzip.open(infile,'rt')
idset = set()

# read snp ids
lins = 0
for line in fh:
	# dostuff
	# expect chr:pos
	line = line.strip('\r\n')
	idset.add(line)
	lins = lins + 1
	if (lins % 1000000 ==0 ):
		print(str(lins)+" lines read")
fh.close()

print("done loading "+infile)

# read snp map
fh2 = gzip.open(infile2, 'rt')

snpmap = {}
nrlines = 0
for line in fh2:
	line = line.strip('\r\n')
	elems = line.split('\t')
	id  = elems[0]+":"+elems[1]
	if id in idset:
		snpmap[id] = elems[2]
	nrlines = nrlines + 1
	if (nrlines % 1000000) == 0:
		print(str(nrlines)+" lines read, "+str(len(snpmap))+" in dict")
		

fh2.close()
print("done reading "+infile2)



fh = gzip.open(infile,'rt')
fho = gzip.open(outfile+".txt.gz", 'wt') 
fho2 = gzip.open(outfile+"-mappings.txt.gz", 'wt')

written = set()
nrlines = 0
nrfound = 0
# read snp ids
for line in fh:
	# dostuff
	line = line.strip('\r\n')
	elems = line.split("_")[0].split(":")
	rsid = snpmap.get(line)
	
	if rsid is not None:
		nrfound = nrfound + 1
		# rs id found
		# rsid = snpmap[line]
		if rsid not in written:
			written.add(rsid)
			fho.write(rsid+"\n")
			fho2.write(elems[0]+"\t"+elems[1]+"\t"+rsid+"\n")
		else:
			add = 2
	                outid2 = rsid+"_"+str(add)
                	while outid2 in written:
                        	outid2 = rsid+"_"+str(add)
	                        add = add + 1	
			fho.write(rsid+"\n")
			fho2.write(elems[0]+"\t"+elems[1]+"\t"+outid2+"\n")
			written.add(outid2)			
	else :
		# write original id
		fho.write(line+"\n")
		fho2.write(elems[0]+"\t"+elems[1]+"\t"+line+"\n")
	nrlines = nrlines + 1
	if (nrlines % 100000) == 0:
                print(str(nrlines)+" lines read, "+str(nrfound)+" found")		
	
	
fh.close()

fho.close()
fho2.close()
