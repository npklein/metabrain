import sys
import gzip 
import os

inputfile=sys.argv[1]
referencefile=sys.argv[2]

print("Using "+referencefile+" to update "+inputfile)

idstoload = set()
snpmap = {}
fh = open(inputfile,'r')
ctr=0
for line in fh:
	elems = line.rstrip().split('\t')
	id = elems[0]+":"+elems[3]
	idstoload.add(id)
	ctr = ctr + 1
	if (ctr % 1000000 == 0):
		print(str(ctr)+ " lines read")
fh.close()

print("done reading 1")

# load snpmap
fh2 = gzip.open(referencefile,'rt')
ctr = 0
for line in fh2:
	elems = line.rstrip().split('\t')
	id = elems[0]+":"+elems[1]
	if id in idstoload:
		snpmap[id]=elems[2]
	ctr = ctr + 1
	if (ctr % 1000000 == 0):
		print(str(ctr)+ " lines read, "+str(len(snpmap)) +" found")

fh2.close()

print("done loading 2")

fho = open(inputfile+'.upd','w')
fh = open(inputfile,'r')
ctr = 0
for line in fh:
	# 1       chr1:752566     0       752566  G       A
	elems = line.rstrip().split('\t')
	id = elems[0]+":"+elems[3]
	out = elems[0]+"\t"+elems[1]+"\t"+elems[2]+"\t"+elems[3]+"\t"+elems[4]+"\t"+elems[5]
	try:
		rs = snpmap.get(id, None)
		if rs is not None:
			rs = snpmap[id]
			out = elems[0]+"\t"+rs+"\t"+elems[2]+"\t"+elems[3]+"\t"+elems[4]+"\t"+elems[5]
	except KeyError:
		print("Could not find: "+id)
	fho.write(out+"\n")
	ctr = ctr + 1
	if (ctr % 1000000 == 0):
		print(str(ctr)+ " lines read")


fh.close()
fho.close()
	

os.rename(inputfile+'.upd',inputfile)
