from optparse import OptionParser
import gzip
import tabix
import math

def main():
    print("python main function")

# modes:
# eqtltobed
# liftovereqtl
# updatetrityperrsid
# updatetrityperpos
#
#

def getFH(infile):
	fh = None
	if infile.endswith("gz"):
		fh = gzip.open(infile,'rt')
	else:
		fh = open(infile,'r')
	return fh

def getFHO(outfile):
	fho = None
	if outfile.endswith("gz"):
		fho = gzip.open(outfile,'wt')
	else:
		fho = open(outfile,'w')
	return fho

def snpmaptobed(infile,bedout):
	fh = getFH(infile)
	fho = getFHO(bedout)
	for line in fh:
		line = line.rstrip()
		elems = line.split("\t")
		fho.write("chr"+elems[0]+"\t"+elems[1]+"\t"+str(int(elems[1])+1)+"\t"+elems[2]+"\n")
	fh.close()
	fho.close()

def eqtlfiletobed(infile,bedout):
	fh    = getFH(infile)
	fho   = getFHO(bedout)
	idset = set()
	fh.readline()
	for line in fh:		
		elems = line.split("\t")
		id = elems[1]
		if id not in idset:
			fho.write("chr"+elems[2]+"\t"+elems[3]+"\t"+str(int(elems[3])+1)+"\t"+elems[1]+"\n")
			idset.add(id)
	fh.close()
	fho.close()
	
def updatetrityperpos(originalbed, liftoverbed, outdir):
	print("not implemented yet")

def updatetrityperrsid(liftoverbed, dbsnpvcf, outdir):
	print("not implemented yet")
	
def liftovereqtl(eqtlfilein, eqtlfileout, geneannotation, liftoverbed, dbsnpvcf):
	
	# read in lifted snps
	snpmap = {}
	bedfh = getFH(liftoverbed)
	for line in bedfh:
		# chr pos1 pos2 id
		elems = line.strip().split("\t")
		chr = elems[0].replace("chr","")
		snpmap[elems[3]] = [chr,elems[1],None]
	bedfh.close()
	
	# read new rsids
	dbsnp = tabix.open(dbsnpvcf
	for k,v in snpmap.items():
		pos = int(v[1])
		records = dbsnp.query(v[0],pos,pos+1)
		for record in records:
			# ideally, we want to compare alleles here...
			
			# skip indels
			ref=record[3]
			alt=record[4].split(",")
			if len(ref) == 1 and len(alt[0]) == 1:
				v[2]=record[2]
	dbsnp.close()
	
	# read gene annotation
	genemap = {}
	genefh=getFH(geneannotation)
	for line in genefh:
		elems = line.strip().split()
		sta = int(elems[2])
		sto = int(elems[3])
		diff = math.floor( (sto-sta)/2 )
		midpoint = sta + diff
		genemap[elems[0]] = [elems[1],midpoint]		
	genefh.close()
	
	# iterate original eqtl file
	# write new version at the same time
	eqtlinfh = getFH(eqtlfilein)
	eqtloutfh = getFHO(eqtlfileout)
	
	eqtloutfh.write(eqtlinfh.readline())
	for line in eqtlinfh:
		elems = line.split("\t")
		
		# update snp
		id = elems[1]
		vals = snpmap.get(id)
		if vals is not None:
			elems[1] = vals[2]
			elems[2] = vals[0]
			elems[3] = vals[1]
			
			# update gene
			# split gene by dot
			geneid = elems[4].split(".")[0]
			genevals = genemap.get(geneid)
			if genevals is not None:
				elems[4] = geneid
				elems[5] = genevals[0]
				elems[6] = genevals[1]
				eqtloutfh.write("\t".join(elems))
		
	eqtlinfh.close()
	eqtloutfh.close()


if __name__ == '__main__':
    main()