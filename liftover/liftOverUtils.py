import argparse
import gzip
import tabix
import math

def main():
	
	parser = argparse.ArgumentParser(description='LiftOver utils.')
	parser.add_argument('--mode', choices=['gwastobed', 'liftovergwas'], help="Type of program to run", required=True)
	
	
	parser.add_argument('--in',help="Input file", dest="input")
	parser.add_argument('--out',help="Output file", dest="output")
	parser.add_argument('--bed',help="Bed file", dest="bed")
	parser.add_argument('--snpmap',help="SNPMap file", dest="snpmap")
	
	args = parser.parse_args()
	
	if args.mode == "gwastobed":
		gwastobed(args.input,args.bed,args.snpmap)
	elif args.mode == "liftovergwas":
		liftovergwas(args.input,args.bed,args.output)
		
	
# modes:
# eqtltobed
# liftovereqtl
# updatetrityperrsid
# updatetrityperpos
# gwastobed
# liftovergwas

# get file handle for reading
def getFH(infile):
	fh = None
	if infile.endswith("gz"):
		fh = gzip.open(infile,'rt')
	else:
		fh = open(infile,'r')
	return fh

# get file handle for writing
def getFHO(outfile):
	fho = None
	if outfile.endswith("gz"):
		fho = gzip.open(outfile,'wt')
	else:
		fho = open(outfile,'w')
	return fho

def readLiftOverBed(liftoverbed):
	snpmap = {}
	bedfh = getFH(liftoverbed)
	for line in bedfh:
		# chr pos1 pos2 id
		elems = line.strip().split("\t")
		chr = elems[0].replace("chr","")
		snpmap[elems[3]] = [chr,elems[1],None]
	bedfh.close()
	return snpmap

def gwastobed(infile,bedout,snpmapref):
	# not quite sure why this would need a snpmap reference
	# maybe for GWAS files where the snp positions are not given?
	# snpmap = None;
	# if snpmapref is not None:
		# snpmap = {}
		# snpset = set()
		# fh = getFH(infile)
		# fh.readline()
		# for line in fh:
			# elems=line.split()
			# snpset.add(elems[0])
		# fh.close()
		# fh = getFH(snpmapref)
		# fh.readline()
		# for line in fh:
			# elems = line.strip().split()
			# id = elems[2]
			# query = elems[0]+":"+elems[1]
			# if query in snpset:
				# snpmap[query] = id
	print("GWAS to BED: "+infile+" --> "+bedout)
	fh = getFH(infile)
	fho = getFHO(bedout)
	fh.readline()
	for line in fh:
		elems = line.strip().split()
		chr = elems[6].lower()
		bp = int(elems[7])
		if not chr.startswith("chr"):
			chr = "chr"+chr
		fho.write(chr+"\t"+str(bp)+"\t"+str(bp+1)+"\t"+elems[0]+"\n")
	fh.close()
	fho.close()
	print("GWAS to BED: "+infile+" --> "+bedout+" DONE.")

def liftovergwas(infile,liftoverbed,outfile):
	# read in lifted snps
	print("LiftOver GWAS: "+infile+" --> "+outfile+" using "+liftoverbed)
	snpmap = readLiftOverBed(liftoverbed)
	fh = getFH(infile)
	fho = getFHO(outfile)
	header = fh.readline()
	header = header.strip().split()
	fho.write("\t".join(header)+"\tIdBeforeLiftOver\n")
	for line in fh:
		elems = line.strip().split()
		query = elems[0]
		update = snpmap.get(query)
		if update is not None:
			elems[0] = update[0]+":"+update[1]
			# SNP     A1      A2      BETA    SE      P       CHR     BP
			elems[6] = update[0]
			elems[7] = update[1]
			fho.write("\t".join(elems)+"\t"+query+"\n")
	fh.close()
	fho.close()
	print("LiftOver GWAS: "+infile+" --> "+outfile+" using "+liftoverbed+" DONE.")


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
	snpmap = readLiftOverBed(liftoverbed)
	
	# read new rsids
	dbsnp = tabix.open(dbsnpvcf)
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