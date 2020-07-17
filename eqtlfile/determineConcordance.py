import gzip
import sys

if len(sys.argv) < 3:
	print("Usage: eqtls1.txt.gz eqtls2.txt.gz output.txt.gz [name1 name2]")
	exit()

file1=sys.argv[1]
file2=sys.argv[2]
output=sys.argv[3]

name1="Dataset1"
name2="Dataset2"

if len(sys.argv) > 5:
	name1 = sys.argv[4]
	name2 = sys.argv[5]

def countSharedAlleles(al1, al2):
	# count number of shared alleles
	shared = 0
	i = 0
	while i < 2:
		j = 0
		while j < 2:
			if al1[i] == al2[j]:
				shared = shared + 1
			j = j + 1
		i = i + 1	
	return shared	

def flipalleles(alleles1, assessed1, alleles2, assessed2):
	al1 = alleles1.split("/")
	al2 = alleles2.split("/")

	shared = countSharedAlleles(al1,al2)
	if shared < 2:
		# try complement
		al2 = complement(al2)
		shared = countSharedAlleles(al1,al2)
	
	if shared == 2:
		if assessed1 == assessed2:
			return 0
		else:
			return 1
	else:
		return -1


eqtlmap={}

total = 0
shared = 0
mismatchingalleles = 0
concordant = 0


fh=gzip.open(file1, 'rt')
fh.readline()
for line in fh:
	elems = line.strip().split("\t")
	snp = elems[1]
	gene = elems[2]
	zscore = float(elems[10])
	# fdr = float(elems[len(elems)-1])
	alleles = elems[8]
	assessed = elems[9]
	eqtl = snp+"-"+gene
	eqtlmap[eqtl] = [zscore, alleles, assessed]
	total = total + 1
fh.close()

fho = gzip.open(output, 'wt')
fho.write("SNP\tGene\tAlleles1\tAssessed1\tZ1\tAlleles2\tAssessed2\tZ2\tflipped\n")
fh=gzip.open(file2, 'rt')
fh.readline()

for line in fh:
	elems = line.strip().split("\t")
	snp = elems[1]
	gene = elems[2]
	zscore = float(elems[10])
	fdr = float(elems[len(elems)-1])
	alleles = elems[8]
	assessed = elems[9]
	eqtl = snp+"-"+gene
	ref = eqtlmap.get(eqtl)
	if ref is not None:
		flip = flipalleles(ref[1],ref[2],alleles,assessed)
		if flip == 1:
			zscore = zscore * -1
		if flip < 0:
			mismatchingalleles = mismatchingalleles + 1
		if flip > -1:
			shared = shared + 1
#			print(snp)
#			print(gene)
#			print(ref[0])
#			print(ref[1])
#			print(ref[2])
#			print(alleles)
#			print(assessed)
#			print(zscore)
#			print(flip)
			fho.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(snp,gene,ref[1],ref[2],ref[0],alleles,assessed,zscore,flip))
			if (ref[0] >= 0 and zscore >=0) or (ref[0] < 0 and zscore < 0):
				concordant = concordant + 1
fh.close()
fho.close()

fho = open(output+"-Summary.txt",'w')
fho.write("Dataset1\tDataset2\tEqtlsInDataset1\tShared\tConcordant\tConcordantOverTotal\tConcordantOverShared\tMismatchingAlleles\n")
percovertotal = concordant / total
percovershared = concordant / shared

outln = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(name1,name2,total,shared,concordant,percovertotal,percovershared,mismatchingalleles)
fho.write(outln+"\n")
print(outln)
fho.close()
